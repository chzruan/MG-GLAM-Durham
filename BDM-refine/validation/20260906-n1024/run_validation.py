"""Production-scale GR validation; native large allocations require Slurm.

Run with micromamba run -n cosemu python3 -B. Completed stage receipts allow
verified reuse. Failed stages retain their logs and outputs; no implicit restart
from potentially incomplete Fortran checkpoints is attempted.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import socket
import struct
import subprocess
import tarfile
import time
import traceback

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(os.environ.get('BDM_VALIDATION_ROOT', Path(__file__).resolve().parent)).resolve()
REPO = ROOT.parents[2]
WORK = ROOT / 'work'
BIN = WORK / 'bin'
AUDIT = REPO / 'BDM-refine/analysis/full-audit-20260906'
REPAIRS = REPO / 'BDM-refine/repairs/20260906'
CONFIG = 'iVirial=1\nMassMin=2.5e12\nRext=0.15\n'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while chunk := f.read(8 * 1024**2):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    staging = path.with_name(path.name + '.tmp')
    with staging.open('w') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    staging.replace(path)


def now():
    return datetime.now(timezone.utc).isoformat()


def copy_member(tar, name, destination, expected=None, executable=False):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        with tar.extractfile(name) as src, destination.open('xb') as out:
            shutil.copyfileobj(src, out, length=8 * 1024**2)
        destination.chmod(0o755 if executable else 0o644)
    if expected is not None and sha(destination) != expected:
        raise ValueError(f'Archive member hash mismatch: {destination}')
    return sha(destination)


def prepare(box):
    """Restore only nine binaries and three small inputs, verifying archives."""
    WORK.mkdir(parents=True, exist_ok=True)
    receipt = ROOT / 'preparation.json'
    build = json.loads((REPAIRS / 'native-build.json').read_text())
    for name, expected in build['source_sha256'].items():
        if sha(REPO / name) != expected:
            raise ValueError(f'Current source differs from audited native build: {name}')
    archives = {
        str(REPAIRS / 'work-artifacts.tar.gz'): '62a4e2cd1b7363c4a029d570801252c324fb871ec1739625a64242bc8fed134c',
        str(AUDIT / 'work-artifacts.tar.gz'): 'b39e39ea17aaf65a53f42a2058986717b89d8ff1513cb02fe65e14294ff9a9c9',
    }
    for path, expected in archives.items():
        if sha(path) != expected:
            raise ValueError(f'Archive hash mismatch: {path}')
    binary_hashes = {}
    with tarfile.open(REPAIRS / 'work-artifacts.tar.gz', 'r:gz') as tar:
        for name, expected in build['binaries_sha256'].items():
            binary_hashes[name] = copy_member(tar, 'work/native-build/' + name, BIN / name, expected, True)
    old = json.loads((AUDIT / 'simulation-n128.json').read_text())
    old_bdm = next(r['sha256'] for r in old['records'] if r['executable'].endswith('/PMP2BDM.exe'))
    with tarfile.open(AUDIT / 'work-artifacts.tar.gz', 'r:gz') as tar:
        binary_hashes['PMP2BDM.old.exe'] = copy_member(tar, 'work/full-build/PMP2BDM.exe',
                                                     BIN / 'PMP2BDM.old.exe', old_bdm, True)
        for name in ['Init.dat', 'PkTable.dat', 'TableSeeds.dat']:
            copy_member(tar, 'reference-inputs/' + name, WORK / 'reference-inputs' / name)
    for name, expected in {
        'PkTable.dat': 'b90b23d14a403544045ec93533402a70db8b5bc3381a3694547dfd742e0476f2',
        'TableSeeds.dat': '31e0ae21219dadc93da86c6cfe8cd1972e406f7bbf80c2a1046bafe5ab506e4d',
    }.items():
        assert sha(WORK / 'reference-inputs' / name) == expected
    spec = dict(box_mpc_h=box, gravity='GR', nrow=1024, ngrid=2048,
                particle_limit_exclusive=1200**3, source_commit='b6e96669b6f7a8cba2f3dcaa9853278aac7edf11',
                old_finder_commit='a8c7715', sources_sha256=build['source_sha256'],
                binaries_sha256=binary_hashes, archive_sha256=archives,
                cosmology=dict(Omega_m=.3089, Omega_Lambda=.6911, h=.6774, sigma8=.8159),
                mass_one_estimate_msun_h=2.775e11 * .3089 * (box / 1024)**3,
                mesh_spacing_mpc_h=box / 2048, initial_redshift=100, realization=1,
                halo_config=CONFIG, main_epochs=[2, 1, 0], pilot_epochs=[0],
                interpretation='Finder validation and same-snapshot before/after comparison; not resolution convergence')
    if receipt.exists():
        assert json.loads(receipt.read_text()) == spec, 'Refuse to change an existing experiment'
    else:
        write_json(receipt, spec)
    print('Prepared verified binaries and inputs', flush=True)


def case_spec(name):
    spec = json.loads((ROOT / 'preparation.json').read_text())
    assert name in ['pilot', 'main']
    spec.update(nrow=512 if name == 'pilot' else 1024, ngrid=1024 if name == 'pilot' else 2048,
                epochs=[0] if name == 'pilot' else [2, 1, 0])
    assert 0 < spec['nrow']**3 < spec['particle_limit_exclusive']
    return spec


def make_init(spec):
    changes = {'Box': spec['box_mpc_h'], 'Nrow': spec['nrow'], 'Ngrid': spec['ngrid'],
               '#outputs': -len(spec['epochs']), 'Steps between checkpoints': 10000,
               'Save snapshots': 1, 'DM power spectrum': 0, 'Find BDM halos': 1,
               'MG_flag': 0, 'MG_test': 0}
    source = (WORK / 'reference-inputs/Init.dat').read_text().splitlines()
    lines = []
    cursor = 0
    while cursor < len(source):
        line = source[cursor]
        key = line.split('=')[0].strip()
        if key == '#outputs':
            previous = abs(int(line.split('=')[1].split()[0]))
            lines.append(f'#outputs = {changes[key]}')
            lines.extend(f' {z:.8f}' for z in spec['epochs'])
            cursor += previous + 1
            continue
        lines.append(f'{key} = {changes[key]}' if key in changes else line)
        cursor += 1
    return '\n'.join(lines) + '\n'


def run_native(binary, cwd, stdin, threads, tag, receipts):
    """Persist logs while running; completed receipts are checked before reuse."""
    binary = BIN / binary
    expected = json.loads((ROOT / 'preparation.json').read_text())['binaries_sha256'][binary.name]
    assert sha(binary) == expected
    output = cwd / (tag + '.log')
    receipt = cwd / (tag + '.json')
    timing = cwd / (tag + '.time')
    identity = dict(binary_sha256=expected, threads=threads, stdin=stdin)
    if receipt.exists():
        previous = json.loads(receipt.read_text())
        if all(previous.get(k) == v for k, v in identity.items()) and previous.get('completed'):
            assert sha(output) == previous['log_sha256']
            receipts.append(previous)
            return previous
        raise RuntimeError(f'Incomplete or incompatible previous stage: {receipt}; inspect before restarting')
    if output.exists() or timing.exists():
        raise RuntimeError(f'Unreceipted stage artifacts require inspection: {tag}')
    env = {**os.environ, 'OMP_NUM_THREADS': str(threads), 'OMP_DYNAMIC': 'FALSE',
           'OMP_PROC_BIND': 'close', 'OMP_PLACES': 'cores',
           'LD_LIBRARY_PATH': os.environ['BDM_AUDIT_NATIVE_LIBS']}
    record = dict(**identity, binary=binary.name, started_at_utc=now(), completed=False,
                  driver_sha256=sha(__file__), job_id=os.environ['SLURM_JOB_ID'], host=socket.gethostname())
    write_json(receipt, record)
    start = time.monotonic()
    process = None
    try:
        with output.open('xb') as stream:
            process = subprocess.Popen(['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(timing), str(binary)],
                                       cwd=cwd, env=env, stdin=subprocess.PIPE, stdout=stream,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            process.communicate(stdin.encode(), timeout=4 * 3600)
        record.update(returncode=process.returncode, elapsed_seconds=time.monotonic() - start,
                      finished_at_utc=now(), log_sha256=sha(output))
        values = timing.read_text().splitlines()[-1].split()
        record.update(cpu_user_seconds=float(values[1]), cpu_system_seconds=float(values[2]),
                      maxrss_kib=int(values[3]))
        assert process.returncode == 0, f'{binary.name} failed; see {output}'
        record['completed'] = True
    except BaseException:
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        record['failure_traceback'] = traceback.format_exc()
        raise
    finally:
        write_json(receipt, record)
        receipts.append(record)
    print(tag, 'completed', round(record['elapsed_seconds'], 2), 's', record['maxrss_kib'], 'KiB', flush=True)
    return record


def read_header(path):
    raw = path.read_bytes()
    size = struct.unpack('>i', raw[:4])[0]
    assert len(raw) == size + 8 and raw[-4:] == raw[:4]
    # HEADER is 45 characters followed by 4-byte reals/integers, then int64 Nparticles.
    offset = 4 + 45
    scale = struct.unpack_from('>f', raw, offset)[0]
    step = struct.unpack_from('>i', raw, offset + 4 * 4)[0]
    nrow, ngrid = struct.unpack_from('>ii', raw, offset + 12 * 4)
    nparticles = struct.unpack_from('>q', raw, offset + 19 * 4)[0]
    return dict(scale_factor=scale, step=step, nrow=nrow, ngrid=ngrid, particles=nparticles)


def load_catalogue(path, spec, strict=True):
    data = np.loadtxt(path, skiprows=8, ndmin=2)
    if not data.size:
        data = np.empty((0, 24))
    assert data.shape[1] == 24
    record = dict(path=str(path.relative_to(ROOT)), sha256=sha(path), rows=len(data),
                  nonfinite_by_column=np.sum(~np.isfinite(data), axis=0).tolist())
    if strict:
        assert np.isfinite(data).all()
        assert np.all((data[:, :3] >= 0) & (data[:, :3] < spec['box_mpc_h']))
        assert np.all(data[:, 6] <= data[:, 7]) and np.all(data[:, 8] > 0)
        assert np.all(data[:, 10] >= 0)
        assert len(np.unique(data[:, 11])) == len(data), 'Duplicate published candidate ID'
    return data, record


def simulate(name, threads):
    spec = case_spec(name)
    case = WORK / name
    case.mkdir(exist_ok=True)
    run = case / 'Run1'
    run.mkdir(exist_ok=True)
    (run / 'CATALOGS').mkdir(exist_ok=True)
    intended = make_init(spec)
    for path, content in [(case / 'Init.dat', intended), (run / 'BDM.config', CONFIG)]:
        if path.exists():
            assert path.read_text() == content
        else:
            path.write_text(content)
    for filename in ['PkTable.dat', 'TableSeeds.dat']:
        path = case / filename
        target = WORK / 'reference-inputs' / filename
        if path.is_symlink():
            assert path.resolve() == target
        elif path.exists():
            assert sha(path) == sha(target)
        else:
            path.symlink_to(target)
    records = []
    report = dict(started_at_utc=now(), job_id=os.environ['SLURM_JOB_ID'], spec=spec,
                  inputs_sha256={p.name: sha(p) for p in [case / 'Init.dat', case / 'PkTable.dat',
                                  case / 'TableSeeds.dat', run / 'BDM.config']},
                  records=records, completed=False)
    try:
        run_native('PMP2init.exe', case, '', threads, 'init', records)
        assert (case / 'Setup.dat').is_file()
        report['setup'] = (case / 'Setup.dat').read_text()
        run_native('PMP2start.exe', run, '1\n', threads, 'start', records)
        initial = read_header(run / 'PMcrd.DAT')
        assert initial['nrow'] == spec['nrow'] and initial['ngrid'] == spec['ngrid']
        assert initial['particles'] == spec['nrow']**3 and abs(initial['scale_factor'] - 1/101) < 1.e-8
        report['initial_header'] = initial
        run_native('PMP2main.exe', run, '2000\n', threads, 'main', records)
        headers = sorted(run.glob('PMcrd.*.DAT'))
        assert len(headers) == len(spec['epochs']), f'Wrong output count: {headers}'
        snapshots = []
        arrays = {}
        for z, header in zip(spec['epochs'], headers):
            info = read_header(header)
            assert abs(info['scale_factor'] - 1 / (1 + z)) < 2.e-7, (z, info)
            assert info['particles'] == spec['nrow']**3
            step = info['step']
            files = sorted(run.glob(f'PMcr*.{step:04d}.DAT'))
            expected_data = 24 * spec['nrow']**3
            assert sum(p.stat().st_size for p in files if p.name.startswith('PMcrs')) == expected_data
            # Catshort names contain the integration step and realization.
            cats = list((run / 'CATALOGS').glob(f'Catshort*{step:04d}.0001.DAT'))
            if len(cats) != 1:
                raise ValueError(f'Cannot identify z={z} catalogue: {cats}')
            data, catalogue = load_catalogue(cats[0], spec)
            arrays[f'new_z{z}'] = data
            snapshots.append(dict(z=z, header=info, files_sha256={p.name: sha(p) for p in files},
                                  catalogue=catalogue))
            print(name, 'verified snapshot', z, step, len(data), 'haloes', flush=True)
        np.savez_compressed(ROOT / f'{name}-inline-catalogues.npz', **arrays)
        report.update(snapshots=snapshots, completed=True, finished_at_utc=now())
    except BaseException:
        report['failure_traceback'] = traceback.format_exc()
        raise
    finally:
        write_json(ROOT / f'{name}-simulation.json', report)


def check_memberships(path, spec):
    """Stream each membership once; memory is O(number of haloes + largest halo)."""
    file_size = path.stat().st_size
    offsets = []
    indices = []
    counts = []
    properties = []
    seen = {}
    duplicate_pairs = []
    with path.open('rb') as f:
        selected, candidates, mass_one = struct.unpack('>qqf', f.read(20))
        assert 0 <= selected <= candidates <= spec['nrow']**3
        for _ in range(selected):
            candidate, count = struct.unpack('>qq', f.read(16))
            values = np.frombuffer(f.read(84), dtype='>f4').astype(np.float64)
            offset = f.tell()
            assert 20 <= count <= spec['nrow']**3 and offset + count * 8 <= file_size
            raw = f.read(count * 8)
            ids = np.frombuffer(raw, dtype='>i8')
            assert np.isfinite(values).all() and ids.min() >= 1 and ids.max() <= spec['nrow']**3
            assert np.all(ids[1:] > ids[:-1]), 'Repeated original particle IDs within a halo'
            assert abs(values[6] - count * mass_one) <= 2 * abs(float(np.spacing(np.float32(values[6]))))
            key = (count, hashlib.sha256(raw).digest())
            if key in seen:
                previous_offset, previous_candidate = seen[key]
                cursor = f.tell()
                f.seek(previous_offset)
                if f.read(count * 8) == raw:
                    duplicate_pairs.append([previous_candidate, candidate])
                f.seek(cursor)
            seen[key] = (offset, candidate)
            offsets.append(offset)
            indices.append(candidate)
            counts.append(count)
            properties.append(values)
        assert f.tell() == file_size
    assert not duplicate_pairs, duplicate_pairs[:10]
    props = np.asarray(properties).reshape(-1, 21)
    index = np.asarray(indices, dtype=np.int64)
    box = spec['box_mpc_h']
    centre = props[:, :3] % box
    tree = cKDTree(centre, boxsize=box)
    violations = []
    examined = 0
    for start in range(0, selected, 2048):
        neighbours = tree.query_ball_point(centre[start:start+2048], props[start:start+2048, 8], workers=1)
        for i, matches in enumerate(neighbours, start):
            for j in matches:
                if i == j or (props[i, 6], -index[i]) <= (props[j, 6], -index[j]):
                    continue
                delta = centre[i] - centre[j]
                delta -= box * np.rint(delta / box)
                examined += 1
                if np.dot(delta, delta) < props[i, 8]**2:
                    violations.append([int(index[i]), int(index[j])])
    assert not violations, violations[:10]
    np.savez_compressed(path.with_suffix('.index.npz'), candidates=index, counts=np.asarray(counts),
                        offsets=np.asarray(offsets), properties=props)
    return dict(selected=selected, candidates=candidates, mass_one=mass_one,
                exact_duplicate_member_sets=0, repeated_original_ids=0, mass_count_mismatches=0,
                host_exclusion_violations=0, higher_priority_neighbours_examined=examined,
                raw_sha256=sha(path), retained_raw=str(path.relative_to(ROOT)),
                total_memberships=int(sum(counts)), index_sha256=sha(path.with_suffix('.index.npz')))


def validate(name, threads):
    spec = case_spec(name)
    simulation = json.loads((ROOT / f'{name}-simulation.json').read_text())
    assert simulation['completed']
    snapshots = simulation['snapshots']
    run = WORK / name / 'Run1'
    records = []
    arrays = {}
    report = dict(started_at_utc=now(), job_id=os.environ['SLURM_JOB_ID'], spec=spec,
                  simulation_receipt_sha256=sha(ROOT / f'{name}-simulation.json'), records=records, completed=False)
    try:
        for snapshot in snapshots:
            z = snapshot['z']
            step = snapshot['header']['step']
            for filename, expected in snapshot['files_sha256'].items():
                assert sha(run / filename) == expected
            inline = ROOT / snapshot['catalogue']['path']
            assert sha(inline) == snapshot['catalogue']['sha256']
            arrays[f'new_z{z}'], _ = load_catalogue(inline, spec)
            tasks = [('members', threads), ('old', threads)]
            if z == 0:
                tasks.extend([('baseline', max(1, threads//2)), ('state', threads)])
            for variant, nthreads in tasks:
                destination = WORK / name / f'replay-z{z}-{variant}-t{nthreads}'
                destination.mkdir(exist_ok=True)
                (destination / 'CATALOGS').mkdir(exist_ok=True)
                config = destination / 'BDM.config'
                if config.exists():
                    assert config.read_text() == CONFIG
                else:
                    config.write_text(CONFIG)
                for filename in snapshot['files_sha256']:
                    target = destination / filename
                    if target.is_symlink():
                        assert target.resolve() == run / filename
                    else:
                        target.symlink_to(run / filename)
                binary = 'PMP2BDM.exe' if variant == 'baseline' else f'PMP2BDM.{variant}.exe'
                stage = run_native(binary, destination, f'{step}\n', nthreads, variant, records)
                cats = list((destination / 'CATALOGS').glob('Catshort*.DAT'))
                assert len(cats) == 1
                data, cat = load_catalogue(cats[0], spec, strict=variant != 'old')
                stage.update(z=z, variant=variant, catalogue=cat)
                if variant == 'old':
                    arrays[f'old_z{z}'] = data
                else:
                    assert cat['sha256'] == snapshot['catalogue']['sha256'], (variant, 'inline/replay mismatch')
                    stage['byte_identical_to_inline'] = True
                if variant == 'members':
                    stage['membership'] = check_memberships(destination / 'repair-members.bin', spec)
                    assert stage['membership']['selected'] == len(data)
                if variant == 'state':
                    assert 'TWO CALLS BITWISE IDENTICAL' in (destination / 'state.log').read_text()
                    stage['two_calls_bitwise_identical'] = True
                write_json(ROOT / f'{name}-validation.json', report)
        np.savez_compressed(ROOT / f'{name}-comparison-catalogues.npz', **arrays)
        report.update(completed=True, finished_at_utc=now(),
                      comparison_sha256=sha(ROOT / f'{name}-comparison-catalogues.npz'))
    except BaseException:
        report['failure_traceback'] = traceback.format_exc()
        raise
    finally:
        write_json(ROOT / f'{name}-validation.json', report)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['prepare', 'simulate', 'validate', 'pilot'])
    parser.add_argument('--case', choices=['pilot', 'main'], default='pilot')
    parser.add_argument('--box', type=float, default=512.)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    if args.phase == 'prepare':
        prepare(args.box)
        return
    assert 'SLURM_JOB_ID' in os.environ, 'Large allocations require Slurm'
    assert 1 <= args.threads <= int(os.environ['SLURM_CPUS_PER_TASK'])
    if args.phase in ['simulate', 'pilot']:
        simulate('pilot' if args.phase == 'pilot' else args.case, args.threads)
    if args.phase in ['validate', 'pilot']:
        validate('pilot' if args.phase == 'pilot' else args.case, args.threads)


if __name__ == '__main__':
    main()
