"""N1024 fixed-FI experiment. Run only on the sealed Slurm allocation.

Use micromamba run -n cosemu python3 -B. Native Fortran is copied unchanged
from numerical-threading, linked against the final production objects.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import socket
import struct
import subprocess
import time
import traceback

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np

ROOT = Path(os.environ.get('BDM_MAIN_THREAD_ROOT', Path(__file__).resolve().parent))
GROUPS = [['d64-t32-p1', 'd64-t64-p2', 'd64-t32-p3', 'd64-t64-p4'],
          ['d32-t32-p1', 'd32-t64-p2']]
TAGS = sum(GROUPS, [])


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while chunk := stream.read(8 * 1024**2):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    staged = path.with_suffix(path.suffix + '.tmp')
    with staged.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    staged.replace(path)


def read_header(path):
    raw = Path(path).read_bytes()
    assert len(raw) == 537 and raw[:4] == raw[-4:] == struct.pack('>i', 529), 'Invalid PM header record'
    base = 49
    return dict(scale_factor=struct.unpack_from('>f', raw, base)[0],
                step=struct.unpack_from('>i', raw, base+16)[0],
                nrow=struct.unpack_from('>i', raw, base+48)[0],
                ngrid=struct.unpack_from('>i', raw, base+52)[0],
                particles=struct.unpack_from('>q', raw, base+76)[0])


def select_snapshot(simulation, plan):
    assert simulation['completed'] and simulation['job_id'] == plan['parent_job_id']
    spec = simulation['spec']
    assert (spec['nrow'], spec['ngrid'], spec['box_mpc_h'], spec['epochs']) == (1024, 2048, 512., [2, 1, 0])
    assert spec['sources_sha256'] == plan['source_sha256']
    assert spec['binaries_sha256'] == plan['production_binaries_sha256']
    assert spec['halo_config'] == plan['config'] and spec['realization'] == 1 and spec['gravity'] == 'GR'
    snapshots = [item for item in simulation['snapshots'] if item['z'] == 0]
    assert len(snapshots) == 1, 'Require exactly one completed z=0 snapshot'
    snapshot = snapshots[0]
    header = snapshot['header']
    assert (header['nrow'], header['ngrid'], header['particles']) == (1024, 2048, 1024**3)
    assert header['particles'] < 1200**3 and header['scale_factor'] == 1. and header['step'] > 0
    step = header['step']
    assert set(snapshot['files_sha256']) == {f'PMcrd.{step:04d}.DAT', f'PMcrs0.{step:04d}.DAT', f'PMcrs1.{step:04d}.DAT'}
    main = [stage for stage in simulation['records'] if stage['binary'] == 'PMP2main.exe']
    assert len(main) == 1 and main[0]['completed'] and main[0]['returncode'] == 0
    assert main[0]['binary_sha256'] == plan['production_binaries_sha256']['PMP2main.exe']
    assert main[0]['threads'] == 64
    return snapshot


def restoration_markers(text):
    assert text.count('BDM FIXED DENSITY THREAD CONTROL COMPLETE') == 1, 'Missing native completion marker'
    fixed = re.findall(r'PROBE FIXED DENSITY VERIFIED\s+(\S+)\s+bit_differences=\s*(\d+)', text)
    particles = re.findall(r'PROBE PARTICLE RESTORATION VERIFIED\s+(\S+)', text)
    timings = re.findall(r'PROBE FINDER COMPLETE\s+(\S+)\s+seconds=\s*([\d.]+)', text)
    assert fixed == [(tag, '0') for tag in TAGS], 'Missing, duplicate, or nonzero fixed-FI restoration'
    assert particles == TAGS, 'Missing or duplicate particle-restoration marker'
    assert [tag for tag, _ in timings] == TAGS, 'Missing finder completion'
    return dict(fixed_density_verified=len(fixed), particle_restoration_verified=len(particles),
                tags=TAGS, exact_original_particle_bits=True), [dict(tag=tag, seconds=float(seconds)) for tag, seconds in timings]


def load_catalogue(path):
    data = np.loadtxt(path, skiprows=8, ndmin=2)
    assert data.shape[1] == 24 and len(data) > 0 and np.isfinite(data).all()
    assert np.all((data[:, :3] >= 0) & (data[:, :3] <= 512))
    assert np.all((data[:, 6] > 0) & (data[:, 6] <= data[:, 7])) and np.all(data[:, 8] > 0)
    assert np.all(data[:, 10] >= 0) and np.array_equal(data[:, 11], np.arange(1, len(data)+1))
    assert np.all((data[:, 20] >= 0) & (data[:, 20] <= data[:, 19]) & (data[:, 19] <= 1))
    return data


def compare_catalogues(first, second):
    first, second = Path(first), Path(second)
    a, b = load_catalogue(first), load_catalogue(second)
    answer = dict(first=str(first), second=str(second), first_sha256=sha(first), second_sha256=sha(second),
                  first_rows=len(a), second_rows=len(b))
    answer['byte_identical'] = answer['first_sha256'] == answer['second_sha256']
    if a.shape == b.shape:
        changed = np.flatnonzero(np.any(a != b, axis=1))
        answer.update(changed_rows=int(len(changed)),
                      row_alignment='published row order; candidate identity across distinct density fields is not assumed',
                      rows=[dict(row_one_based=int(i+1), columns_one_based=(np.flatnonzero(a[i] != b[i])+1).tolist(),
                                 first=a[i].tolist(), second=b[i].tolist()) for i in changed[:100]],
                      max_absolute_difference_by_column=np.max(np.abs(a-b), axis=0).tolist(),
                      masses_counts_velocities_identical=bool(np.array_equal(a[:, [3, 4, 5, 6, 13]], b[:, [3, 4, 5, 6, 13]])))
    return answer


def require_fixed_groups(comparisons):
    assert len(comparisons) == 4 and all(c['byte_identical'] for c in comparisons), 'Fixed-FI catalogue byte equality failed'


def read_peaks(path):
    with path.open('rb') as stream:
        count, grid, box = struct.unpack('>qif', stream.read(16))
    dtype = np.dtype([('candidate', '>i8'), ('x', '>f4'), ('y', '>f4'), ('z', '>f4'),
                      ('density', '>f4'), ('seed_radius', '>f4')])
    assert count > 0 and grid == 2048 and box == 512 and path.stat().st_size == 16 + count*dtype.itemsize
    data = np.fromfile(path, dtype=dtype, offset=16)
    assert np.array_equal(data['candidate'], np.arange(1, count+1))
    assert all(np.isfinite(data[column]).all() for column in dtype.names[1:])
    return data, dict(count=count, ngrid=grid, box_mpc_h=box, sha256=sha(path))


def normal_comparisons(validation, snapshot, plan, run):
    references = []
    inline = validation / snapshot['catalogue']['path']
    assert sha(inline) == snapshot['catalogue']['sha256']
    references.append(('main-inline-t64', inline))
    available, pending = [], []
    for directory in sorted((validation / 'work/main').glob('replay-z0-*-t*')):
        if not any(f'-{variant}-t' in directory.name for variant in ['baseline', 'members', 'state']):
            continue
        receipts = []
        for receipt in sorted(directory.glob('*.json')):
            value = json.loads(receipt.read_text())
            if value.get('binary') in ['PMP2BDM.exe', 'PMP2BDM.members.exe', 'PMP2BDM.state.exe']:
                receipts.append((receipt, value))
        if not receipts or not all(value.get('completed') for _, value in receipts):
            pending.append(directory.name)
            continue
        for receipt, value in receipts:
            assert value['binary_sha256'] == plan['production_binaries_sha256'][value['binary']]
            assert int(value['stdin'].strip()) == snapshot['header']['step'] and value['returncode'] == 0
            expected = {**snapshot['files_sha256'], 'BDM.config': plan['config_sha256']}
            assert {Path(name).name: digest for name, digest in value['inputs_sha256'].items()} == expected
            outputs = [(Path(name), digest) for name, digest in value['outputs_sha256'].items()
                       if Path(name).name.startswith('Catshort')]
            assert len(outputs) == 1
            path, digest = outputs[0]
            assert sha(path) == digest
            references.append((directory.name, path))
            available.append(dict(path=str(receipt), sha256=sha(receipt),
                                  scope='completed native stage; membership validation may still be running'))
    return dict(collected_at_utc=now(), available_receipts=available, pending_directories=pending,
                independent_of_job='11948467',
                comparisons=[dict(reference=tag, **compare_catalogues(path, run / (control+'.DAT')))
                             for tag, path in references for control in ['d64-t64-p2', 'd32-t64-p2']])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan-sha256', required=True)
    args = parser.parse_args()
    assert 'SLURM_JOB_ID' in os.environ and int(os.environ['SLURM_CPUS_PER_TASK']) == 64
    assert sha(ROOT / 'plan.json') == args.plan_sha256
    plan = json.loads((ROOT / 'plan.json').read_text())
    assert sha(__file__) == plan['runner_sha256']
    assert not (ROOT / 'results.json').exists(), 'Existing results require inspection; never overwrite or silently resume'
    run = ROOT / 'work/run'
    run.mkdir(exist_ok=False)
    (run / 'CATALOGS').mkdir()
    record = dict(started_at_utc=now(), completed=False, fixed_density_controls_passed=False,
                  job_id=os.environ['SLURM_JOB_ID'], host=socket.gethostname(), plan_sha256=args.plan_sha256,
                  driver_sha256=sha(__file__), protocol=plan['protocol'], source_sha256=plan['source_sha256'],
                  resources=plan['resources'], native_runtime_libraries=os.environ['BDM_AUDIT_NATIVE_LIBS'])
    process = None
    start = time.monotonic()
    write_json(ROOT / 'results.json', record)
    try:
        for filename, expected in plan['frozen_files_sha256'].items():
            assert sha(ROOT / filename) == expected, filename
        build = json.loads((ROOT / 'build.json').read_text())
        assert build['completed'] and build['source_sha256'] == plan['source_sha256']
        executable = ROOT / plan['executable']
        assert sha(executable) == build['binary_sha256']
        record.update(binary_sha256=sha(executable), build_receipt_sha256=sha(ROOT / 'build.json'))
        validation = Path(plan['validation_root'])
        simulation_path = validation / 'main-simulation.json'
        simulation_digest = sha(simulation_path)
        simulation = json.loads(simulation_path.read_text())
        assert sha(simulation_path) == simulation_digest
        snapshot = select_snapshot(simulation, plan)
        original = validation / 'work/main/Run1'
        step = snapshot['header']['step']
        assert read_header(original / f'PMcrd.{step:04d}.DAT') == snapshot['header']
        record.update(simulation_receipt_sha256=simulation_digest, snapshot=snapshot,
                      command=[str(executable)], stdin=f'{step} 64 32\n')
        for name, expected in snapshot['files_sha256'].items():
            path = original / name
            if name.startswith('PMcrs'):
                assert path.stat().st_size == 12*1024**3
            assert sha(path) == expected, f'Changed snapshot: {name}'
            (run / name).symlink_to(path)
        config = original / 'BDM.config'
        assert sha(config) == simulation['inputs_sha256']['BDM.config'] == plan['config_sha256']
        assert config.read_text() == plan['config']
        (run / 'BDM.config').write_bytes(config.read_bytes())
        record['config_sha256'] = sha(config)
        write_json(ROOT / 'results.json', record)
        env = dict(os.environ, LD_LIBRARY_PATH=os.environ['BDM_AUDIT_NATIVE_LIBS'],
                   OMP_NUM_THREADS='64', OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close', OMP_PLACES='cores')
        env.pop('LIBRARY_PATH', None)
        output, timing = run / 'probe.log', run / 'probe.time'
        native_started = time.monotonic()
        with output.open('xb') as log:
            process = subprocess.Popen(['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(timing), str(executable)],
                                       cwd=run, env=env, stdin=subprocess.PIPE, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            process.communicate(record['stdin'].encode(), timeout=max(1, 2400-(time.monotonic()-start)))
        record.update(returncode=process.returncode, elapsed_seconds=time.monotonic()-native_started,
                      log_sha256=sha(output))
        values = timing.read_text().splitlines()[-1].split()
        record.update(cpu_user_seconds=float(values[1]), cpu_system_seconds=float(values[2]), maxrss_kib=int(values[3]))
        assert process.returncode == 0, 'Native diagnostic failed'
        text = output.read_text()
        record['restoration'], record['native_finder_timings'] = restoration_markers(text)
        record['density_difference'] = re.findall(r'PROBE DENSITY DIFFERENCES[^\n]+', text)
        density = re.findall(r'PROBE DENSITY DIFFERENCES cells=(\d+) maximum_absolute=\s*(\S+)', text)
        assert len(density) == 1
        record['density_field_difference'] = dict(cells_changed=int(density[0][0]), max_absolute=float(density[0][1]))
        assert 0 <= int(density[0][0]) <= 2048**3 and np.isfinite(float(density[0][1]))
        record['fixed_density_comparisons'] = [compare_catalogues(run / (group[0]+'.DAT'), run / (other+'.DAT'))
                                               for group in GROUPS for other in group[1:]]
        record['fixed_density_groups'] = [dict(density_field=group[0].split('-')[0], tags=group,
                byte_identical=all(sha(run / (tag+'.DAT')) == sha(run / (group[0]+'.DAT')) for tag in group))
                for group in GROUPS]
        record['between_density_comparison'] = compare_catalogues(run / 'd64-t64-p2.DAT', run / 'd32-t64-p2.DAT')
        a, summary_a = read_peaks(run / 'd64.peaks.bin')
        b, summary_b = read_peaks(run / 'd32.peaks.bin')
        record['peaks'] = dict(d64=summary_a, d32=summary_b)
        if len(a) == len(b):
            record['peaks'].update(changed_records=int(np.count_nonzero(a != b)),
                positions_identical=all(np.array_equal(a[c], b[c]) for c in ['candidate', 'x', 'y', 'z']),
                densities_changed=int(np.count_nonzero(a['density'] != b['density'])),
                seed_radii_changed=int(np.count_nonzero(a['seed_radius'] != b['seed_radius'])))
        del a, b
        record['normal_comparisons'] = normal_comparisons(validation, snapshot, plan, run)
        record['density_tapes'] = {}
        for threads in [64, 32]:
            path = run / f'd{threads}.density.bin'
            with path.open('rb') as stream:
                assert struct.unpack('>qqi', stream.read(20)) == (2048, 1024**3, threads)
            assert path.stat().st_size == 20 + 4*2048**3
            record['density_tapes'][path.name] = dict(bytes=path.stat().st_size, path=str(path), sha256=sha(path))
        # Check inputs again after the native run; restored memory alone does
        # not establish unchanged files across a long shared-filesystem replay.
        assert sha(simulation_path) == simulation_digest
        for name, expected in snapshot['files_sha256'].items():
            assert sha(original / name) == expected
        assert sha(config) == plan['config_sha256']
        record['input_files_verified_after_run'] = True
        record['outputs_sha256'] = {p.name: record['density_tapes'][p.name]['sha256']
            if p.name in record['density_tapes'] else sha(p)
            for p in sorted(run.iterdir()) if p.is_file() and not p.is_symlink()}
        require_fixed_groups(record['fixed_density_comparisons'])
        record.update(completed=True, fixed_density_controls_passed=True, finished_at_utc=now(),
                      total_elapsed_seconds=time.monotonic()-start)
    except BaseException:
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        record['failure_traceback'] = traceback.format_exc()
        record['failed_at_utc'] = now()
        raise
    finally:
        write_json(ROOT / 'results.json', record)
    print('Fixed-density diagnostic complete:', record['fixed_density_controls_passed'], flush=True)


if __name__ == '__main__':
    main()
