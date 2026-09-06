"""Run the native fixed-density control on Slurm and preserve exact evidence."""
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

ROOT = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while chunk := f.read(8 * 1024**2):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    staging = path.with_suffix(path.suffix + '.tmp')
    with staging.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    staging.replace(path)


def compare_catalogues(first, second):
    a = np.loadtxt(first, skiprows=8, ndmin=2)
    b = np.loadtxt(second, skiprows=8, ndmin=2)
    assert a.shape[1] == b.shape[1] == 24 and np.isfinite(a).all() and np.isfinite(b).all()
    answer = dict(first=str(first), second=str(second), first_sha256=sha(first), second_sha256=sha(second),
                  first_rows=len(a), second_rows=len(b), byte_identical=sha(first) == sha(second))
    if a.shape == b.shape:
        different = np.flatnonzero(np.any(a != b, axis=1))
        answer.update(changed_rows=len(different),
                      rows=[dict(row_one_based=int(i+1), columns_one_based=(np.flatnonzero(a[i] != b[i])+1).tolist(),
                                 first=a[i].tolist(), second=b[i].tolist()) for i in different[:100]],
                      max_absolute_difference_by_column=np.max(np.abs(a-b), axis=0).tolist(),
                      masses_counts_velocities_identical=bool(np.array_equal(a[:, [3,4,5,6,13]], b[:, [3,4,5,6,13]])))
    return answer


def read_peaks(path):
    with path.open('rb') as f:
        count, grid, box = struct.unpack('>qif', f.read(16))
    dtype = np.dtype([('candidate', '>i8'), ('x', '>f4'), ('y', '>f4'), ('z', '>f4'),
                      ('density', '>f4'), ('seed_radius', '>f4')])
    assert path.stat().st_size == 16 + count * dtype.itemsize
    return np.fromfile(path, dtype=dtype, offset=16), dict(count=count, ngrid=grid, box_mpc_h=box, sha256=sha(path))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--original-repo', type=Path, required=True)
    args = parser.parse_args()
    assert 'SLURM_JOB_ID' in os.environ, 'All native allocations belong on Slurm'
    assert int(os.environ['SLURM_CPUS_PER_TASK']) >= 32
    build = json.loads((ROOT / 'build.json').read_text())
    assert build['completed']
    executable = ROOT / 'work/build/BDM-thread-probe.exe'
    assert sha(executable) == build['binary_sha256']
    for name, expected in build['probe_source_sha256'].items():
        assert sha(ROOT / name) == expected
    validation = args.original_repo.resolve() / 'BDM-refine/validation/20260906-n1024'
    simulation_receipt = validation / 'pilot-simulation.json'
    simulation = json.loads(simulation_receipt.read_text())
    assert simulation['completed'] and simulation['spec']['nrow'] == 512 and simulation['spec']['ngrid'] == 1024
    snapshot = next(item for item in simulation['snapshots'] if item['z'] == 0)
    assert snapshot['header']['particles'] == 512**3 < 1200**3
    assert snapshot['header']['step'] == 157
    original = validation / 'work/pilot/Run1'
    run = ROOT / 'work/run'
    run.mkdir(exist_ok=False)
    (run / 'CATALOGS').mkdir()
    record = dict(started_at_utc=datetime.now(timezone.utc).isoformat(), completed=False,
                  job_id=os.environ['SLURM_JOB_ID'], host=socket.gethostname(),
                  build_receipt_sha256=sha(ROOT / 'build.json'), driver_sha256=sha(__file__),
                  binary_sha256=sha(executable), simulation_receipt_sha256=sha(simulation_receipt),
                  source_sha256=build['source_sha256'], snapshot=snapshot,
                  command=[str(executable)], stdin='157 32 16\n',
                  protocol='DENSIT once at32; identical saved FI -> BDM(0) at16,32,16,32; DENSIT once at16; identical second FI -> BDM(0) at16,32',
                  native_runtime_libraries=os.environ['BDM_AUDIT_NATIVE_LIBS'])
    output = run / 'probe.log'
    timing = run / 'probe.time'
    process = None
    try:
        for name, expected in snapshot['files_sha256'].items():
            assert sha(original / name) == expected
            (run / name).symlink_to(original / name)
        config = original / 'BDM.config'
        assert config.read_text() == 'iVirial=1\nMassMin=2.5e12\nRext=0.15\n'
        (run / 'BDM.config').write_bytes(config.read_bytes())
        record['config_sha256'] = sha(config)
        write_json(ROOT / 'results.json', record)
        env = dict(os.environ, LD_LIBRARY_PATH=os.environ['BDM_AUDIT_NATIVE_LIBS'],
                   OMP_NUM_THREADS='32', OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close', OMP_PLACES='cores')
        env.pop('LIBRARY_PATH', None)
        start = time.monotonic()
        with output.open('xb') as log:
            process = subprocess.Popen(['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(timing), str(executable)],
                                       cwd=run, env=env, stdin=subprocess.PIPE, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            process.communicate(record['stdin'].encode(), timeout=245)
        record.update(returncode=process.returncode, elapsed_seconds=time.monotonic()-start,
                      log_sha256=sha(output))
        values = timing.read_text().splitlines()[-1].split()
        record.update(cpu_user_seconds=float(values[1]), cpu_system_seconds=float(values[2]), maxrss_kib=int(values[3]))
        assert process.returncode == 0, 'Native diagnostic failed'
        text = output.read_text()
        assert text.count('BDM FIXED DENSITY THREAD CONTROL COMPLETE') == 1
        assert text.count('PROBE FIXED DENSITY VERIFIED') == text.count('PROBE PARTICLE RESTORATION VERIFIED') == 6
        record['density_difference'] = re.findall(r'PROBE DENSITY DIFFERENCES[^\n]+', text)
        record['native_finder_timings'] = re.findall(r'PROBE FINDER COMPLETE[^\n]+', text)
        groups = [['d32-t16-p1', 'd32-t32-p2', 'd32-t16-p3', 'd32-t32-p4'], ['d16-t16-p1', 'd16-t32-p2']]
        comparisons = []
        for group in groups:
            for other in group[1:]:
                comparisons.append(compare_catalogues(run / (group[0]+'.DAT'), run / (other+'.DAT')))
        record['fixed_density_comparisons'] = comparisons
        record['fixed_density_controls_passed'] = all(c['byte_identical'] for c in comparisons)
        record['between_density_comparison'] = compare_catalogues(run / 'd32-t32-p2.DAT', run / 'd16-t32-p2.DAT')
        # Paths in the original driver are relative to its validation root.
        inline = validation / snapshot['catalogue']['path']
        old16 = next((validation / 'work/pilot/replay-z0-baseline-t16/CATALOGS').glob('Catshort*.DAT'))
        record['original_comparison'] = compare_catalogues(inline, old16)
        record['against_original'] = [compare_catalogues(reference, run / (tag+'.DAT'))
                                      for reference in [inline, old16] for tag in ['d32-t32-p2', 'd16-t32-p2']]
        peaks32, summary32 = read_peaks(run / 'd32.peaks.bin')
        peaks16, summary16 = read_peaks(run / 'd16.peaks.bin')
        record['peaks'] = {'d32': summary32, 'd16': summary16}
        if len(peaks32) == len(peaks16):
            changed = np.flatnonzero(peaks32 != peaks16)
            record['peaks']['changed_records'] = len(changed)
            record['peaks']['positions_identical'] = all(np.array_equal(peaks32[c], peaks16[c]) for c in ['x','y','z'])
            record['peaks']['seed_radii_changed'] = int(np.count_nonzero(peaks32['seed_radius'] != peaks16['seed_radius']))
            source_index = np.load(validation / 'work/pilot/replay-z0-members-t32/repair-members.index.npz')
            tracked = [source_index['candidates'][item['row_one_based']-1] for item in record['original_comparison'].get('rows', [])]
            record['peaks']['original_affected_candidates'] = [
                dict(candidate=int(candidate),
                     d32={name:float(peaks32[name][candidate-1]) for name in ['x','y','z','density','seed_radius']},
                     d16={name:float(peaks16[name][candidate-1]) for name in ['x','y','z','density','seed_radius']}) for candidate in tracked]
        record['outputs_sha256'] = {p.name:sha(p) for p in sorted(run.iterdir()) if p.is_file() and not p.is_symlink()}
        record.update(completed=True, finished_at_utc=datetime.now(timezone.utc).isoformat())
    except BaseException:
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        record['failure_traceback'] = traceback.format_exc()
        raise
    finally:
        write_json(ROOT / 'results.json', record)
    print('Fixed-density diagnostic complete:', record['fixed_density_controls_passed'], flush=True)


if __name__ == '__main__':
    main()
