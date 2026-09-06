"""Prove the original production membership checker detects a violating pair.

Run the small pilot on the login node; the full three-epoch test requires a
sized shared Slurm allocation. Original tapes, indexes and catalogues are read
only. Each temporary corrupt copy differs in exactly one halo centre, with the
catalogue row changed consistently so earlier property checks still pass.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import resource
import shutil
import struct
import tempfile
import time

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np

ROOT = Path(os.environ.get('BDM_REVIEW_ROOT',Path(__file__).resolve().parent.parent)).resolve()
HERE = ROOT/'host-control'
REPO = ROOT.parents[2]
VALIDATION = REPO/'BDM-refine/validation/20260906-n1024'
definition = importlib.util.spec_from_file_location('validated_driver', VALIDATION/'run_validation.py')
v = importlib.util.module_from_spec(definition)
definition.loader.exec_module(v)
# This affects only the returned relative pathname, not the checking arithmetic.
v.ROOT = REPO


def verify_rejection(path, spec, catalogue, host_row, moved_row, index):
    expected_pair = [int(index['candidates'][host_row]),int(index['candidates'][moved_row])]
    coordinates = np.asarray(index['properties'][host_row,:3],dtype='>f4')
    position = int(index['offsets'][moved_row])-84
    with path.open('r+b') as stream:
        stream.seek(position)
        original = stream.read(12)
        stream.seek(position)
        stream.write(coordinates.tobytes())
    corrupt_catalogue = catalogue.copy()
    corrupt_catalogue[moved_row,:3] = coordinates.astype(np.float64)
    try:
        v.check_memberships(path,spec,corrupt_catalogue)
    except AssertionError as exc:
        pairs = exc.args[0] if exc.args else None
        assert isinstance(pairs,list) and expected_pair in pairs, repr(exc)
    else:
        raise AssertionError('Production checker accepted the known violating pair')
    assert not path.with_suffix('.index.npz').exists(), 'Invalid data acquired a successful index'
    return dict(expected_host_and_lower_priority_candidate=expected_pair,
        host_row_zero_based=host_row,moved_row_zero_based=moved_row,
        property_byte_offset=position,original_coordinate_bytes_hex=original.hex(),
        injected_coordinate_bytes_hex=coordinates.tobytes().hex(),
        host_mass=float(index['properties'][host_row,6]),moved_mass=float(index['properties'][moved_row,6]),
        host_radius_mpc_h=float(index['properties'][host_row,8]),
        separation_mpc_h=0.,detected_pairs=pairs,expected_host_assertion_failed=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    if not args.pilot:
        assert os.environ.get('SLURM_JOB_ID'), 'Production-scale checker needs its sized allocation'
    output = args.output or HERE/('pilot.json' if args.pilot else 'results.json')
    assert not output.exists(), 'Preserve existing control evidence'
    simulation = json.loads((VALIDATION/'main-simulation.json').read_text())
    validation = json.loads((VALIDATION/'main-validation.json').read_text())
    assert simulation['completed'] and validation['completed']
    assert v.sha(VALIDATION/'main-simulation.json') == validation['simulation_receipt_sha256']
    spec = simulation['spec']
    epochs = [0] if args.pilot else [2,1,0]
    records = []
    start = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_SELF)
    scratch_root = ROOT/'work/host-control'
    scratch_root.mkdir(parents=True,exist_ok=True)
    report = dict(completed=False,pilot=args.pilot,job_id=os.environ.get('SLURM_JOB_ID'),
        started_at_utc=datetime.now(timezone.utc).isoformat(),checker_module_sha256=v.sha(VALIDATION/'run_validation.py'),
        checker_function_sha256=hashlib.sha256(inspect.getsource(v.check_memberships).encode()).hexdigest(),
        driver_sha256=v.sha(__file__),simulation_receipt_sha256=v.sha(VALIDATION/'main-simulation.json'),
        validation_receipt_sha256=v.sha(VALIDATION/'main-validation.json'),
        catalogue_arrays_sha256=v.sha(VALIDATION/'main-verified-comparison-catalogues.npz'),
        checker_change='Only returned path base set to repository root; original checker function unmodified',records=records)
    try:
        with np.load(VALIDATION/'main-verified-comparison-catalogues.npz',allow_pickle=False) as catalogues:
            for z in epochs:
                record = next(r for r in validation['records'] if r['variant']=='members' and r['z']==z)
                membership = record['membership']
                original = VALIDATION/membership['retained_raw']
                index_path = original.with_suffix('.index.npz')
                assert v.sha(index_path)==membership['index_sha256']
                with np.load(index_path,allow_pickle=False) as arrays:
                    index = {k:arrays[k] for k in arrays.files}
                catalogue = catalogues[f'new_z{z}']
                assert len(catalogue)==membership['selected']
                with tempfile.TemporaryDirectory(prefix=f'z{z}-',dir=scratch_root) as temporary:
                    directory = Path(temporary)
                    if args.pilot:
                        # Two original rows are sufficient to size/control the test
                        # without reading a billion-particle membership tape here.
                        base = directory/'pilot.bin'
                        with original.open('rb') as src,base.open('xb') as dst:
                            _,candidates,mass_one=struct.unpack('>qqf',src.read(20))
                            dst.write(struct.pack('>qqf',2,candidates,mass_one))
                            offsets=[]
                            for j in range(2):
                                src.seek(int(index['offsets'][j])-100)
                                raw=src.read(100+8*int(index['counts'][j]))
                                offsets.append(dst.tell()+100)
                                dst.write(raw)
                        index = {k:a[:2].copy() for k,a in index.items()}
                        index['offsets'] = np.asarray(offsets,dtype=np.int64)
                        catalogue = catalogue[:2].copy()
                        catalogue[:,11] = np.arange(1,3)
                        baseline = v.check_memberships(base,spec,catalogue)
                    else:
                        assert v.sha(original)==membership['raw_sha256']
                        base=original
                        baseline=v.check_memberships(base,spec,catalogue)
                        assert baseline['raw_sha256']==membership['raw_sha256']
                    assert baseline['host_exclusion_violations']==0
                    host_row=int(np.argmax(index['properties'][:,6]))
                    moved_row=int(np.argmin(index['properties'][:,6]))
                    if moved_row==host_row:
                        host_row,moved_row=0,1
                    assert (index['properties'][host_row,6],-index['candidates'][host_row]) > (
                        index['properties'][moved_row,6],-index['candidates'][moved_row])
                    corrupted=directory/'violating.bin'
                    shutil.copyfile(base,corrupted)
                    injection=verify_rejection(corrupted,spec,catalogue,host_row,moved_row,index)
                    records.append(dict(z=z,selected=len(catalogue),baseline=baseline,injection=injection,
                        original_raw=str(original.relative_to(REPO)),original_raw_sha256=membership['raw_sha256'],
                        original_index_sha256=membership['index_sha256']))
                assert v.sha(index_path)==membership['index_sha256']
                if not args.pilot:
                    assert v.sha(original)==membership['raw_sha256']
                print(f'z={z}: {len(catalogue)} rows passed; injected host pair was rejected',flush=True)
        report.update(completed=True,original_data_and_indexes_unchanged=True,
                      temporary_copies_removed=True,finished_at_utc=datetime.now(timezone.utc).isoformat())
    finally:
        after=resource.getrusage(resource.RUSAGE_SELF)
        report.update(elapsed_seconds=time.monotonic()-start,
            total_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            maxrss_kib=after.ru_maxrss)
        v.write_json(output,report)


if __name__=='__main__':
    main()
