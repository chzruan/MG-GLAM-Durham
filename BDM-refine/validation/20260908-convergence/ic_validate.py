"""Production IC mode matching and all-row fine-mesh/timestep controls."""
import hashlib
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

from common import (MATRIX, ROOT, WORK, now, read_json_snapshot, sha,
                    verify_json_snapshot, verify_manifest, write_json)
from ic.configure import read_receipt


def source_identity(entry=None):
    """Identify the code actually executed, including frozen Slurm zipapps."""
    entry=Path(sys.argv[0] if entry is None else entry)
    names={'ic_validate.py':'__main__.py','common.py':'common.py',
           'ic/configure.py':'ic/configure.py'}
    if entry.is_file() and zipfile.is_zipfile(entry):
        with zipfile.ZipFile(entry) as archive:
            return {name:hashlib.sha256(archive.read(member)).hexdigest()
                    for name,member in names.items()}
    folder=Path(__file__).resolve().parent
    return {name:sha(folder/name) for name in names}


def main():
    started=now();sources=source_identity()
    snapshots={name:read_json_snapshot(ROOT/f'{name}-ic.json') for name in MATRIX}
    reports={name:payload for name,(payload,digest) in snapshots.items()}
    modes={};receipts={}
    for name,report in reports.items():
        assert report['completed'] and report['header']['particles']==MATRIX[name][0]**3<1200**3
        verify_manifest(report['outputs'])
        run=WORK/'cases'/name/'Run1'
        modes[name]=sha(run/'matched_modes.bin')
        receipts[name]=read_receipt(run/'matched_ic_receipt.txt')
    assert len(set(modes.values()))==1
    assert len({r['alpha'] for r in receipts.values()})==1
    ratio=float(receipts['T']['vcons'])/float(receipts['F']['vcons'])
    assert float(receipts['T']['a_velocity'])>float(receipts['F']['a_velocity'])
    max_position=0.;max_velocity_error=0.;count=0;position_sumsq=0.
    pure_tolerance=2*float(np.spacing(np.float32(4097)))*256/4096
    clamp_e=float(np.float32(2049)-np.float32(np.float32(2049)-np.float32(.001)))*.125
    clamp_f=float(np.float32(4097)-np.float32(np.float32(4097)-np.float32(.001)))*.0625
    native_tolerance=pure_tolerance+max(clamp_e,clamp_f)
    excess_count=nonlocal_excess=unexplained_excess=0;max_displacement=0.
    page=6*1024**2*4
    for part in range(4):
        paths=[WORK/f'cases/{name}/Run1/PMcrs{part}.DAT' for name in ['E','F','T']]
        with paths[0].open('rb') as es,paths[1].open('rb') as fs,paths[2].open('rb') as ts:
            while raw_e:=es.read(page):
                raw_f=fs.read(page);raw_t=ts.read(page)
                assert len(raw_e)==len(raw_f)==len(raw_t)==page
                e=np.frombuffer(raw_e,dtype='>f4').reshape(6,-1).astype(np.float64)
                f=np.frombuffer(raw_f,dtype='>f4').reshape(6,-1).astype(np.float64)
                t=np.frombuffer(raw_t,dtype='>f4').reshape(6,-1).astype(np.float64)
                assert np.isfinite(e).all() and np.isfinite(f).all() and np.isfinite(t).all()
                assert np.all((e[:3]>1)&(e[:3]<2049)) and np.all((f[:3]>1)&(f[:3]<4097))
                assert np.array_equal(f[:3],t[:3]), 'T initial positions differ from F'
                assert np.array_equal(2*e[3:],f[3:]), 'Fine-mesh IC physical velocities differ'
                delta=(e[:3]-1)*.125-(f[:3]-1)*.0625
                delta-=256*np.rint(delta/256)
                max_position=max(max_position,float(np.max(np.abs(delta))))
                position_sumsq+=float(np.sum(delta*delta))
                excess=np.abs(delta)>pure_tolerance
                excess_count+=int(np.count_nonzero(excess))
                xe=(e[:3]-1)*.125;xf=(f[:3]-1)*.0625
                near_e=np.minimum(xe,256-xe)<=native_tolerance
                near_f=np.minimum(xf,256-xf)<=native_tolerance
                nonlocal_excess+=int(np.count_nonzero(excess&~(near_e&near_f)))
                sentinel=(e[:3]==float(np.float32(2049)-np.float32(.001)))|(f[:3]==float(np.float32(4097)-np.float32(.001)))
                unexplained_excess+=int(np.count_nonzero(excess&~sentinel))
                displacement=np.abs(f[3:]*float(receipts['F']['xcons'])/float(receipts['F']['vcons'])*.0625)
                max_displacement=max(max_displacement,float(np.max(displacement)))
                residual=np.abs(t[3:]-ratio*f[3:])
                scale=np.maximum(np.abs(t[3:]),np.abs(ratio*f[3:]))
                assert np.all(residual<=4*np.finfo('f4').eps*scale+np.finfo('f4').tiny)
                relative=np.divide(residual,scale,out=np.zeros_like(residual),where=scale>0)
                max_velocity_error=max(max_velocity_error,float(np.max(relative)))
                count+=e.shape[1]
            assert not fs.read(1) and not ts.read(1)
    assert count==1024**3
    valid=max_position<=native_tolerance and nonlocal_excess==0 and unexplained_excess==0 and max_displacement<256
    for name,(payload,digest) in snapshots.items():
        verify_json_snapshot(ROOT/f'{name}-ic.json',digest)
    assert source_identity()==sources, 'Validator source changed during the check'
    write_json(ROOT/'ic-production-validation.json',dict(completed=valid,started_at_utc=started,completed_at_utc=now(),
        source_sha256=sources['ic_validate.py'],support_source_sha256=sources,
        input_receipt_sha256={name:digest for name,(payload,digest) in snapshots.items()},
        all_seven_mode_samples_identical=True,mode_sample_sha256=modes,shared_alpha=receipts['E']['alpha'],
        fine_rows_examined=count,F_T_initial_positions_all_identical=True,
        E_F_initial_physical_velocities_all_identical=True,E_F_max_position_difference_mpc_h=max_position,
        E_F_pure_float32_position_bound_mpc_h=pure_tolerance,
        E_F_native_periodic_mapper_bound_mpc_h=native_tolerance,
        E_F_position_component_rms_mpc_h=float(np.sqrt(position_sumsq/(3*count))),
        E_F_components_exceeding_pure_roundoff_bound=excess_count,
        E_F_nonedge_excess_components=nonlocal_excess,E_F_excess_without_exact_clamp_sentinel=unexplained_excess,
        max_inferred_displacement_mpc_h=max_displacement,F_T_staggered_velocity_ratio=ratio,
        F_T_max_relative_velocity_rounding_error=max_velocity_error,
        scope='Every fine-particle row for grid and timestep controls; common packed-mode sample across all seven runs',
        mode_matching_basis='Fixed master row RNG stride, per-plane seeds, master alpha and common physical origin; GNU/Intel full small-field controls',
        coordinate_caveat='Native BLOCKS upper-edge guard subtracts 1e-3 mesh units. This is an explicitly bounded periodic serialization effect, distinct from unmatched Fourier modes or pure float32 rounding.',
        boundary_source='PMP2start.f90:408-421; generated matched source retains the same guard'))
    assert valid, 'Production IC coordinate diagnostics exceed the native mapper bound/locality assumptions'
    print('All production IC controls passed; fine rows examined:',count,flush=True)


if __name__=='__main__':main()
