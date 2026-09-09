"""Small native-versus-table evolution control (16^3 particles / 32^3 mesh)."""
import json
from pathlib import Path
import re
import shutil
import struct
import sys

import numpy as np

from common import BIN, ROOT, WORK, checker, file_manifest, now, write_json
from campaign import initialize, link, run_native, spec_for


def particles(path,n):
    # Read only populated rows; unused direct-record padding is not simulation data.
    data=np.memmap(path,mode='r',dtype='>f4',shape=(6,1024**2))
    return np.array(data[:,:n])


def main():
    root=WORK/'schedule-preflight';root.mkdir(parents=True,exist_ok=False)
    spec=spec_for('C');spec.update(nrow=16,ngrid=32,box_mpc_h=32.,case='tiny')
    runs={}
    for label in ['native','table','half']:
        test=dict(spec,half_steps=label=='half')
        case=root/label;run=initialize(case,test,1,allow_login=True)
        if label=='table':
            for name in ['PMcrd.DAT','PMcrs0.DAT']:link(run/name,runs['native']/name)
        else:
            run_native(BIN/'PMP2start.native.exe',run,'ic',1,
                       [case/'Setup.dat',case/'PkTable.dat',case/'TableSeeds.dat'],
                       lambda:[run/'PMcrd.DAT',run/'PMcrs0.DAT'],stdin='1\n',allow_login=True)
        binary='PMP2main.native.exe' if label=='native' else 'PMP2main.schedule.exe'
        run_native(BIN/binary,run,'main',1,
                   [case/'Setup.dat',run/'PMcrd.DAT',run/'PMcrs0.DAT',run/'campaign_schedule.dat'],
                   lambda:list(run.glob('PMcr*.[0-9][0-9][0-9][0-9].DAT')),stdin='1000\n',allow_login=True)
        runs[label]=run
    controls=[]
    for step in [92,109,158]:
        names=[f'PMcrd.{step:04d}.DAT',f'PMcrs0.{step:04d}.DAT']
        left,right=runs['native'],runs['table']
        assert (left/names[0]).read_bytes()==(right/names[0]).read_bytes()
        lhs=particles(left/names[1],16**3);rhs=particles(right/names[1],16**3)
        assert np.array_equal(lhs.view('u4'),rhs.view('u4'))
        half_header=checker().read_header(runs['half']/f'PMcrd.{step*2:04d}.DAT')
        normal_header=checker().read_header(left/names[0])
        assert half_header['scale_factor']==normal_header['scale_factor']
        controls.append(dict(normal_step=step,half_step=2*step,native_table_headers_bit_identical=True,
                             all_native_table_particle_words_bit_identical=True,exact_half_output_epoch=True))
    normal_ic=particles(runs['native']/'PMcrs0.DAT',16**3)
    half_ic=particles(runs['half']/'PMcrs0.DAT',16**3)
    assert np.array_equal(normal_ic[:3],half_ic[:3])
    assert not np.array_equal(normal_ic[3:],half_ic[3:])
    write_json(ROOT/'schedule-preflight.json',dict(completed=True,completed_at_utc=now(),controls=controls,
                half_initial_positions_identical=True,half_initial_velocities_differ=True,
                simulation_scope='16^3 particles, 32^3 mesh, one thread; all three requested epochs',
                evidence=file_manifest([p for run in runs.values() for p in run.glob('*.json')]),
                padding_excluded='Unused direct-record rows are excluded; all populated particle words are checked.'))
    print('Native/table evolution bit-identical; 316-step repeat reaches identical epochs',flush=True)


if __name__=='__main__':main()
