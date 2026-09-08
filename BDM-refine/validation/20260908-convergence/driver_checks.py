"""Small membership-matching and real replay publication/resume controls."""
import json
from pathlib import Path
import sys

import numpy as np

import analyze
import campaign
from common import ROOT, WORK, file_manifest, load_module, now, sha, write_json
import replays


def matching(root):
    expected=np.arange(1,65,dtype=np.int64)
    assert np.array_equal(analyze.project_ids(np.arange(1,513),8,4),expected)
    assert analyze.project_ids(np.array([439]),8,4).tolist()==[64]
    samples=[]
    for label,n,groups,centres in [
        ('left',4,[np.arange(1,33),np.arange(33,65)],[[.1,2,2],[10.05,10,10]]),
        ('right',8,[np.arange(257,513),np.arange(1,257)],[[10,10,10],[31.99,2,2]])]:
        raw=root/(label+'.bin');offsets=[]
        with raw.open('wb') as stream:
            for ids in groups:offsets.append(stream.tell());stream.write(ids.astype('>i8').tobytes())
        props=np.zeros((2,21));props[:,:3]=centres;props[:,6:8]=1.e13;props[:,8]=.5
        sample=dict(name=label,z=0,raw=raw,properties=props,counts=np.array(list(map(len,groups))),offsets=np.array(offsets),
                    report=dict(spec=dict(nrow=n,box_mpc_h=32.),
                                variants=dict(v3=dict(membership=dict(mass_one=1.e13/len(groups[0]))))))
        samples.append(sample)
    pairs,info=analyze.match(*samples)
    assert pairs[:,:2].astype(int).tolist()==[[0,1],[1,0]],pairs
    assert np.all(pairs[:,2]==1) and np.all(pairs[:,4:6]==1)
    assert abs(pairs[0,6]-.11)<1.e-12
    statistics=analyze.compare(*samples,pairs,'control',floor=20)
    json.dumps(analyze.finite_json(statistics),allow_nan=False)
    return dict(nested_lattice_corner=True,full_fine_lattice_projects_exactly=True,
                permuted_halo_ids_matched_by_membership=True,periodic_boundary_pair=True,
                empty_statistical_bins_json_safe=True)


def replay_publication(root):
    sys.path.insert(0,str(ROOT/'replay'))
    helper=load_module('adapter_fixtures',ROOT/'replay/test_adapter.py')
    source=root/'snapshot';helper.fixture(source,16)
    source_files=[source/'PMcrd.0001.DAT',source/'PMcrs0.0001.DAT']
    spec=dict(case='fixture',nrow=4,ngrid=16,box_mpc_h=32.,epochs=[0],half_steps=False,
              cosmology=dict(Omega_m=.3),shape_repair_commit='fixture')
    simulation=dict(completed=True,spec=spec,outputs=file_manifest(source_files),snapshots=[dict(
        redshift=0,header=dict(step=1,scale_factor=1.),header_path=str(source_files[0]),data_paths=[str(source_files[1])])])
    original_root,original_work,original_run=replays.ROOT,replays.WORK,replays.run_native
    replays.ROOT=root/'reports';replays.ROOT.mkdir();replays.WORK=root
    # Tiny native fixtures are explicitly recorded with job_id=None, never a fake Slurm ID.
    def native(*args,**kwargs):return campaign.run_native(*args,**kwargs,allow_login=True)
    replays.run_native=native
    try:
        write_json(replays.ROOT/'fixture-simulation.json',simulation)
        replays.replay('fixture',1,epoch=0,finder='v3',analysis_ngrid=16)
        replays.replay('fixture',1,epoch=0,finder='pair',analysis_ngrid=16)
        receipt=replays.ROOT/'replay-fixture-z0-ng16.json'
        assert json.loads(receipt.read_text())['pair_completed']
        before=sha(receipt)
        replays.replay('fixture',1,epoch=0,finder='v3',analysis_ngrid=16)
        replays.replay('fixture',1,epoch=0,finder='pair',analysis_ngrid=16)
        assert sha(receipt)==before,'Completed comparison was downgraded or rewritten'
        density=root/'replays/fixture-z0-ng16/density.bin'
        with density.open('r+b') as stream:
            original=stream.read(1);stream.seek(0);stream.write(bytes([original[0]^1]))
        try:replays.replay('fixture',1,epoch=0,finder='pair',analysis_ngrid=16)
        except ValueError:pass
        else:raise AssertionError('Changed completed density was accepted')
        with density.open('r+b') as stream:stream.write(original)
        stop=root/'legacy-stop.sh'
        stop.write_text("#!/bin/bash\nprintf '%s\\n' 'Partial legacy catalogue header' > CATALOGS/Catshort.partial.DAT\n")
        stop.chmod(0o755)
        def legacy_stop(binary,*args,**kwargs):
            return native(stop if Path(binary).name=='PMP2replay.legacy.exe' else binary,*args,**kwargs)
        replays.run_native=legacy_stop
        write_json(replays.ROOT/'stopped-simulation.json',simulation)
        replays.replay('stopped',1,epoch=0,finder='pair',analysis_ngrid=16)
        saved=json.loads((replays.ROOT/'replay-stopped-z0-ng16.json').read_text())
        assert saved['completed'] and not saved['pair_completed'] and not saved['variants']['legacy']['completed']
        assert (replays.ROOT/'v3-validation-stopped-z0-ng16.json').is_file()
        return dict(real_native_pair=True,completed_pair_preserved_on_both_resume_modes=True,
                    changed_retained_density_rejected=True,legacy_exit_zero_partial_catalogue_preserves_v3=True)
    finally:replays.ROOT,replays.WORK,replays.run_native=original_root,original_work,original_run


def main():
    root=WORK/'driver-controls';root.mkdir(parents=True,exist_ok=False)
    result=dict(completed=False,started_at_utc=now(),matching=matching(root))
    result['replay']=replay_publication(root)
    result.update(completed=True,completed_at_utc=now())
    write_json(ROOT/'driver-controls.json',result)
    print('Membership matching and replay driver controls passed',flush=True)


if __name__=='__main__':main()
