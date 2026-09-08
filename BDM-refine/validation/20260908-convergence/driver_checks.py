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


def statistics():
    masses=np.array([2.51e12,3.e12,4.e12,8.e12,1.1e13,2.e13])
    props=np.zeros((len(masses),21));props[:,6]=masses;props[:,7]=masses*1.2
    props[:,8]=.5;props[:,11]=120.;props[:,16:18]=.7
    sample=dict(name='reference',z=0,properties=props,counts=(masses/1.e10).astype(int),
                report=dict(spec=dict(box_mpc_h=256.),variants=dict(v3=dict(membership=dict(mass_one=1.e10)))))
    other=dict(sample,name='left',properties=props.copy())
    other['properties'][3,11]=0.;other['properties'][4,11]=132.;props[2,11]=0.
    pairs=np.zeros((len(masses),7));pairs[:,0]=pairs[:,1]=np.arange(len(masses));pairs[:,2:6]=1.
    density=analyze.abundance(sample)
    assert density['counts'][0]==2 and not density['publication_complete_bins'][0]
    assert np.isnan(density['dn_dlog10m'][0]),'Partially published first bin looked complete'
    result=analyze.compare(other,sample,pairs,'control',floor=300)
    assert not result['valid_mass_bins'][0] and np.isnan(result['abundance_ratio'][0])
    assert result['valid_mass_bins'][1]
    stats=result['matched_statistics'];assert 'so_mass' not in stats and 'aperture_total_mass' in stats
    for index in [2,3]:
        bin_index=np.searchsorted(analyze.EDGES,np.log10(masses[index]),side='right')-1
        assert stats['vmax'][bin_index]['count']==0
        assert result['vmax_resolution'][bin_index]['either_unresolved_count']==1
    bin_index=np.searchsorted(analyze.EDGES,np.log10(masses[4]),side='right')-1
    assert abs(stats['vmax'][bin_index]['q16_median_q84'][1]-10.)<1.e-12
    json.dumps(analyze.finite_json(result),allow_nan=False)
    return dict(partial_publication_bin_excluded=True,common_mass_floor_explicit=True,
                aperture_mass_not_mislabeled_so=True,either_unresolved_vmax_excluded_and_counted=True,
                known_resolved_vmax_shift_preserved=True)


def replay_publication(root):
    import struct

    sys.path.insert(0,str(ROOT/'replay'))
    helper=load_module('adapter_fixtures',ROOT/'replay/test_adapter.py')
    source=root/'snapshot';helper.fixture(source,16)
    source_files=[source/'PMcrd.0001.DAT',source/'PMcrs0.0001.DAT']
    spec=dict(case='fixture',nrow=4,ngrid=16,box_mpc_h=32.,epochs=[0],half_steps=False,
              cosmology=dict(Omega_m=.3),shape_repair_commit='fixture')
    simulation=dict(completed=True,spec=spec,outputs=file_manifest(source_files),snapshots=[dict(
        redshift=0,header=dict(step=1,scale_factor=1.),header_path=str(source_files[0]),data_paths=[str(source_files[1])])])
    original_root,original_work,original_run=replays.ROOT,replays.WORK,replays.run_native
    original_config,original_bin=replays.CONFIG,replays.BIN
    replays.ROOT=root/'reports';replays.ROOT.mkdir();replays.WORK=root
    # Tiny native fixtures are explicitly recorded with job_id=None, never a fake Slurm ID.
    def native(*args,**kwargs):return campaign.run_native(*args,**kwargs,allow_login=True)
    replays.run_native=native
    try:
        write_json(replays.ROOT/'fixture-simulation.json',simulation)
        replays.replay('fixture',1,epoch=0,finder='v3',analysis_ngrid=16)
        receipt=replays.ROOT/'replay-fixture-z0-ng16.json'
        v3_before=sha(receipt)
        replays.replay('fixture',1,epoch=0,finder='v3',analysis_ngrid=16)
        assert sha(receipt)==v3_before,'Completed v3 validation was rewritten'
        replays.replay('fixture',1,epoch=0,finder='pair',analysis_ngrid=16)
        assert json.loads(receipt.read_text())['pair_completed']
        before=sha(receipt)
        replays.replay('fixture',1,epoch=0,finder='v3',analysis_ngrid=16)
        replays.replay('fixture',1,epoch=0,finder='pair',analysis_ngrid=16)
        assert sha(receipt)==before,'Completed comparison was downgraded or rewritten'
        rejections={}
        def reject_resume(label):
            try:replays.replay('fixture',1,epoch=0,finder='pair',analysis_ngrid=16)
            except ValueError as error:rejections[label]=str(error)
            else:raise AssertionError(f'Changed completed replay was accepted: {label}')
            assert sha(receipt)==before,f'Rejected reuse rewrote published science: {label}'

        density=root/'replays/fixture-z0-ng16/density.bin'
        with density.open('r+b') as stream:
            original=stream.read(1);stream.seek(0);stream.write(bytes([original[0]^1]))
        try:reject_resume('retained_density')
        finally:
            with density.open('r+b') as stream:stream.write(original)

        # A new, independently valid simulation receipt with the same spec must
        # not inherit catalogues belonging to the previous particle realization.
        with source_files[1].open('r+b') as stream:
            original_particle=stream.read(4);stream.seek(0)
            stream.write(struct.pack('>f',struct.unpack('>f',original_particle)[0]+.125))
        try:
            simulation['outputs']=file_manifest(source_files)
            write_json(replays.ROOT/'fixture-simulation.json',simulation)
            reject_resume('new_valid_snapshot_same_spec')
        finally:
            with source_files[1].open('r+b') as stream:stream.write(original_particle)
            simulation['outputs']=file_manifest(source_files)
            write_json(replays.ROOT/'fixture-simulation.json',simulation)

        replays.CONFIG=original_config.replace('2.5e12','3.0e12')
        assert replays.CONFIG!=original_config
        try:reject_resume('requested_configuration')
        finally:replays.CONFIG=original_config
        config=root/'replays/fixture-z0-ng16/legacy/BDM.config'
        config.write_text(original_config+'! Changed frozen configuration\n')
        try:reject_resume('frozen_legacy_configuration')
        finally:config.write_text(original_config)

        # Redirect only this fixture's lookup; campaign executables stay intact.
        changed_bin=root/'changed-bin';changed_bin.mkdir()
        for variant in ['v3','legacy']:
            (changed_bin/f'PMP2replay.{variant}.exe').symlink_to(original_bin/f'PMP2replay.{variant}.exe')
        replays.BIN=changed_bin
        try:
            for variant in ['v3','legacy']:
                binary=changed_bin/f'PMP2replay.{variant}.exe';binary.unlink()
                binary.write_text(f'Changed {variant} binary: resume must reject before execution\n')
                try:reject_resume(f'{variant}_binary')
                finally:
                    binary.unlink();binary.symlink_to(original_bin/binary.name)
        finally:replays.BIN=original_bin
        replays.replay('fixture',1,epoch=0,finder='pair',analysis_ngrid=16)
        assert sha(receipt)==before,'Restored inputs did not preserve the completed comparison'
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
        stopped_before=sha(replays.ROOT/'replay-stopped-z0-ng16.json')
        replays.replay('stopped',1,epoch=0,finder='v3',analysis_ngrid=16)
        assert sha(replays.ROOT/'replay-stopped-z0-ng16.json')==stopped_before
        return dict(real_native_pair=True,completed_v3_preserved_on_resume=True,
                    completed_pair_preserved_on_both_resume_modes=True,
                    changed_retained_density_rejected=True,resume_identity_negative_controls=rejections,
                    restored_inputs_preserve_completed_pair=True,
                    legacy_exit_zero_partial_catalogue_preserves_v3=True,
                    valid_v3_resume_ignores_failed_legacy_stage=True)
    finally:
        replays.ROOT,replays.WORK,replays.run_native=original_root,original_work,original_run
        replays.CONFIG,replays.BIN=original_config,original_bin


def main():
    root=WORK/'driver-controls';root.mkdir(parents=True,exist_ok=False)
    result=dict(completed=False,started_at_utc=now(),matching=matching(root),statistics=statistics())
    result['replay']=replay_publication(root)
    result.update(completed=True,completed_at_utc=now())
    write_json(ROOT/'driver-controls.json',result)
    print('Membership matching and replay driver controls passed',flush=True)


if __name__=='__main__':main()
