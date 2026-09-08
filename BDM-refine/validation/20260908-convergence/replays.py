"""Legacy/v3 comparisons on identical particles and a shared analysis density."""
import json
from pathlib import Path
import struct

import numpy as np

from common import BIN, CHECKER_SHA, CONFIG, ROOT, WORK, checker, file_manifest, now, sha, verify_manifest, write_json
from campaign import fixed_text, link, run_native


def replay(name,threads,epoch=None,finder='pair',analysis_ngrid=2048):
    simulation=json.loads((ROOT/f'{name}-simulation.json').read_text())
    assert simulation['completed']
    spec=dict(simulation['spec'],evolution_ngrid=simulation['spec']['ngrid'],ngrid=analysis_ngrid)
    check=checker()
    for snapshot in simulation['snapshots']:
        z=snapshot['redshift']
        if epoch is not None and z!=epoch:continue
        tag=f'{name}-z{z}-ng{analysis_ngrid}'
        parent=WORK/'replays'/tag;parent.mkdir(parents=True,exist_ok=True)
        source=[Path(snapshot['header_path']),*map(Path,snapshot['data_paths'])]
        for p in source:
            verify_manifest({str(p):simulation['outputs'][str(p)]})
        density=parent/'density.bin'
        step=snapshot['header']['step']
        aggregate=ROOT/f'replay-{tag}.json'
        if aggregate.exists():
            saved=json.loads(aggregate.read_text())
            if saved.get('completed') and (saved.get('pair_completed') or finder=='v3'):
                assert saved['spec']==spec and saved['threads']==threads
                assert saved['membership_checker_sha256']==CHECKER_SHA
                verify_manifest(saved['science_files'])
                verify_manifest({saved['density']['path']:{key:saved['density'][key] for key in ['sha256','bytes']}})
                print(f'Verified completed {tag}; preserving its full comparison',flush=True)
                continue
        reports={};arrays={}
        density_evidence=None
        def publish():
            # Different content gets a distinct file: a failed pair cannot
            # invalidate the already published v3 arrays or validation receipt.
            scope='pair' if 'legacy' in arrays else 'v3'
            packed=parent/f'catalogues-{scope}.npz'
            staged=parent/f'catalogues-{scope}.tmp.npz'
            np.savez_compressed(staged,**arrays);staged.replace(packed)
            permanent=[packed]
            for variant in arrays:
                folder=parent/variant
                permanent.extend([folder/'finder.json',folder/'finder.log',folder/'finder.time',
                                  *list((folder/'CATALOGS').glob('Catshort*.DAT'))])
                if variant=='v3':permanent.extend([folder/'repair-members.bin',folder/'repair-members.index.npz'])
            report=dict(completed=True,scope=scope,requested_scope=finder,
                        pair_completed='legacy' in arrays,spec=spec,redshift=z,step=step,
                        membership_checker_sha256=CHECKER_SHA,threads=threads,science_files=file_manifest(permanent),
                        density=density_evidence,variants=reports,
                        catalogue_arrays=str(packed),catalogue_arrays_sha256=sha(packed),completed_at_utc=now())
            write_json(aggregate,report)
            if scope=='v3':write_json(ROOT/f'v3-validation-{tag}.json',report)
        for variant in ['v3','legacy']:
            if finder=='v3' and variant=='legacy':continue
            folder=parent/variant;folder.mkdir(exist_ok=True);(folder/'CATALOGS').mkdir(exist_ok=True)
            fixed_text(folder/'BDM.config',CONFIG)
            for p in source:link(folder/p.name,p)
            inputs=[*source,folder/'BDM.config']
            if variant=='legacy':inputs.append(density)
            def outputs():
                result=list((folder/'CATALOGS').glob('Catshort*.DAT'))
                if variant=='v3':result.extend([folder/'repair-members.bin',density])
                return result
            try:
                record=run_native(BIN/f'PMP2replay.{variant}.exe',folder,'finder',threads,inputs,outputs,
                                  arguments=[step,analysis_ngrid,threads,'write' if variant=='v3' else 'read',density],
                                  timeout_seconds=900 if variant=='legacy' else None)
            except Exception as error:
                if variant=='v3':raise
                reports['legacy']=dict(completed=False,error=repr(error),stage_receipt=str(folder/'finder.json'),
                                       interpretation='Legacy failure retained; audited v3 remains independently validated.')
                publish()
                print(f'Legacy failed at {tag}: {error}; preserving v3 and continuing',flush=True)
                continue
            if 'REPLAY COMPLETE' not in (folder/'finder.log').read_text():
                if variant=='v3':raise RuntimeError('Audited finder did not reach its completion marker')
                reports['legacy']=dict(completed=False,error='Missing completion marker despite process exit zero',
                                       stage_receipt=str(folder/'finder.json'))
                publish();continue
            assert density.stat().st_size==28+4*analysis_ngrid**3
            with density.open('rb') as stream:
                header=struct.unpack('>qqiff',stream.read(28))
            assert header[:3]==(analysis_ngrid,spec['nrow']**3,step)
            assert header[3]==snapshot['header']['scale_factor'] and header[4]==spec['box_mpc_h']
            try:
                cats=list((folder/'CATALOGS').glob('Catshort*.DAT'));assert len(cats)==1
                data,cat=check.load_catalogue(cats[0],spec,strict=variant=='v3')
            except Exception as error:
                if variant=='v3':raise
                reports['legacy']=dict(completed=False,error=repr(error),stage_receipt=str(folder/'finder.json'))
                publish();continue
            evidence=dict(catalogue=cat,stage_receipt=str(folder/'finder.json'),
                          elapsed_seconds=record['elapsed_seconds'],maxrss_kib=record['maxrss_kib'])
            if variant=='v3':
                evidence['diagnostics']=check.native_diagnostics(folder/'finder.log')
                evidence['membership']=check.check_memberships(folder/'repair-members.bin',spec,data)
                density_evidence=dict(path=str(density),**record['outputs'][str(density)],
                                      retained=True,shared_bit_identical_density='legacy' in arrays)
            else:
                assert record['identity']['inputs'][str(density)]['sha256']==density_evidence['sha256']
                assert sha(density)==density_evidence['sha256']
                density_evidence['shared_bit_identical_density']=True
            arrays[variant]=data;reports[variant]=evidence
            publish()
        print(f'Validated {tag}: '+', '.join(f'{v}={len(a)} haloes' for v,a in arrays.items()),flush=True)
