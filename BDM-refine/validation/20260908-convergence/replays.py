"""Legacy/v3 comparisons on identical particles and a shared analysis density."""
import json
from pathlib import Path
import struct

import numpy as np

from common import BIN, CHECKER_SHA, CONFIG, ROOT, WORK, checker, file_manifest, now, sha, verify_manifest, write_json
from campaign import fixed_text, link, run_native


def _verify_aggregate(saved,spec,threads,snapshot,source_manifest,parent,analysis_ngrid):
    """Reuse old aggregate receipts only when their native stage identities still apply."""
    if (saved['spec']!=spec or saved['threads']!=threads or
        saved['membership_checker_sha256']!=CHECKER_SHA or
        saved['redshift']!=snapshot['redshift'] or saved['step']!=snapshot['header']['step']):
        raise ValueError(f'Changed replay request: {parent}')
    verify_manifest(saved['science_files'])
    density=parent/'density.bin'
    if Path(saved['density']['path']).resolve()!=density.resolve():
        raise ValueError(f'Changed shared density path: {parent}')
    density_manifest={key:saved['density'][key] for key in ['sha256','bytes']}
    verify_manifest({str(density):density_manifest})
    verified_files=dict(saved['science_files'])
    verified_files[str(density.resolve())]=density_manifest
    variants=['v3','legacy'] if saved.get('pair_completed') else ['v3']
    for variant in variants:
        folder=parent/variant;receipt=folder/'finder.json';config=folder/'BDM.config'
        if (str(receipt.resolve()) not in saved['science_files'] or
            Path(saved['variants'][variant]['stage_receipt']).resolve()!=receipt.resolve()):
            raise ValueError(f'Unverified finder stage receipt: {receipt}')
        stage=json.loads(receipt.read_text())
        # The requested text matters even if the previously frozen file is intact.
        if not config.is_file() or config.read_text()!=CONFIG:
            raise ValueError(f'Changed requested finder configuration: {config}')
        inputs=dict(source_manifest);inputs.update(file_manifest([config]))
        if variant=='legacy':inputs[str(density.resolve())]=density_manifest
        identity=dict(binary_sha256=sha(BIN/f'PMP2replay.{variant}.exe'),threads=threads,stdin='',
                      arguments=list(map(str,[snapshot['header']['step'],analysis_ngrid,threads,
                                              'write' if variant=='v3' else 'read',density])),inputs=inputs)
        if not stage.get('completed') or stage['identity']!=identity:
            raise ValueError(f'Changed completed finder stage identity: {receipt}')
        # Reconcile already verified hashes without rereading the large particle
        # and density files for each variant. No aggregate format change is needed.
        for kind in ['outputs','evidence']:
            for path,expected in stage[kind].items():
                if verified_files.get(path)!=expected:
                    raise ValueError(f'Finder stage differs from aggregate evidence: {path}')
        if variant=='v3' and stage['outputs'].get(str(density.resolve()))!=density_manifest:
            raise ValueError(f'Shared density is not the verified v3 output: {density}')


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
        source_manifest={}
        for p in source:
            expected=simulation['outputs'][str(p)]
            verify_manifest({str(p):expected})
            source_manifest[str(p.resolve())]=expected
        density=parent/'density.bin'
        step=snapshot['header']['step']
        aggregate=ROOT/f'replay-{tag}.json'
        if aggregate.exists():
            saved=json.loads(aggregate.read_text())
            if saved.get('completed') and (saved.get('pair_completed') or finder=='v3'):
                _verify_aggregate(saved,spec,threads,snapshot,source_manifest,parent,analysis_ngrid)
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
            if scope=='v3':
                independent=ROOT/f'v3-validation-{tag}.json'
                if independent.exists():
                    prior=json.loads(independent.read_text())
                    # Pair upgrades and legacy failures must not change the
                    # independently published v3 receipt used by analysis.
                    for key in ['spec','redshift','step','membership_checker_sha256','threads','science_files']:
                        if prior[key]!=report[key]:raise ValueError('Changed published v3 evidence: '+key)
                    if prior['density']['sha256']!=report['density']['sha256']:
                        raise ValueError('Changed published v3 density')
                else:write_json(independent,report)
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
