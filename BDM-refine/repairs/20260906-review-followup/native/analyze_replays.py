"""Summarize verified replays; run via micromamba run -n cosemu python3 -B.

The full-catalogue stage uses a one-core shared allocation. Match by the
unchanged density-peak candidate identity on exactly the same FI, not by the
published row number, which changes when the mass/host cuts change.
"""
import argparse
import json
import os
from pathlib import Path
import re
import resource
import time

import numpy as np

import run_replays as r


def routine(source, name, kind='subroutine'):
    qualifiers=r'(?:pure\s+logical\s+)?' if kind=='function' else ''
    match=re.search(r'^\s*'+qualifiers+kind+r'\s+'+name+r'\b.*?^\s*end\s+'+kind+r'\s+'+name+r'\b',
                    source,re.I|re.M|re.S)
    assert match,name
    return match.group().strip()


def statistics(values):
    values=np.asarray(values,dtype=np.float64)
    assert np.isfinite(values).all() and len(values)>0
    return dict(count=len(values),mean=float(values.mean()),
                percentiles={str(p):float(np.percentile(values,p)) for p in [0,16,50,84,100]})


def compare(first,second,output,box):
    a=Path(first['membership']['retained_raw'])
    b=Path(second['membership']['retained_raw'])
    if not a.is_absolute():a=r.REPO/a
    if not b.is_absolute():b=r.REPO/b
    assert r.sha(a)==first['membership']['raw_sha256']
    assert r.sha(b)==second['membership']['raw_sha256']
    assert first['density']['sha256']==second['density']['sha256']
    assert first['config_sha256']==second['config_sha256']
    assert first['membership']['candidates']==second['membership']['candidates']
    with np.load(a.with_suffix('.index.npz')) as saved:
        x={key:saved[key] for key in saved.files}
    with np.load(b.with_suffix('.index.npz')) as saved:
        y={key:saved[key] for key in saved.files}
    for indexed in [x,y]:
        props=indexed['properties']
        assert np.all((props[:,17]>=0)&(props[:,17]<=props[:,16])&(props[:,16]<=1))
        length=np.linalg.norm(props[:,18:21],axis=1)
        assert np.all((length==0)|np.isclose(length,1,rtol=2.e-7,atol=0))
    common,ia,ib=np.intersect1d(x['candidates'],y['candidates'],assume_unique=True,return_indices=True)
    assert len(common)>0
    xp,yp=x['properties'][ia],y['properties'][ib]
    above=(xp[:,6]>=10**12.5)&(yp[:,6]>=10**12.5)
    fields={'Mbound':6,'Mtotal':7,'reported_aperture_radius':8,'Vmax':11,'Rrms':15}
    changes={}
    for name,column in fields.items():
        valid=above&(xp[:,column]>0)&(yp[:,column]>0)
        changes[name]=statistics(100*(yp[valid,column]/xp[valid,column]-1))
    drift=np.linalg.norm(yp[:,3:6]-xp[:,3:6],axis=1)
    position=yp[:,:3]-xp[:,:3]
    position-=box*np.rint(position/box)
    identical=0
    overlapping=[]
    with a.open('rb') as fa,b.open('rb') as fb:
        for i,j in zip(ia,ib):
            fa.seek(int(x['offsets'][i]));fb.seek(int(y['offsets'][j]))
            left=fa.read(8*int(x['counts'][i]));right=fb.read(8*int(y['counts'][j]))
            assert len(left)==8*x['counts'][i] and len(right)==8*y['counts'][j]
            identical+=int(left==right)
            aa=np.frombuffer(left,dtype='>i8');bb=np.frombuffer(right,dtype='>i8')
            overlap=np.intersect1d(aa,bb,assume_unique=True).size
            overlapping.append(overlap/max(len(aa),len(bb)))
    output.update(reference_rows=len(x['candidates']),v3_rows=len(y['candidates']),matched_candidates=len(common),
                  reference_only=len(x['candidates'])-len(common),v3_only=len(y['candidates'])-len(common),
                  count_change_percent=100*(len(y['candidates'])/len(x['candidates'])-1),
                  identical_bound_sets_for_matched_candidates=identical,
                  changed_bound_sets_for_matched_candidates=len(common)-identical,
                  matched_both_above_log10_mass_12p5=int(above.sum()),
                  percentage_changes_for_both_above_mass_cut=changes,
                  all_matched_position_shift_mpc_h=statistics(np.linalg.norm(position,axis=1)),
                  all_matched_bulk_velocity_shift_km_s=statistics(drift),
                  all_matched_membership_overlap_fraction=statistics(overlapping),
                  reference_catalogue=first['catalogue'],v3_catalogue=second['catalogue'],
                  reference_membership=first['membership'],v3_membership=second['membership'],
                  shared_density=first['density'])
    output['finite_ordered_shapes_and_normalized_axes']=True
    return x,y


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--pilot',action='store_true',help='Small parser/reduction controls only; no full catalogue work')
    args=parser.parse_args()
    if args.pilot:
        start=time.monotonic()
        assert statistics([1,2,3])['percentiles']['50']==2
        assert routine('subroutine Foo\n x=1\nend subroutine Foo\n','Foo').strip().endswith('Foo')
        r.write_json(r.HERE/'analysis-pilot.json',dict(completed=True,driver_sha256=r.sha(__file__),
                     elapsed_seconds=time.monotonic()-start,maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                     scope='Small reductions and verbatim routine extraction. Full streaming/checker memory pilot is host-control job11949318.'))
        return
    assert os.environ.get('SLURM_JOB_ID') and int(os.environ['SLURM_CPUS_PER_TASK'])==1
    assert not (r.HERE/'comparison.json').exists()
    started=time.monotonic()
    plan=json.loads((r.HERE/'plan.json').read_text())
    build=json.loads((r.HERE/'build.json').read_text())
    groups={n:json.loads((r.HERE/f'results-t{n}.json').read_text()) for n in [32,64]}
    for group in groups.values():
        assert group['completed'] and group['optimization_byte_identical']
        assert group['input_snapshots_unchanged'] and group['input_density_tapes_unchanged']
        assert group['plan_sha256']==r.sha(r.HERE/'plan.json')
        assert group['build_sha256']==r.sha(r.HERE/'build.json')
        for stage in group['stages']:
            assert stage['completed'] and stage['particle_state_and_workspace_restored']
            assert stage['repeated_catalogue_byte_identical']
            for path,expected in stage['outputs_sha256'].items():assert r.sha(path)==expected,path
    stages={n:{s['case']['label']:s for s in group['stages']} for n,group in groups.items()}
    thread_control={}
    for variant in ['reference','optimized-v2','v3']:
        first,second=stages[32]['z0-'+variant],stages[64]['z0-'+variant]
        assert first['density']['sha256']==second['density']['sha256']
        assert first['catalogue']['sha256']==second['catalogue']['sha256']
        assert first['membership']['raw_sha256']==second['membership']['raw_sha256']
        assert first['unbinding']['sha256']==second['unbinding']['sha256']
        thread_control[variant]=dict(catalogue_byte_identical=True,membership_byte_identical=True,
                                    unbinding_diagnostics_byte_identical=True)
    # Source proof of candidate alignment: FindMaxima is unchanged and each
    # physical pair uses identical config, density and candidate count.
    source_paths={key:Path(value['source_path']) for key,value in build['variants'].items()}
    for key,path in source_paths.items():
        assert r.sha(path)==build['variants'][key]['source_sha256']
    for name,kind in [('FindMaxima','subroutine'),('SetOverdensity','subroutine'),('IsDensityMaximum','function')]:
        definitions={key:routine(path.read_text(),name,kind) for key,path in source_paths.items()}
        assert len(set(definitions.values()))==1,name
    result=dict(completed=False,started_at_utc=r.now(),job_id=os.environ['SLURM_JOB_ID'],
                plan_sha256=r.sha(r.HERE/'plan.json'),build_sha256=r.sha(r.HERE/'build.json'),
                driver_sha256=r.sha(__file__),fixed_density_32_64_thread_control=thread_control,
                interpretation='V2 and V3 matched by the identical initial density-peak candidate ID. '
                  'Unmatched rows are unmatched published candidate IDs, not a claim of new or lost physical objects. '
                  'Both matched masses must exceed 10^12.5 Msun/h for property percentage summaries. '
                  'Each property additionally requires both values >0 to exclude unresolved sentinels. '
                  'Membership overlap is |A intersection B| / max(|A|,|B|), not Jaccard. '
                  'The reported radius retains Rext and is not the unextended SO radius.',
                performance={},epochs={},unbinding={})
    arrays={}
    for n in [32,64]:
        first,second=stages[n]['z0-reference'],stages[n]['z0-optimized-v2']
        result['performance'][str(n)]={key:dict(reference_seconds=first['stage_seconds'][key],
            optimized_v2_seconds=second['stage_seconds'][key],
            median_speedup=float(np.median(first['stage_seconds'][key])/np.median(second['stage_seconds'][key])))
            for key in ['List','AddBuffer','ParametersDistinct']}
        result['performance'][str(n)]['first_pass_finder_without_tapes']=dict(
            reference_seconds=first['finder_seconds'][0],optimized_v2_seconds=second['finder_seconds'][0],
            speedup=first['finder_seconds'][0]/second['finder_seconds'][0])
    for z in [2,1,0]:
        first,second=stages[64][f'z{z}-reference'],stages[64][f'z{z}-v3']
        epoch={}
        left,right=compare(first,second,epoch,plan['spec']['box_mpc_h'])
        result['epochs'][str(z)]=epoch
        result['unbinding'][str(z)]={key:s['unbinding'] for key,s in [('v2',first),('v3',second)]}
        for label,stage,index in [('v2',first,left),('v3',second,right)]:
            arrays[f'{label}_z{z}']=np.loadtxt(stage['catalogue']['path'],skiprows=8,ndmin=2)
            arrays[f'{label}_candidates_z{z}']=index['candidates']
    np.savez_compressed(r.HERE/'comparison-catalogues.npz',**arrays)
    result.update(completed=True,finished_at_utc=r.now(),elapsed_seconds=time.monotonic()-started,
                  maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  comparison_catalogues_sha256=r.sha(r.HERE/'comparison-catalogues.npz'))
    r.write_json(r.HERE/'comparison.json',result)
    print('All fixed-FI thread, computational identity and physical comparison controls passed',flush=True)


if __name__=='__main__':
    main()
