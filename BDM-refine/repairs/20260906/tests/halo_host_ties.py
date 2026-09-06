"""Native equal-mass host reproductions and independent host-priority oracles.

micromamba run -n cosemu python3 -B halo_host_ties.py --compiler gfortran
Outputs go only to --output, and all temporary compilation/input files disappear.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
BASELINE='3740f1b4102ba8e7b85ee781823aab5c64bfdd62'


def extract(source,name):
    return re.search(rf'^\s*subroutine\s+{name}\b.*?^\s*end\s+subroutine\s+{name}\b[^\n]*',source,re.I|re.M|re.S).group()


def build(work,compiler):
    source=(REPO/'PMP2linker.f90').read_text()
    before=subprocess.check_output(['git','show',f'{BASELINE}:PMP2linker.f90'],cwd=REPO,text=True)
    modules=[re.search(rf'^module\s+{name}\b.*?^end module\s+{name}',source,re.I|re.M|re.S).group()
             for name in ['BdmDuplicateRules','Structures']]
    stub='''module Tools
real :: Box=32.,AEXPN=.8
contains
real function seconds()
call cpu_time(seconds)
end function
end module
module LinkerList
use Structures
use Tools
use BdmDuplicateRules
contains
'''
    compilation={}
    for variant in ['before','repaired']:
        chosen=before if variant=='before' else source
        routines=extract(chosen,'RemoveDuplicates')+'\n'+extract(source,'MergeNumericalDuplicates')+'\n'+extract(source,'ListMaxima')
        # The historical host pass calls Limits; the repaired pass computes
        # double query bounds directly. Include the unchanged routine for both.
        routines+='\n'+extract(source,'Limits')
        generated='\n'.join(modules)+'\n'+stub+routines+'\nend module\n'+(HERE/'halo_host_ties_cases.f90').read_text()
        filename=f'{variant}.f90';(work/filename).write_text(generated)
        for mode in ['checked','optimized']:
            if compiler=='gfortran':
                flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if mode=='checked' else ['-O3','-fno-fast-math']
                common=['-fopenmp','-ffree-line-length-none']
            else:
                flags=['-O0','-check','bounds','-fpe0','-fp-model','precise'] if mode=='checked' else ['-O3','-fp-model','precise']
                common=['-qopenmp','-extend-source']
            binary=f'{variant}-{mode}'
            command=[compiler,*flags,*common,filename,'-o',binary]
            env=os.environ.copy()
            if compiler=='ifx':
                env['LD_LIBRARY_PATH']=env['BDM_AUDIT_NATIVE_LIBS'];env.pop('LIBRARY_PATH',None)
            result=subprocess.run(command,cwd=work,env=env,capture_output=True,text=True)
            assert result.returncode==0,result.stderr
            compilation[binary]=dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
    return compilation,hashlib.sha256(source.encode()).hexdigest(),hashlib.sha256(before.encode()).hexdigest()


def record(index,position,radius=.1,count=20,mass_one=1.,ids=None):
    return dict(candidate=index,position=np.asarray(position,dtype=np.float32).astype(float).tolist(),
                radius=float(np.float32(radius)),mass=float(np.float32(count*mass_one)),
                ids=list(range(1000*index,1000*index+count)) if ids is None else list(ids))


def cases():
    fixture=json.loads((HERE/'halo_equal_mass_native.json').read_text())
    native=[dict(candidate=c['candidate'],position=c['properties'][:3],radius=c['properties'][8],
                 mass=c['properties'][6],ids=c['ids']) for c in fixture['candidates']]
    result=[dict(name='native_three_pairs',box=fixture['box'],mass_one=fixture['mass_one'],records=native)]
    for left,right in [(261,263),(673,695),(1086,1087)]:
        pair=[dict(r) for r in native if r['candidate'] in (left,right)]
        shift=np.asarray([.001,.001,.001])-np.asarray(pair[0]['position'])
        for r in pair:r['position']=np.mod(np.asarray(r['position'])+shift,128.).astype(np.float32).astype(float).tolist()
        result.append(dict(name=f'native_periodic_corner_{left}_{right}',box=128.,mass_one=fixture['mass_one'],records=pair))
    def add(name,records,**kwargs):result.append(dict(name=name,box=32.,mass_one=1.,records=records,**kwargs))
    add('disjoint_low_mass',[record(1,[5.,8.,8.],radius=.021024,count=64),record(2,[5.15,8.,8.],radius=.021024,count=64)])
    add('overlapping_outskirts',[record(1,[5.,8.,8.],radius=.1),record(2,[5.15,8.,8.],radius=.1)])
    add('equal_mass_unequal_radius_reverse_containment',[record(1,[5.,8.,8.],radius=.02),record(2,[5.15,8.,8.],radius=.2)])
    add('equal_mass_unequal_radius_forward_containment',[record(1,[5.,8.,8.],radius=.2),record(2,[5.15,8.,8.],radius=.02)])
    add('immutable_equal_mass_chain',[record(i,[1.+.15*(i-1),8.,8.],radius=.2) for i in [1,2,3]])
    add('immutable_unequal_mass_chain',[record(i,[1.+.15*(i-1),8.,8.],radius=.2,count=10*i) for i in [1,2,3]])
    add('mass_priority_over_index',[record(1,[5.,8.,8.],radius=.2,count=19),record(2,[5.1,8.,8.],radius=.2,count=20)])
    for name,x in [('strict_on',np.float32(1.25)),('strict_inside',np.nextafter(np.float32(1.25),np.float32(1.))),
                   ('strict_outside',np.nextafter(np.float32(1.25),np.float32(2.)))]:
        add(name,[record(1,[1.,8.,8.],radius=.25),record(2,[x,8.,8.],radius=.25)])
    add('periodic_float64_cutoff',[record(1,[.1,8.,8.],radius=.2000006,count=20),record(2,[31.9,8.,8.],radius=.1,count=19)])
    add('periodic_corner_containment',[record(1,[.1,.1,.1],radius=.4),record(2,[31.9,31.9,31.9],radius=.4)])
    add('shifted_query_boundary',[record(1,[28.,4.,4.],radius=4.100000858306885,count=20),
                                  record(2,[0.10000038892030716,4.,4.],radius=.1,count=19)])
    add('exact_set_fallback',[record(1,[5.,8.,8.],radius=.1),record(2,[7.,8.,8.],radius=.1,ids=range(1000,1020))])
    add('stable_index_with_reverse_input',[record(i,[5.+.01*(i%3),8.,8.],radius=.1) for i in [105,42,7]])
    rng=np.random.default_rng(620906)
    graph=[]
    for i in range(1,65):
        position=np.mod(rng.normal(0.,1.,3)+([0.,0.,0.] if i<33 else [16.,16.,16.]),32.)
        graph.append(record(i,position,radius=rng.uniform(.2,2.5),count=int(rng.integers(4,25))))
    add('independent_cluster_graph',graph)
    return result


def oracle(case):
    records=case['records'];box=float(np.float32(case['box']));mass_one=float(np.float32(case['mass_one']))
    removed=set();edges=[]
    for low in records:
        if low['mass']<=mass_one:continue
        for high in records:
            priority=high['mass']>low['mass'] or (high['mass']==low['mass'] and high['candidate']<low['candidate'])
            if not priority:continue
            displacement=[float(a)-float(b) for a,b in zip(low['position'],high['position'])]
            displacement=[d-box*round(d/box) for d in displacement]
            distance2=math.fsum(d*d for d in displacement)
            if distance2<high['radius']**2:
                removed.add(low['candidate']);edges.append([high['candidate'],low['candidate']])
    keep=[];sets=set()
    for r in sorted(records,key=lambda r:r['candidate']):
        if r['mass']<=mass_one or r['candidate'] in removed:continue
        key=tuple(r['ids'])
        if key in sets:continue
        sets.add(key);keep.append(r['candidate'])
    return keep,edges


def run(work,case,variant,mode,threads,compiler):
    records=case['records']
    lines=[f"{max(r['candidate'] for r in records)} {len(records)} {case['box']:.17g} 3.5 {case['mass_one']:.17g}"]
    for r in records:
        lines.extend([f"{r['candidate']} {r['mass']:.17g} {r['radius']:.17g} {len(r['ids'])}",
                      ' '.join(f'{x:.17g}' for x in r['position']),' '.join(map(str,r['ids']))])
    env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE','OMP_PROC_BIND':'false'}
    if compiler=='ifx':env['LD_LIBRARY_PATH']=env['BDM_AUDIT_NATIVE_LIBS']
    result=subprocess.run([str(work/f'{variant}-{mode}')],input='\n'.join(lines)+'\n',cwd=work,
                          env=env,capture_output=True,text=True,timeout=15)
    assert result.returncode==0,(case['name'],variant,mode,result.stdout,result.stderr)
    keep=[int(match) for match in re.findall(r'^KEEP\s+(\d+)',result.stdout,re.M)]
    return dict(case=case['name'],variant=variant,mode=mode,threads=threads,keep=keep)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--compiler',choices=['gfortran','ifx'],default='gfortran')
    parser.add_argument('--output',type=Path);args=parser.parse_args();fixtures=cases();records=[]
    with tempfile.TemporaryDirectory(prefix='bdm-host-ties-') as directory:
        work=Path(directory);compilation,source_sha,before_sha=build(work,args.compiler)
        before=[]
        for mode in ['checked','optimized']:
            for case in fixtures:
                expected,edges=oracle(case)
                for threads in [1,2,4,8]:
                    result=run(work,case,'repaired',mode,threads,args.compiler)
                    assert result['keep']==expected,(result,expected,edges)
                    records.append(result)
                if case['name'] in ['native_three_pairs','periodic_float64_cutoff','shifted_query_boundary']:
                    result=run(work,case,'before',mode,1,args.compiler)
                    assert result['keep']!=expected,(result,expected)
                    result.update(expected_after_repair=expected);before.append(result)
            print(args.compiler,mode,len(fixtures),'host scenarios passed at 1/2/4/8 threads',flush=True)
    native=fixtures[0]['records'];native_pairs=[]
    for a,b in zip(native[::2],native[1::2]):
        native_pairs.append(dict(candidates=[a['candidate'],b['candidate']],counts=[len(a['ids']),len(b['ids'])],
                                 shared=len(set(a['ids'])&set(b['ids'])),keepers=oracle(dict(records=[a,b],box=128.,mass_one=fixtures[0]['mass_one']))[0]))
    report=dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),compiler=args.compiler,
                source_sha256=source_sha,before_source_sha256=before_sha,compilation=compilation,
                test_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'halo_host_ties_cases.f90',HERE/'halo_equal_mass_native.json']},
                native_reproductions=native_pairs,before_failures=before,passed=len(records),results=records,
                oracle='Immutable all-pairs higher-bound-mass / lower-index priority plus exact-set fallback; strict double minimum-image containment')
    if args.output:args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(len(records),'host-priority regressions passed; historical failures reproduced')


if __name__=='__main__':main()
