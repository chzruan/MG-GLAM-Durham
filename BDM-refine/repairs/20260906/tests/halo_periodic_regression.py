"""Check original-row image arithmetic against independent periodic oracles.

micromamba run -n cosemu python3 -B halo_periodic_regression.py --compiler gfortran
For ifx, initialize the Intel runtime and export BDM_AUDIT_NATIVE_LIBS before
entering cosemu, as recorded in the repair validation scripts.
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
SOURCE=REPO/'PMP2linker.f90'
COVERAGE_CASES=('coverage','parallel_list','half_box')


def extract(source,name,kind='subroutine'):
    pattern=rf'^\s*(?:pure\s+)?(?:real\*8\s+)?{kind}\s+{name}\b.*?^\s*end\s+{kind}\s+{name}\b[^\n]*'
    found=re.search(pattern,source,re.I|re.M|re.S)
    if not found:raise ValueError(name)
    return found.group()


def build(work,compiler):
    source=SOURCE.read_text()
    modules=[re.search(rf'^module\s+{name}\b.*?^end module\s+{name}',source,re.I|re.M|re.S).group()
             for name in ['BdmDuplicateRules','Structures']]
    stub='''module Tools
real :: Box=32.,AEXPN=.8,Om=.3,OmL=.7,ASTEP=.004,hubble=.7
integer :: NGRID=128,NROW=32,ISTEP=1,Nrealization=1
integer*8 :: Nparticles=0
character*45 :: HEADER='Periodic image precision regression'
real,allocatable :: FI(:,:,:),Xpar(:),Ypar(:),Zpar(:),VX(:),VY(:),VZ(:)
contains
real function seconds()
use omp_lib, only: omp_get_wtime
seconds=real(omp_get_wtime())
end function
real function Memory(n)
integer*8 :: n
Memory=0.
end function
end module
module LinkerList
use Structures
use Tools
use BdmDuplicateRules
contains
'''
    names=['GetHalo','BdmHaloGather','BdmHaloSortRadii','BdmHaloSortIds',
           'BdmHaloSphericalPotential','BdmHaloMembershipInit','ParametersDistinct',
           'List','Limits','EigenValues','FindDistinctCandidates','BdmParticlePosition',
           'AddBuffer','PrepareParticleSearch','SizeList']
    generated='\n'.join(modules)+'\n'+stub+'\n'.join(extract(source,n) for n in names)
    generated+='\n'+extract(source,'BdmParticleCoordinate','function')+'\nend module\n'
    (work/'source.f90').write_text(generated+(HERE/'halo_periodic_cases.f90').read_text())
    builds={}
    if compiler=='gfortran':
        modes={'checked':['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'],
               'optimized':['-O3','-fno-fast-math']}
        common=['-fopenmp','-ffree-line-length-none']
    else:
        modes={'checked':['-O0','-check','bounds','-fpe0','-fp-model','precise'],
               'optimized':['-O3','-fp-model','precise']}
        common=['-qopenmp']
    for mode,flags in modes.items():
        command=[compiler,*flags,*common,'source.f90','-o',mode]
        env=os.environ.copy()
        if compiler=='ifx':
            env['LD_LIBRARY_PATH']=env['BDM_AUDIT_NATIVE_LIBS']
            env.pop('LIBRARY_PATH',None)
        result=subprocess.run(command,cwd=work,text=True,capture_output=True,env=env)
        builds[mode]=dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
        if result.returncode:raise RuntimeError(result.stderr)
    return builds


def coverage_inputs():
    # Fine offsets straddle list-cell boundaries and both box faces. These are
    # stored float32 inputs; the oracle never treats their decimal labels as exact.
    points=[]
    for x in [0.,.1,.5,np.nextafter(np.float32(.5),np.float32(0)),np.nextafter(np.float32(.5),np.float32(1)),31.9]:
        for y in [0.,.1,16.,31.9]:
            for z in [0.,.1,16.,31.9]:points.append([x,y,z])
    for q in range(128):points.append([((q*13)%127)*32/127,((q*31)%127)*32/127,((q*53)%127)*32/127])
    phase=np.zeros((len(points),6),dtype=np.float32);phase[:,:3]=points
    queries=[];radii=[]
    centers=np.asarray([[31.9,16,16],[.1,16,16],[31.9,31.9,31.9],[0,0,0],[.5,.5,.5],[16,16,16]],dtype=np.float32)
    for centre in centers:
        for radius in [.2,.5,1.,7.5]:queries.append(centre);radii.append(radius)
        for point in phase[::7,:3]:
            displacement=[float(a)-float(b) for a,b in zip(point,centre)]
            displacement=[d-32*round(d/32) for d in displacement]
            radius=math.sqrt(math.fsum(d*d for d in displacement))
            if 0<radius<7.5:
                for ratio in [1.-2.e-12,1.,1.+2.e-12]:queries.append(centre);radii.append(radius*ratio)
    return phase,np.asarray(queries,dtype=np.float32),np.asarray(radii,dtype=np.float64)


def cloud_inputs(corner=False):
    points=[]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                # Dyadic positions permit exact translation between the face,
                # corner, and interior; the negative mean must wrap canonically.
                offset=np.asarray([(a-2)/128,(b-2)/128,(c-2)/128])
                centre=np.asarray([0.,0.,0.] if corner else [0.,16.,16.])
                points.append(np.mod(centre+offset,32.))
    phase=np.zeros((64,6),dtype=np.float32);phase[:,:3]=points
    phase[:,3:]=[7.25,-3.5,11.]
    query=np.asarray([[0.,0.,0.] if corner else [0.,16.,16.]],dtype=np.float32)
    return phase,query,np.asarray([.2])


def run(work,compiler,mode,threads,name):
    if name in COVERAGE_CASES:phase,queries,radii=coverage_inputs()
    else:phase,queries,radii=cloud_inputs(corner=name=='corner')
    if name=='parallel_list':
        phase=np.tile(phase,(64,1)) # exceed the production 10,000-row serial shortcut
        queries=np.asarray([[31.9,16.,16.],[31.9,31.9,31.9],[0.,0.,0.],[.5,.5,.5],[16.,16.,16.]],dtype=np.float32)
        radii=np.asarray([.2,.4,.5,7.5,7.5])
    if name=='half_box':
        queries=np.asarray([[0.,0.,0.],[.1,16.,16.],[31.9,31.9,31.9],[16.,16.,16.]],dtype=np.float32)
        radii=np.full(len(queries),float(np.nextafter(np.float32(16.),np.float32(0.))))
    with tempfile.TemporaryDirectory(prefix='bdm-periodic-case-') as directory:
        directory=Path(directory)
        raw=np.asarray([len(phase),len(queries)],dtype=np.int32).tobytes()
        raw+=phase.T.copy().tobytes()+queries.copy().tobytes()+radii.tobytes()
        (directory/'input.bin').write_bytes(raw)
        env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE','OMP_PROC_BIND':'false'}
        if compiler=='ifx':env['LD_LIBRARY_PATH']=env['BDM_AUDIT_NATIVE_LIBS']
        action='coverage' if name in COVERAGE_CASES else 'halo'
        if name=='half_box':action='half_box'
        result=subprocess.run([str(work/mode),action],cwd=directory,env=env,capture_output=True,text=True,timeout=30)
        assert result.returncode==0,(name,result.stdout,result.stderr)
        raw=(directory/'images.bin').read_bytes();count=int(np.frombuffer(raw[:8],dtype=np.int64)[0]);offset=8
        stored=np.frombuffer(raw[offset:offset+count*12],dtype=np.float32).reshape(3,count).T.astype(float);offset+=count*12
        identity=np.frombuffer(raw[offset:offset+count*8],dtype=np.int64);offset+=count*8
        exact=np.frombuffer(raw[offset:],dtype=np.float64).reshape(count,3)
        original=phase[identity-1,:3].astype(float)
        shifts=np.rint((stored-original)/32.)
        expected=original+shifts*32.
        assert np.array_equal(exact,expected)
        extra_rounding=float(np.max(np.abs(stored-expected)))
        raw=(directory/'gathers.bin').read_bytes();offset=0
        for centre,radius in zip(queries,radii):
            n=int(np.frombuffer(raw[offset:offset+8],dtype=np.int64)[0]);offset+=8
            ids=np.frombuffer(raw[offset:offset+n*8],dtype=np.int64);offset+=n*8
            distance=np.frombuffer(raw[offset:offset+n*8],dtype=np.float64);offset+=n*8
            assert len(set(ids))==n,'one physical row counted through multiple images'
            displacement=phase[:,:3].astype(float)-centre.astype(float)
            displacement-=32.*np.rint(displacement/32.)
            distance2=displacement[:,0]**2+displacement[:,1]**2+displacement[:,2]**2
            expected_ids=np.flatnonzero(distance2<=radius**2)+1
            assert np.array_equal(np.sort(ids),expected_ids),(name,centre,radius,ids,expected_ids)
            expected_distance=np.sqrt(distance2[ids-1])
            assert np.allclose(distance,expected_distance,rtol=2e-15,atol=1e-14)
        checks=['exact image coordinates from stored original rows','original-ID minimum-image membership oracle',
                'face/corner/list-boundary cutoffs','no duplicate images in a supported sphere']
        record=dict(compiler=compiler,mode=mode,threads=threads,case=name,particle_count=len(phase),
                    query_count=len(queries),image_count=count,max_extra_float32_ghost_error=extra_rounding,checks=checks)
        if name not in COVERAGE_CASES:
            centres=np.fromfile(directory/'centres.bin',dtype=np.float32).reshape(7,len(queries)).T.astype(float)
            local=phase[:,:3].astype(float)-queries[0].astype(float)
            local-=32.*np.rint(local/32.)
            expected_centre=np.mod(queries[0].astype(float)+np.mean(local,axis=0),32.).astype(np.float32)
            assert np.array_equal(centres[0,:3],expected_centre)
            assert np.all((centres[:,:3]>=0)&(centres[:,:3]<32.))
            raw=(directory/'haloes.bin').read_bytes()
            values=np.frombuffer(raw[:24],dtype=np.float32).astype(float)
            status=int(np.frombuffer(raw[24:28],dtype=np.int32)[0])
            n=int(np.frombuffer(raw[28:36],dtype=np.int64)[0])
            ids=np.frombuffer(raw[36:],dtype=np.int64)
            assert n==len(phase) and np.array_equal(ids,np.arange(1,len(phase)+1))
            local=phase[:,:3].astype(float)-centres[0,:3]
            local-=32.*np.rint(local/32.)
            a=float(np.float32(.8));om=float(np.float32(.3));hubble_a=100*a*math.sqrt(om/a**3+1-om)
            kinetic=.5*1.e10*math.fsum(float(hubble_a**2)*math.fsum(float(v)**2 for v in point) for point in local)
            rrms=math.sqrt(math.fsum(math.fsum(float(v)**2 for v in point) for point in local)/n)
            assert np.isclose(values[3],kinetic,rtol=6e-8)
            assert np.isclose(values[5],rrms,rtol=6e-8)
            record.update(canonical_centre=centres[0,:3].tolist(),halo_properties=values.tolist(),halo_status=status)
            checks.extend(['negative mean wraps into canonical domain on every centring pass',
                           'periodic GetHalo membership, Hubble kinetic energy and RMS radius'])
        return record


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--compiler',choices=['gfortran','ifx'],default='gfortran')
    parser.add_argument('--output',type=Path);args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='bdm-periodic-build-') as directory:
        work=Path(directory);builds=build(work,args.compiler);results=[]
        for mode in ['checked','optimized']:
            for threads in [1,2,4]:
                for name in ['coverage','parallel_list','half_box','face','corner']:
                    results.append(run(work,args.compiler,mode,threads,name))
                    print(args.compiler,mode,threads,name,'passed',flush=True)
        # Exact dyadic translations produce identical physical halo properties.
        cloud=[r for r in results if r['case'] not in COVERAGE_CASES]
        assert all(r['halo_properties']==cloud[0]['halo_properties'] for r in cloud)
        report=dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),
                    source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                    test_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'halo_periodic_cases.f90']},
                    compilation=builds,passed=len(results),results=results,
                    translation_check='All face/corner/thread/compiler-mode physical properties agree exactly within this compiler')
        if args.output:args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(len(results),'periodic precision regressions passed')


if __name__=='__main__':main()
