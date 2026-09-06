"""Independent SO, binding, energy, and population regressions.

Run with micromamba run -n cosemu python3 -B halo_regression.py.
Temporary build/run products are removed automatically; --output records evidence.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

# Array oracles are small; do not create a BLAS thread pool on the login host.
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
SOURCE=REPO/'PMP2linker.f90'


def extract(source,name,kind='subroutine'):
    pattern=rf'^\s*(?:pure\s+)?(?:real(?:\*8)?\s+)?{kind}\s+{name}\b.*?^\s*end\s+{kind}\s+{name}\b[^\n]*'
    found=re.search(pattern,source,re.I|re.M|re.S)
    if not found: raise ValueError(name)
    return found.group()


def build(work):
    source=SOURCE.read_text()
    modules=[re.search(rf'^module\s+{name}\b.*?^end module\s+{name}',source,re.I|re.M|re.S).group()
             for name in ['BdmDuplicateRules','Structures']]
    # Branch development can precede the particle agent's interface commit.
    # On the integrated source these stubs are absent; production declarations
    # are used. No production GetHalo arithmetic is instrumented or changed.
    declarations=[]
    for name,declaration in [
        ('OriginalParticleId','integer*8,allocatable :: OriginalParticleId(:)'),
        ('HaloSearchRadius','real :: HaloSearchRadius=0.'),
        ('ParticleSearchRadius','real :: ParticleSearchRadius=0.')]:
        if not re.search(rf'\b{name}\b',modules[1],re.I):declarations.append(declaration)
    modules[1]=re.sub(r'end module structures','\n'.join(declarations)+'\nend module Structures',modules[1],flags=re.I)
    stub='''module Tools
real :: Box=32.,AEXPN=.8,Om=.3,OmL=.7,ASTEP=.004,hubble=.7
integer :: NGRID=128,NROW=32,ISTEP=1,Nrealization=1
integer*8 :: Nparticles=0
character*45 :: HEADER='Controlled halo repair regression'
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
           'List','Limits','EigenValues','BdmParticlePosition']
    generated='\n'.join(modules)+'\n'+stub+'\n'.join(extract(source,n) for n in names)
    generated+='\n'+extract(source,'BdmParticleCoordinate','function')+'\nend module\n'
    (work/'source.f90').write_text(generated+(HERE/'halo_cases.f90').read_text())
    builds={}
    for mode,flags in [('checked',['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']),
                       ('optimized',['-O3','-fno-fast-math'])]:
        command=['gfortran',*flags,'-fopenmp','-ffree-line-length-none','source.f90','-o',mode]
        result=subprocess.run(command,cwd=work,text=True,capture_output=True)
        builds[mode]=dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
        if result.returncode:raise RuntimeError(result.stderr)
    return builds,declarations


def expected_so(radius,mass,threshold,maximum):
    """Independent ascending scan of constant-mass intervals; outermost SO root."""
    edges=np.sort(radius[radius<=maximum])
    count=np.arange(1,len(edges)+1,dtype=np.float64)
    roots=np.cbrt(count*mass/threshold)
    upper=np.r_[edges[1:],maximum]
    valid=(count>=10)&(roots>=edges)&(roots<upper)
    return float(roots[valid][-1]) if np.any(valid) else 0.


def spherical_pair_potential(radius,mass):
    """Independent all-pairs sum for modest binding fixtures, no self term."""
    values=radius.tolist()
    result=[math.fsum(1/max(a,b) for j,b in enumerate(values) if i!=j)
            for i,a in enumerate(values)]
    return 4.333e-9*mass/float(np.float32(.8))*np.asarray(result)


def exact_geometric_potential(position,mass):
    # math.dist + fsum gives an independent compensated scalar-pair oracle;
    # it also avoids allocating an O(N^2) distance matrix on the login host.
    points=position.tolist()
    result=[math.fsum(1/math.dist(a,b) for j,b in enumerate(points) if i!=j)
            for i,a in enumerate(points)]
    return 4.333e-9*mass/float(np.float32(.8))*np.asarray(result)


def reference_unbind(phase,rso,mass):
    """Different algorithm: dense pair mask, recomputed from survivor rows."""
    radius=np.linalg.norm(phase[:,:3]-5.,axis=1)
    keep=radius<=rso
    hubble_a=100*np.sqrt(float(np.float32(.3))/float(np.float32(.8))**3+1-float(np.float32(.3)))*float(np.float32(.8))
    passes=0
    while np.any(keep):
        indices=np.flatnonzero(keep)
        survivors=phase[indices]
        velocity=survivors[:,3:]-survivors[:,3:].mean(axis=0)+hubble_a*(survivors[:,:3]-5.)
        energy=.5*np.sum(velocity**2,axis=1)-spherical_pair_potential(radius[indices],mass)
        passes+=1
        if np.all(energy<=0):break
        keep[indices[energy>0]]=False
    return keep,passes


def run(work,name,mode='checked',threads=1,argument=None):
    with tempfile.TemporaryDirectory(prefix='bdm-halo-case-') as run_dir:
        env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE','OMP_PROC_BIND':'false',
             'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
        command=[str(work/mode),name]+([] if argument is None else [str(argument)])
        completed=subprocess.run(command,cwd=run_dir,env=env,capture_output=True,text=True,timeout=30)
        assert completed.returncode==0,(name,mode,completed.stdout,completed.stderr)
        outputs=[]
        for line in completed.stdout.splitlines():
            if line.startswith('HALO '):outputs.append([float(v) for v in line.split()[1:]])
        measured=np.asarray(outputs)
        assert np.all(np.isfinite(measured)),(name,measured)
        phase=np.fromfile(Path(run_dir)/'phase.bin',dtype=np.float32).reshape(6,-1).T.astype(np.float64)
        ids=np.fromfile(Path(run_dir)/'members.bin',dtype=np.int64)
        # Original rows map to themselves; the production radial and linked
        # list orders differ from the ascending final identity order.
        selected=np.zeros(len(phase),dtype=bool)
        selected[ids-1]=True
        assert np.all(ids[:-1]<ids[1:]),ids
        assert len(np.unique(ids))==len(ids)
        row=measured[0]
        status=int(row[1]); mass=row[20]
        fields=['candidate','status','Mbound','Mtotal','Rvir','Ekin','Epot','Vmax','Rmax',
                'Vx','Vy','Vz','Rrms','Xoff','lambda','b_over_a','c_over_a','axis_x','axis_y','axis_z',
                'particle_mass','Rext','dLogR']
        result=dict(case=name,mode=mode,threads=threads,argument=argument,
                    properties=dict(zip(fields,row.tolist())),bound_count=len(ids),
                    bound_ids_sha256=hashlib.sha256(ids.tobytes()).hexdigest(),checks=[])
        radius=np.linalg.norm(phase[:,:3]-5.,axis=1)
        threshold=1.150e12*float(np.float32(.3))*200.
        search_max=16. if name.startswith('so') else 7.5
        so=expected_so(radius,mass,threshold,search_max)
        if name in ['truncated','aperture_truncated']:
            assert status&2 and row[2]==0
            result['checks'].append('explicit search-domain rejection')
            return result
        if name=='singular':
            assert status&16 and row[2]==0
            result['checks'].append('explicit coincident-centre rejection')
            return result
        aperture=so+.25*min(row[21]/(so/.25)**float(np.float32(.2)),.75)
        assert np.isclose(row[4],aperture,rtol=6e-8,atol=1e-10),(name,row[4],aperture)
        assert np.isclose(row[3],np.sum(radius<=aperture)*mass,rtol=6e-8),(name,row[3])
        assert np.isclose(row[2],len(ids)*mass,rtol=6e-8)
        result['so_radius_reference']=so
        result['checks'].extend(['independent sorted-particle SO crossing','aperture mass uses Rext aperture',
                                 'bound mass matches exact unique original membership'])
        if name.startswith('so'):
            assert abs(so-1.)<3.e-5
            assert abs(len(ids)-np.sum(radius<=so))==0
        else:
            expected,passes=reference_unbind(phase,so,mass)
            assert np.array_equal(selected,expected),(name,np.sum(selected),np.sum(expected))
            result['reference_unbinding_passes']=passes
            result['checks'].append('independent dense-pair spherical unbinding fixed point')
        if not len(ids):
            assert status&8
            assert row[5]==row[6]==row[7]==0.
            result['checks'].append('empty survivor set has finite documented outputs')
            return result
        survivors=phase[selected]
        bulk=survivors[:,3:].mean(axis=0)
        hubble_a=100*np.sqrt(float(np.float32(.3))/float(np.float32(.8))**3+1-float(np.float32(.3)))*float(np.float32(.8))
        velocity=survivors[:,3:]-bulk+hubble_a*(survivors[:,:3]-5.)
        ekin=.5*mass*np.sum(velocity**2)
        rrms=np.sqrt(np.mean(np.sum((survivors[:,:3]-5.)**2,axis=1)))
        assert np.allclose(row[9:12],bulk,rtol=6e-8,atol=1e-7)
        assert np.isclose(row[5],ekin,rtol=6e-8)
        assert np.isclose(row[12],rrms,rtol=6e-8)
        offset=survivors[:,:3]-5.
        xoff=np.linalg.norm(np.mean(offset,axis=0))/aperture
        angular=np.mean(np.cross(offset,velocity),axis=0)
        spin=np.linalg.norm(angular)*float(np.float32(.8))*np.sqrt(np.sum(velocity**2)/len(ids))/(len(ids)*mass)*1.632e8
        assert np.isclose(row[13],xoff,rtol=6e-8,atol=1e-15)
        assert np.isclose(row[14],spin,rtol=6e-8,atol=1e-12)
        result['checks'].append('bound-only drift, kinetic energy, RMS radius, centre offset, and spin')
        if not name.startswith('so'):
            potential=spherical_pair_potential(radius[selected],mass)
            energy=.5*np.sum(velocity**2,axis=1)-potential
            assert np.all(energy<=0.)
            epot=.5*mass*np.sum(potential)
            assert np.isclose(row[6],epot,rtol=6e-8),(name,row[6],epot)
            geometric=exact_geometric_potential(survivors[:,:3],mass)
            geometric_energy=.5*np.sum(velocity**2,axis=1)-geometric
            result.update(max_survivor_spherical_energy=float(energy.max()),
                          geometric_positive_energy_count=int(np.sum(geometric_energy>0)),
                          geometric_pair_energy=float(.5*mass*np.sum(geometric)),
                          spherical_over_geometric_potential_ratio=float(epot/(.5*mass*np.sum(geometric))))
            result['checks'].append('distinct-pair shell energy, including inner shell and excluding self')
            assert np.all(geometric_energy<=0.),(name,geometric_energy.max())
            if 'shell' in name:assert abs(result['spherical_over_geometric_potential_ratio']-1.)<.1
        sorted_bound=np.sort(radius[selected])
        positive=sorted_bound>0
        profile=np.arange(1,len(sorted_bound)+1)[positive]*mass/sorted_bound[positive]
        rmax=sorted_bound[positive][np.argmax(profile)]
        if rmax<=.025:
            assert status&4 and row[7]==row[8]==0.
            result['checks'].append('unresolved Vmax has finite zero sentinel')
        else:
            vmax=np.sqrt(4.333e-9*profile.max()/float(np.float32(.8)))/np.sqrt(1.-.025/rmax)
            assert np.isclose(row[7],vmax,rtol=6e-8)
            assert np.isclose(row[8],rmax,rtol=6e-8)
            result['checks'].append('Vmax uses its own saved bound-profile radius')
        if name=='parallel':
            assert np.all(measured[:,1:]==measured[0,1:]),measured
            result['checks'].append('parallel candidate outputs and memberships agree')
        return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='bdm-halo-build-') as directory:
        work=Path(directory)
        builds,stubs=build(work)
        results=[]
        for mode in ['checked','optimized']:
            for name in ['cold','central','compact','cold_shell','hot_shell','mixed_shell','bulk_shell',
                         'singular','truncated','aperture_truncated']:
                results.append(run(work,name,mode))
                print(mode,name,'passed',flush=True)
            for spacing in [.04,.02,.01,.005]:
                results.append(run(work,'so',mode,argument=spacing))
            results.append(run(work,'so_extended',mode,argument=.02))
            for threads in [1,2,4]:results.append(run(work,'parallel',mode,threads))
        report=dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),
                    source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                    test_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'halo_cases.f90']},
                    build_interface_stubs=stubs,compilation=builds,results=results,
                    passed=len(results),
                    scope='Discrete SO and iterative spherical Newtonian binding; exact non-spherical binding is not claimed')
        if args.output:args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(f'{len(results)} halo repair regressions passed')


if __name__=='__main__':main()
