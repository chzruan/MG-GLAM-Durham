"""Verify outermost SO convergence and large heap indices after review fixes.

Run with micromamba run -n cosemu python3 -B halo_followup.py.
Historical halo_review_results.json is read-only; new evidence has its own file.
"""
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import halo_review as review

HERE=Path(__file__).resolve().parent


def program_transform(program):
    blocks=[]
    for case,routine in [('heap_index','BdmHaloSortRadii'),('heap_ids_index','BdmHaloSortIds')]:
        # Compile the actual production local-index declaration, then exercise
        # its large-child arithmetic without allocating a billion-row array.
        declaration=re.search(r'^\s*integer(?:\*8)?\s*::\s*first,last,parent,child,n\s*$',
                              review.extract(routine),re.I|re.M).group().strip()
        blocks.append(f'''  if(trim(which)=='{case}')then
    block
      {declaration}
      read(*,*)n,parent
      child=2*parent
      write(*,*)'HEAP_CHILD',n,parent,child,child<=n
    end block
    stop
  endif
''')
    blocks.append('''  if(trim(which)=='sort_ids')then
    read(*,*)n
    allocate(OriginalParticleId(n))
    read(*,*)OriginalParticleId
    call BdmHaloSortIds(OriginalParticleId)
    write(*,*)'SORTED_IDS',OriginalParticleId
    stop
  endif
''')
    start=program.index("  if(trim(which)=='heap_index')then")
    end=program.index("  if(trim(which)=='potential')then",start)
    return program[:start]+''.join(blocks)+program[end:]


def exact_outer_root(phase,mass,cap):
    # An ascending interval oracle, independent of the production contraction.
    radius=np.linalg.norm(np.asarray(phase,dtype=np.float32)[:,:3].astype(float)-5.,axis=1)
    radius=np.sort(radius[radius<=cap])
    threshold=1.150e12*float(np.float32(.3))*200.
    valid=[]
    for i,r in enumerate(radius):
        count=i+1
        if count<10:continue
        root=math.cbrt(count*float(np.float32(mass))/threshold)
        outer=radius[i+1] if i+1<len(radius) else math.nextafter(cap,math.inf)
        if r<=root<outer and root<=cap:valid.append((root,count))
    return max(valid,default=(0.,0))


def check_bound_statistics(phase,result,mass):
    population=np.asarray(phase,dtype=np.float32).astype(float)[np.asarray(result['ids'])-1]
    values=result['values']
    if len(population)==0:return
    mass=float(np.float32(mass));a=float(np.float32(.8));om=float(np.float32(.3))
    offset=population[:,:3]-5.
    bulk=np.array([math.fsum(population[:,i+3])/len(population) for i in range(3)])
    hv=100*a*math.sqrt(om/a**3+1-om)
    velocity=population[:,3:]-bulk+hv*offset
    speed2=np.sum(velocity*velocity,axis=1)
    kinetic=.5*mass*math.fsum(speed2)
    radius=np.linalg.norm(offset,axis=1)
    potential=4.333e-9/a*mass*mass*math.fsum(
        1/max(radius[i],radius[j]) for i in range(len(radius)) for j in range(i+1,len(radius)))
    assert np.allclose(values[8:11],bulk,rtol=6e-8,atol=1e-9)
    assert math.isclose(values[4],kinetic,rel_tol=6e-8,abs_tol=1e-9)
    assert math.isclose(values[5],potential,rel_tol=6e-8,abs_tol=1e-9)
    assert math.isclose(values[11],math.sqrt(math.fsum(radius*radius)/len(radius)),rel_tol=6e-8)
    angular=np.mean(np.cross(offset,velocity),axis=0)
    spin=np.linalg.norm(angular)*a*math.sqrt(math.fsum(speed2)/len(radius))/(len(radius)*mass)*1.632e8
    assert math.isclose(values[13],spin,rel_tol=6e-8,abs_tol=1e-12)
    result['independent_bound_statistics']=['bulk','kinetic','distinct-pair shell energy','RMS radius','spin']


def main():
    records=[]
    with tempfile.TemporaryDirectory(prefix='bdm-halo-followup-') as scratch:
        work=Path(scratch)
        compilation=review.build(work,program_transform)
        inner=review.sphere(10,.1)+5.
        outer=review.sphere(200,1.2)+5.
        two_shell=np.c_[np.r_[inner,outer],np.zeros((210,3))]
        # Each outer row lies between the SO radii for q and q+1 particles.
        # A naive contraction discards only one row per pass; this fixture
        # forces the bounded contraction's sorted fallback to finish the job.
        mass=1.e12;threshold=1.150e12*float(np.float32(.3))*200.
        delayed=np.c_[review.sphere(80,.1)+5.,np.zeros((80,3))]
        for q in range(11,81):
            r=.5*(math.cbrt(q*float(np.float32(mass))/threshold)+math.cbrt((q+1)*float(np.float32(mass))/threshold))
            delayed[q-1,:3]=[5.+r,5.,5.]
        cold_position=review.sphere(64,.2)+5.
        cold_velocity=np.c_[-300.*(cold_position[:,1]-5.),300.*(cold_position[:,0]-5.),np.zeros(64)]
        cold_velocity += [30000.,-20000.,10000.]
        hot_position=review.sphere(16,.2)+5.
        hot_velocity=np.tile([[1.e6,0.,0.],[-1.e6,0.,0.]],(8,1))+[30000.,-20000.,10000.]
        contaminated=np.r_[np.c_[cold_position,cold_velocity],np.c_[hot_position,hot_velocity]]
        for mode in ['checked','optimized']:
            for cell in [.5,1.,2.,4.]:
                result=review.run_halo(work,'two_so_crossings_fixed',two_shell,mode,cell,mass,identity=True)
                expected,count=exact_outer_root(two_shell,mass,min(15*cell,math.nextafter(16.,0.)))
                assert count==210 and len(result['ids'])==count
                assert math.isclose(result['values'][3],expected,rel_tol=6e-8)
                result.update(expected_so=expected,expected_count=count)
                check_bound_statistics(two_shell,result,mass)
                records.append(result)
            result=review.run_halo(work,'slow_contraction_fallback',delayed,mode,1.,mass,identity=True)
            expected,count=exact_outer_root(delayed,mass,15.)
            assert count==10 and len(result['ids'])==10
            assert math.isclose(result['values'][3],expected,rel_tol=6e-8)
            records.append(result)
            result=review.run_halo(work,'bound_statistics_with_hot_contaminants',contaminated,mode,1.,mass,identity=True)
            assert result['ids']==list(range(1,65))
            check_bound_statistics(contaminated,result,mass)
            records.append(result)
            for n in [0,1,9,10]:
                phase=np.c_[review.sphere(n,.1)+5.,np.zeros((n,3))] if n else np.empty((0,6))
                result=review.run_halo(work,f'population_{n}',phase,mode,identity=True)
                assert len(result['ids'])==(10 if n==10 else 0)
                records.append(result)
            for case in ['heap_index','heap_ids_index']:
                result=subprocess.run([str(work/mode),case],input='1200000000 1100000000\n',text=True,capture_output=True)
                parts=result.stdout.split()
                assert result.returncode==0 and parts[-2:] == ['2200000000','F'],result.stdout
                records.append(dict(case=case,mode=mode,stdout=result.stdout))
            result=subprocess.run([str(work/mode),'sort_ids'],input='5\n9000000003 2 9000000001 1 9000000002\n',text=True,capture_output=True)
            assert [int(i) for i in result.stdout.split()[1:]] == [1,2,9000000001,9000000002,9000000003]
            records.append(dict(case='sort_int64_values',mode=mode,stdout=result.stdout))
    report=dict(source_sha256=hashlib.sha256(review.SOURCE.encode()).hexdigest(),
                original_review='93b5f90',compilation=compilation,experiments=records,passed=len(records),
                test_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in [Path(__file__),HERE/'halo_review.py',HERE/'halo_review_cases.f90']})
    (HERE/'halo_followup_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'{len(records)} halo review follow-up experiments passed.')


if __name__=='__main__':main()
