"""Independent review probes of integrated halo repair arithmetic.

Run with micromamba run -n cosemu python3 -B halo_review.py.
Does not modify production or earlier regression sources/results.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SOURCE = (REPO/'PMP2linker.f90').read_text()


def extract(name):
    return re.search(rf'^\s*(?:pure\s+)?(?:real\s*\*\s*8\s+)?(?:subroutine|function)\s+{name}\b.*?^\s*end\s+(?:subroutine|function)\s+{name}\b[^\n]*',
                     SOURCE, re.I|re.M|re.S).group()


def build(work, transform_program=None):
    structures = re.search(r'^module\s+Structures\b.*?^end module\s+Structures',
                           SOURCE, re.I|re.M|re.S).group()
    stub = '''module Tools
real :: Box=32.,AEXPN=.8
integer :: NGRID=128,NROW=32
integer*8 :: Nparticles=0
real,allocatable :: Xpar(:),Ypar(:),Zpar(:),VX(:),VY(:),VZ(:)
contains
real function seconds()
seconds=0.
end function
real function Memory(n)
integer*8 :: n
Memory=0.
end function
end module
module LinkerList
use Structures
use Tools
contains
'''
    names = ['GetHalo','BdmHaloGather','BdmHaloSortRadii','BdmHaloSortIds',
             'BdmHaloSphericalPotential','BdmHaloMembershipInit','List','Limits','EigenValues']
    for name in ['BdmParticlePosition', 'BdmParticleCoordinate']:
        if re.search(rf'\b(?:subroutine|function)\s+{name}\b', SOURCE, re.I):
            names.append(name)
    text = structures+'\n'+stub+'\n'.join(extract(n) for n in names)+'\nend module\n'
    program=(HERE/'halo_review_cases.f90').read_text()
    if transform_program is not None:
        program=transform_program(program)
    text += program
    (work/'source.f90').write_text(text)
    records = {}
    for mode, flags in [('checked',['-O0','-fcheck=all','-ffpe-trap=invalid,zero,overflow']),
                        ('optimized',['-O3'])]:
        cmd = ['gfortran',*flags,'-fopenmp','-ffree-line-length-none','source.f90','-o',mode]
        result = subprocess.run(cmd,cwd=work,text=True,capture_output=True)
        assert result.returncode == 0, result.stderr
        records[mode] = dict(command=cmd,stderr=result.stderr)
    return records


def sphere(n, radius):
    q=np.arange(n,dtype=float)
    z=1.-2.*(q+.5)/n
    phi=q*2.399963229728653
    return radius*np.c_[np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z]


def run_halo(work, name, phase, mode, cell=1., mass=1.e12, identity=False):
    data=np.asarray(phase,dtype=np.float32)
    with tempfile.TemporaryDirectory(prefix='bdm-halo-review-case-') as scratch:
        wd=Path(scratch)
        with (wd/'input.dat').open('w') as out:
            out.write(f'{len(data)} {cell} {mass:.17g} .3 200 .8 0 .2\n5 5 5\n')
            for j,row in enumerate(data):
                particle_id=j+1 if identity else 9000000000+len(data)-j
                out.write(' '.join(f'{v:.17g}' for v in row)+f' {particle_id}\n')
        result=subprocess.run([str(work/mode),'halo'],cwd=wd,capture_output=True,text=True,
                              env={**os.environ,'OMP_NUM_THREADS':'1'},timeout=15)
        assert result.returncode==0,(name,mode,result.stdout,result.stderr)
        value_line=next(l for l in result.stdout.splitlines() if l.strip().startswith('VALUES'))
        values=[float(v) for v in value_line.split()[1:]]
        id_line=next(l for l in result.stdout.splitlines() if l.strip().startswith('IDS'))
        ids=[int(v) for v in id_line.split()[1:]]
        assert np.isfinite(values).all(),(name,values)
        assert 'EMPTY_RECALL_PASS' in result.stdout
        assert ids==sorted(set(ids))
        assert all(1<=i<=len(data) for i in ids) if identity else all(i>2**31 for i in ids)
        if ids:
            assert np.isclose(values[1],len(ids)*float(np.float32(mass)),rtol=6.e-8)
        return dict(case=name,mode=mode,cell=cell,particle_count=len(data),values=values,ids=ids,
                    stdout=result.stdout,stderr=result.stderr)


def main():
    outputs=[]
    with tempfile.TemporaryDirectory(prefix='bdm-halo-review-build-') as scratch:
        work=Path(scratch)
        compilation=build(work)
        inner=sphere(10,.1)+5.
        outer=sphere(200,1.2)+5.
        two_crossing=np.c_[np.r_[inner,outer],np.zeros((210,3))]
        for mode in ['checked','optimized']:
            for cell in [1.,2.]:
                outputs.append(run_halo(work,'two_so_crossings',two_crossing,mode,cell))
            for n in [0,1,9,10]:
                phase=np.c_[sphere(n,.1)+5.,np.zeros((n,3))] if n else np.empty((0,6))
                outputs.append(run_halo(work,f'population_{n}',phase,mode))
            phase=np.c_[sphere(10,.1)+5.,sphere(10,1.e7)]
            phase[0,:3]=5.
            # Nine outer velocities arranged in cancelling triples; only the
            # exactly central, motionless row can remain after contaminants go.
            phase[0,3:]=0.
            phase[1:,3:]=np.tile([[1.e7,0.,0.],[-5.e6,1.e7,0.],[-5.e6,-1.e7,0.]],(3,1))
            outputs.append(run_halo(work,'one_central_survivor',phase,mode))
            for radii in [[],[0.],[.2],[0.,.2],[0.,0.],[.1,.1,.2,.2,.5]]:
                stdin=str(len(radii))+'\n'+' '.join(map(str,radii))+'\n'
                result=subprocess.run([str(work/mode),'potential'],input=stdin,text=True,capture_output=True)
                assert result.returncode==0,(radii,result.stderr)
                line=next(l for l in result.stdout.splitlines() if l.strip().startswith('POTENTIAL'))
                _,flag,energy=line.split()
                singular=len(radii)>1 and radii[1]==0
                assert (flag=='T')==singular
                if not singular:
                    expected=3.*4.*math.fsum(1/max(a,b) for i,a in enumerate(radii) for b in radii[i+1:])
                    assert math.isclose(float(energy),expected,rel_tol=2e-14,abs_tol=1e-14)
                outputs.append(dict(case='shell_potential',mode=mode,radii=radii,stdout=result.stdout))
            result=subprocess.run([str(work/mode),'heap_index'],input='1200000000 1100000000\n',text=True,capture_output=True)
            outputs.append(dict(case='heap_index_overflow',mode=mode,stdout=result.stdout,stderr=result.stderr))
    roots={str(n):(n*float(np.float32(1.e12))/(1.150e12*float(np.float32(.3))*200.))**(1/3)
           for n in [10,210]}
    report=dict(reviewed_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
                halo_repair_commit='8b85c9fca71e3e7e8a5cbc54a5e1d1b264c6107b',
                source_sha256=hashlib.sha256(SOURCE.encode()).hexdigest(),compilation=compilation,
                test_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in [Path(__file__),HERE/'halo_review_cases.f90']},
                exact_so_roots=roots,experiments=outputs)
    (HERE/'halo_review_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(exact_so_roots=roots,so_probes=[r for r in outputs if r['case']=='two_so_crossings'],
                         heap_probes=[r for r in outputs if r['case']=='heap_index_overflow']),indent=2))


if __name__=='__main__':main()
