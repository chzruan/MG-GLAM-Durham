"""Compile the actual recentering, unbinding and removal routines on 128 particles.

Only a read-only membership tap is inserted into GetHalo. No finder arithmetic
is reimplemented. The small Tools module supplies the simulation globals.
Run with micromamba run -n cosemu python3 source_reproduction.py SOURCE BUILD.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess


def extract(source, name):
    pattern = rf'^\s*subroutine\s+{name}\b.*?^\s*end\s+subroutine\s+{name}\b[^\n]*'
    match = re.search(pattern, source, re.I | re.M | re.S)
    if match is None:
        raise ValueError(f'Missing routine {name}')
    return match.group()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('source', type=Path)
    p.add_argument('build', type=Path)
    p.add_argument('--fixed', action='store_true')
    a = p.parse_args()
    a.build.mkdir(parents=True, exist_ok=True)
    source = a.source.read_text()
    structures = re.search(r'^Module\s+Structures\b.*?^end Module Structures',
                           source, re.I | re.M | re.S).group()
    structures = structures.replace('end Module Structures',
                                    'logical :: member(128,3)=.false.\nend Module Structures')
    routines = '\n'.join(extract(source, name) for name in
                        ['FindDistinctCandidates', 'GetHalo', 'RemoveDuplicates',
                         'Limits', 'EigenValues', 'ListMaxima'])
    if a.fixed:
        for name in ['MergeNumericalDuplicates']:
            routines += '\n' + extract(source, name)
        # Internal functions of RemoveDuplicates/MergeNumericalDuplicates are
        # included by extraction; the standalone module holds shared rules.
        shared = re.search(r'^module BdmDuplicateRules\b.*?^end module BdmDuplicateRules',
                           source, re.I | re.M | re.S).group()
    else:
        shared = ''
    old = 'if(ee <= 0.)Ncount = Ncount +1'
    assert old in routines
    routines = routines.replace(old,
        'if(ee <= 0.) member(jp,ip)=.true.\n                           ' + old)
    stub = '''module Tools
      real :: Box=32., AEXPN=0.8, tstart=0., tfinish=0.
      integer :: NGRID=128
      real :: Xpar(128),Ypar(128),Zpar(128),VX(128),VY(128),VZ(128)
      contains
      real function seconds()
      call cpu_time(seconds)
      end function
      end module
    '''
    program = '''
    program reproduce
      use LinkerList
      implicit none
      integer :: q,a,b,c,i,j,k,h,expected
      integer*8 :: ip
      real :: origin,saved(10,3)
      logical :: keep(3)
      Nmaxima=3; Np=128; MassOne=1.e10; Om0=0.3; Ovdens=200.
      dLogR=0.02; Rext=0.6; SlopeR=1.667; Cell=0.5
      Nmx=-1;Nmy=-1;Nmz=-1;Nbx=65;Nby=65;Nbz=65
      allocate(Mvir(3),Rvir(3),Mtotal(3),Xoff(3),xMaxx(3),yMaxx(3),zMaxx(3))
      allocate(VxMaxx(3),VyMaxx(3),VzMaxx(3),EpotM(3),EkinM(3),LambdaM(3))
      allocate(RadRms(3),VmaxM(3),RmaxM(3),Xax(3),Yax(3),Zax(3),Axba(3),Axca(3))
      allocate(Lst(128),Label(-1:65,-1:65,-1:65))
      Label=0;Lst=0;q=0
      do h=1,2
        origin=5.+7.*(h-1)
        do a=1,4
        do b=1,4
        do c=1,4
          q=q+1
          Xpar(q)=origin+(a-2.5)*0.0125
          Ypar(q)=origin+(b-2.5)*0.0125
          Zpar(q)=origin+(c-2.5)*0.0125
          VX(q)=100.; VY(q)=20.; VZ(q)=-10.
          i=ceiling(Xpar(q)/Cell)-1
          j=ceiling(Ypar(q)/Cell)-1
          k=ceiling(Zpar(q)/Cell)-1
          Lst(q)=Label(i,j,k);Label(i,j,k)=q
        end do
        end do
        end do
      end do
      Xoff=100.;xMaxx=[4.99,5.01,12.];yMaxx=[5.,5.,12.];zMaxx=[5.,5.,12.]
      call FindDistinctCandidates
      if(xMaxx(1)/=xMaxx(2))error stop 'recentring did not converge'
      do ip=1,3
        call GetHalo(xMaxx(ip),yMaxx(ip),zMaxx(ip),VxMaxx(ip),VyMaxx(ip),VzMaxx(ip),ip)
      end do
      if(count(member(:,1))/=64) error stop 'unexpected bound count'
      if(any(member(:,1).neqv.member(:,2)))error stop 'bound particle sets differ'
      if(any(member(:,1).and.member(:,3)))error stop 'independent set overlaps'
      if(Mvir(1)/=Mvir(2))error stop 'identical sets have different masses'
      print *, 'SAME_BOUND_PARTICLE_IDS=64; DISTINCT_OBJECT_IDS=64'
      saved(1,:)=xMaxx;saved(2,:)=yMaxx;saved(3,:)=zMaxx
      saved(4,:)=VxMaxx;saved(5,:)=VyMaxx;saved(6,:)=VzMaxx
      saved(7,:)=Mvir;saved(8,:)=Mtotal;saved(9,:)=Rvir;saved(10,:)=Mvir/MassOne
      call ListMaxima
      call RemoveDuplicates
      expected=EXPECTED_SURVIVORS
      if(count(Mvir>MassOne)/=expected)error stop 'wrong survivor count'
      if(Mvir(1)/=saved(7,1).or.Mvir(3)/=saved(7,3))error stop 'representative mass changed'
      if(any(xMaxx/=saved(1,:)).or.any(yMaxx/=saved(2,:)).or.any(zMaxx/=saved(3,:))) &
          error stop 'positions changed'
      if(any(VxMaxx/=saved(4,:)).or.any(VyMaxx/=saved(5,:)).or.any(VzMaxx/=saved(6,:))) &
          error stop 'velocities changed'
      if(any(Mtotal/=saved(8,:)).or.any(Rvir/=saved(9,:)))error stop 'aperture changed'
      print *, 'SURVIVORS=',count(Mvir>MassOne),'; SURVIVING_FIELDS_UNCHANGED'
    end program
    '''.replace('EXPECTED_SURVIVORS', '2' if a.fixed else '3')
    text = (shared + '\n' + structures + '\n' + stub + '\nmodule LinkerList\n'
            'use Structures\nuse Tools\n' + ('use BdmDuplicateRules\n' if a.fixed else '') +
            'contains\n' + routines + '\nend module\n' + program)
    (a.build/'reproduction.f90').write_text(text)
    cmd = ['gfortran', '-O0', '-g', '-fopenmp', '-fcheck=all',
           '-ffree-line-length-none', 'reproduction.f90', '-o', 'reproduce']
    subprocess.run(cmd, cwd=a.build, check=True, capture_output=True, text=True)
    results = []
    for threads in [1,2,4]:
        result = subprocess.run(['./reproduce'], cwd=a.build,
                   env={**os.environ, 'OMP_NUM_THREADS':str(threads)},
                   capture_output=True, text=True, check=True)
        results.append(dict(threads=threads, stdout=result.stdout, stderr=result.stderr))
        print(threads, result.stdout)
    (a.build/'result.json').write_text(json.dumps(dict(
        source=str(a.source.resolve()), source_sha256=hashlib.sha256(source.encode()).hexdigest(),
        compile=cmd, fixed=a.fixed, results=results), indent=2)+'\n')


if __name__ == '__main__':
    main()
