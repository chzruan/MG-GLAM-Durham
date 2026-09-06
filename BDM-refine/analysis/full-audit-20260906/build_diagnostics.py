"""Build labelled audit binaries without modifying the production source.

Invoke under the same Intel modules as build_full.sh, using cosemu Python.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
BUILD=ROOT/'work/full-build'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--bounds-only',action='store_true')
    args=parser.parse_args()
    original=(REPO/'PMP2linker.f90').read_text()
    instrumented=original.replace('end Module Structures',
        'integer*8, allocatable :: AuditOriginal(:)\nend Module Structures')
    instrumented=instrumented.replace('      Call AddBuffer\n',
        "      Call AddBuffer\n      open(92,file='audit-members.bin',access='stream',form='unformatted',status='replace')\n")
    instrumented=instrumented.replace('      Call RemoveDuplicates\n',
        """      Call RemoveDuplicates
      close(92)
      open(91,file='audit-survivors.bin',access='stream',form='unformatted',status='replace')
      write(91) Nmaxima,MassOne,Mvir,Rvir,xMaxx,yMaxx,zMaxx
      close(91)
      deallocate(AuditOriginal)
""")
    instrumented=instrumented.replace('       ALLOCATE(VXbb(iPartMax),VYbb(iPartMax),VZbb(iPartMax))',
        '''       ALLOCATE(VXbb(iPartMax),VYbb(iPartMax),VZbb(iPartMax))
       allocate(AuditOriginal(iPartMax))
       do ic=1,Np
         AuditOriginal(ic)=ic
       enddo''')
    instrumented=instrumented.replace('                      Xbb(ip) = x',
        '                      AuditOriginal(ip)=ic\n                      Xbb(ip) = x')
    instrumented=instrumented.replace('      Real*8 :: Tensor(3,3)',
        '      Real*8 :: Tensor(3,3)\n      integer*8,allocatable :: AuditIDs(:)\n      allocate(AuditIDs(Np))')
    instrumented=instrumented.replace('if(ee <= 0.)Ncount = Ncount +1',
        'if(ee <= 0.)AuditIDs(Ncount+1)=AuditOriginal(jp)\n                           if(ee <= 0.)Ncount = Ncount +1')
    instrumented=instrumented.replace('               Mbound   = Ncount*MassOne',
        '''!$OMP CRITICAL (audit_member_output)
               write(92) ip,int(Ncount,8),AuditIDs(:Ncount)
!$OMP END CRITICAL (audit_member_output)
               Mbound   = Ncount*MassOne''')
    assert instrumented.count('AuditOriginal(ip)=ic')==1
    assert instrumented.count('write(92) ip')==1
    # A single controlled algorithm change, to measure the late-config finding.
    configured=original.replace('      Call ReadParameters(ISTEP)\n','').replace('      Call SetParameters\n','')
    configured=configured.replace('      If(mDENSIT==1)Call DENSIT',
        '      Call ReadParameters(ISTEP)\n      Call SetParameters\n      If(mDENSIT==1)Call DENSIT')
    common=['-g','-traceback','-qopenmp','-march=core-avx2','-shared-intel',
            '-mcmodel=medium','-convert','big_endian']
    fast=['-O3','-ftz','-unroll','-mfma','-fp-model','fast=1']
    objects=['PMP2mod_tools','PMP2mod_fft5','PMP2mod_random','PMP2mod_density','PMP2mod_power',
             'PMP2mod_analyze','PMP2MG_subroutines','PMP2mod_MGbackground','PMP2MGsolver_fR',
             'PMP2extradof','PMP2MGsolver_DGP','PMP2MGsolver_sym','PMP2MGsolver_kmf','PMP2MGsolver_csf','PMP2bdm']
    report=json.loads((ROOT/'diagnostic-build.json').read_text()) if args.bounds_only else {}
    env={**os.environ,'LD_LIBRARY_PATH':os.environ['BDM_AUDIT_NATIVE_LIBS']}
    env.pop('LIBRARY_PATH',None)  # Avoid cosemu's incompatible OpenMP link library.
    variants=[
        ('members',instrumented,fast,'Read-only original-row membership and survivor taps'),
        ('configured',configured,fast,'Trial only: read BDM config and set scales before density/peak work'),
        ('checked',original,['-O0','-fp-model','precise','-check','all','-fpe0'],
         'Initial mixed-instrumentation experiment: check all enables MemorySanitizer; runtime-library warnings are not attributed to the finder')]
    if args.bounds_only:
        prefiltered=original.replace('             iMax = 1                   ! look for all 26 neighbors',
            '             if(FI(M1,M2,M3)<=Ovdens/3.)cycle\n             iMax = 1                   ! look for all 26 neighbors')
        assert prefiltered!=original
        variants=[('bounds',original,['-O0','-fp-model','precise','-check','bounds','-fpe0'],
                   'Original finder with bounds checks, precise arithmetic and -fpe0 compilation; other objects, including the entry point, retain production flags'),
                  ('prefiltered',prefiltered,fast,
                   'Performance trial only: skip below-threshold cells before 26-neighbour comparisons; inherited physics defects remain')]
    for variant,source,flags,meaning in variants:
        path=BUILD/f'finder-{variant}.f90';path.write_text(source)
        compile_cmd=['ifx',*common,*flags,'-c',str(path),'-o',f'finder-{variant}.o']
        link_cmd=['ifx',*common,*flags,'-o',f'PMP2BDM.{variant}.exe',
                  *[x+'.o' for x in objects],f'finder-{variant}.o']
        logs=[]
        for command in [compile_cmd,link_cmd]:
            p=subprocess.run(command,cwd=BUILD,env=env,text=True,capture_output=True)
            logs.append(dict(command=command,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))
            if p.returncode:
                raise RuntimeError(p.stderr)
        report[variant]=dict(meaning=meaning,source_sha256=hashlib.sha256(source.encode()).hexdigest(),
            binary_sha256=hashlib.sha256((BUILD/f'PMP2BDM.{variant}.exe').read_bytes()).hexdigest(),build=logs)
    (ROOT/'diagnostic-build.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
