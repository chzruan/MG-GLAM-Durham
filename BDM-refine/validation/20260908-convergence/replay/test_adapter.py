"""Small native integration and coordinate-unit controls; no Slurm required.

Execute using micromamba run -n cosemu python3 -B. Scratch is removed only
after the complete JSON receipt has been written successfully.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile

import numpy as np

import build_adapter as builder

GRID_TEST = '''program GridTest
 use Tools
 use ConvergenceGrid
 implicit none
 integer :: target,j,unit
 character(32) :: arg
 call get_command_argument(1,arg)
 read(arg,*)NGRID
 call get_command_argument(2,arg)
 read(arg,*)target
 NROW=2;Nparticles=8;Box=256.;AEXPN=0.5
 allocate(Xpar(8),Ypar(8),Zpar(8),VX(8),VY(8),VZ(8))
 do j=1,8
   read(*,*)Xpar(j),Ypar(j),Zpar(j),VX(j),VY(j),VZ(j)
 enddo
 call ChangeAnalysisGrid(target)
 open(newunit=unit,file='converted.bin',form='unformatted',access='stream',status='new')
 do j=1,8
   write(unit)Xpar(j),Ypar(j),Zpar(j),VX(j),VY(j),VZ(j)
 enddo
 close(unit)
end program
'''


def fixture(path, source):
    path.mkdir()
    (path/'CATALOGS').mkdir()
    (path/'BDM.config').write_text('! Shared requested settings; retain legacy first-call behavior\n'
                                 'iVirial = 1 ! virial definition\n'
                                 'MassMin = 1 ! retain the controlled halo\n'
                                 'Rext = 0.0 ! no aperture extension\n')
    header = b'Convergence replay control'.ljust(45) + struct.pack(
        '>4fif6f4i3fq100f', 1., .02, 1., .004, 1, 1., *([0.] * 6),
        4, source, 1, 1, .3, .7, .7, 64, *([0.] * 99 + [32.]))
    marker = struct.pack('>i', len(header))
    (path/'PMcrd.0001.DAT').write_bytes(marker+header+marker)
    common = np.array([[6.+.5*(i-1.5), 6.+.5*(j-1.5), 6.+.5*(k-1.5),
                        .0001*(i-1.5), .0001*(j-1.5), .0001*(k-1.5)]
                       for k in range(4) for j in range(4) for i in range(4)], dtype=np.float32)
    native = common.copy()
    native[:,:3] = 1.+(common[:,:3]-1.)*(source/16.)
    native[:,3:] = common[:,3:]*(source/16.)
    with (path/'PMcrs0.0001.DAT').open('xb') as stream:
        stream.truncate(6*1024**2*4)
        for axis in range(6):
            stream.seek(axis*1024**2*4)
            stream.write(native[:,axis].astype('>f4').tobytes())
    return common


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    build = json.loads((args.build_dir/'build.json').read_text())
    assert build['completed']
    common_dir = Path(build['common_dir'])
    scratch = Path(tempfile.mkdtemp(prefix='bdm-common-grid-controls-'))
    report = dict(completed=False, script_sha256=builder.sha(__file__),
                  build_sha256=builder.sha(args.build_dir/'build.json'), commands=[], grid_controls=[],
                  fixture_controls=[], scratch=str(scratch))
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    if env.get('BDM_AUDIT_NATIVE_LIBS'):
        env['LD_LIBRARY_PATH'] = env['BDM_AUDIT_NATIVE_LIBS']
    env.pop('LIBRARY_PATH', None)

    def run(command, cwd, stdin=None, success=True):
        process = subprocess.run(list(map(str,command)), cwd=cwd, env=env, input=stdin,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=180)
        report['commands'].append(dict(command=list(map(str,command)), returncode=process.returncode,
                                      output=process.stdout))
        if success:
            assert process.returncode == 0, process.stdout
        else:
            assert process.returncode != 0, 'Expected a failed control'
        return process

    try:
        # Exercise real production Tools plus the actual adapter with tiny arrays.
        for name in ['tools.mod']:
            shutil.copy2(common_dir/name,scratch/name)
        shutil.copy2(builder.HERE/'common_grid.f90',scratch/'common_grid.f90')
        (scratch/'grid_test.f90').write_text(GRID_TEST)
        run([build['compiler'],*builder.FLAGS,'-o','grid-test.exe',common_dir/'PMP2mod_tools.o',
             'common_grid.f90','grid_test.f90'],scratch)
        for source in [8,16,32,1024,2048,4096]:
            for target in [source,source//2,source*2]:
                # Include origin, upper image, adjacent representable positions,
                # and exponent transitions where adding the origin loses one bit.
                values=np.array([1.,np.nextafter(np.float32(1.),np.float32(2.)),
                                 1.5,np.nextafter(np.float32(3.),np.float32(4.)),
                                 source/2.+1.,np.nextafter(np.float32(source),np.float32(source+1)),
                                 np.nextafter(np.float32(source+1),np.float32(0)),source+1.],dtype=np.float32)
                data=np.column_stack([values,np.roll(values,1),np.roll(values,3),
                                      np.linspace(-.25,.25,8,dtype=np.float32),np.zeros(8),
                                      -np.ones(8,dtype=np.float32)]).astype(np.float32)
                text='\n'.join(' '.join(format(float(x),'.9g') for x in row) for row in data)+'\n'
                run([scratch/'grid-test.exe',source,target],scratch,text)
                actual=np.fromfile(scratch/'converted.bin',dtype='>f4').reshape(8,6).astype(np.float64)
                (scratch/'converted.bin').unlink()
                before=np.mod(data[:,:3].astype(np.float64)-1.,source)*256./source
                after=(actual[:,:3]-1.)*256./target
                error=np.abs(after-before);error=np.minimum(error,256.-error)
                bound=float(np.spacing(np.float32(target+1)))*256./target
                assert np.max(error)<=bound
                assert np.all(actual[:,:3]>=1.) and np.all(actual[:,:3]<target+1.)
                before_v=data[:,3:].astype(np.float64)*51200./source
                after_v=actual[:,3:]*51200./target
                assert np.array_equal(before_v,after_v)
                if source==target:
                    # Compare identity away from the upper periodic image.
                    mask=data[:,:3] != source+1.
                    assert np.array_equal(actual[:,:3][mask],data[:,:3][mask])
                report['grid_controls'].append(dict(source=source,target=target,
                    max_position_error_mpc_h=float(error.max()),bound_mpc_h=bound,
                    physical_velocity_bitwise_equal=True,canonical_boundary=True))
        invalid=np.tile(np.array([1.,1.,1.,0.,0.,0.],dtype=np.float32),(8,1))
        invalid[0,0]=0.
        text='\n'.join(' '.join(map(str,row)) for row in invalid)+'\n'
        run([scratch/'grid-test.exe',16,16],scratch,text,success=False)
        run([scratch/'grid-test.exe',16,3],scratch,text,success=False)

        # Link original standalone entries with the same frozen runtime objects.
        for variant in ['legacy','v3']:
            folder=scratch/('original-'+variant);folder.mkdir()
            adapter=args.build_dir/variant
            for path in common_dir.glob('*.mod'):
                if path.name not in ['linkerlist.mod','structures.mod','bdmduplicaterules.mod']:
                    shutil.copy2(path,folder/path.name)
            if variant=='legacy':
                finder_obj=adapter/'finder.o'
                for name in ['linkerlist.mod','structures.mod','bdmduplicaterules.mod']:
                    shutil.copy2(adapter/name,folder/name)
                entry=subprocess.check_output(['git','show','a8c7715:PMP2bdm.f90'],cwd=args.repo,text=True)
            else:
                (folder/'finder.f90').write_text((args.repo/'PMP2linker.f90').read_text())
                run([build['compiler'],*builder.FLAGS,'-c','finder.f90'],folder)
                finder_obj=folder/'finder.o'
                entry=(args.repo/'PMP2bdm.f90').read_text()
            (folder/'original_entry.f90').write_text(entry)
            run([build['compiler'],*builder.FLAGS,'-o','original.exe',
                 *[common_dir/(name+'.o') for name in builder.OBJECTS],finder_obj,'original_entry.f90'],folder)
        tape=scratch/'shared-fi.bin'
        catalogues={}; members={}
        for source in [16,8,32]:
            for variant in ['v3','legacy']:
                folder=scratch/f'{variant}-{source}'
                fixture(folder,source)
                original_hashes={p.name:builder.sha(p) for p in folder.glob('PMcr*.DAT')}
                mode='write' if not tape.exists() else 'read'
                result=run([args.build_dir/variant/'replay.exe',1,16,1,mode,tape],folder)
                assert 'REPLAY COMPLETE' in result.stdout
                paths=list((folder/'CATALOGS').glob('Catshort*.DAT'));assert len(paths)==1
                rows=np.loadtxt(paths[0],skiprows=8,ndmin=2)
                assert rows.shape[0]>=1 and rows.shape[1]==24
                catalogues[(variant,source)]=paths[0].read_bytes()
                if source!=16:
                    assert catalogues[(variant,source)]==catalogues[(variant,16)]
                if variant=='v3':
                    raw=(folder/'repair-members.bin').read_bytes()
                    selected,candidates,mass=struct.unpack_from('>qqf',raw)
                    assert selected==1 and candidates>=selected and mass>0
                    candidate,count=struct.unpack_from('>qq',raw,20)
                    assert count==64
                    ids=struct.unpack_from('>64q',raw,120)
                    assert ids==tuple(range(1,65)) and len(raw)==120+8*64
                    members[source]=raw
                    if source!=16:assert members[source]==members[16]
                assert original_hashes=={p.name:builder.sha(p) for p in folder.glob('PMcr*.DAT')}
                report['fixture_controls'].append(dict(variant=variant,source_grid=source,analysis_grid=16,
                    rows=len(rows),catalogue_sha256=builder.sha(paths[0]),snapshot_unchanged=True,
                    same_physical_input_same_variant_byte_identical=True))
        for variant in ['legacy','v3']:
            folder=scratch/('standalone-'+variant);fixture(folder,16)
            run([scratch/('original-'+variant)/'original.exe'],folder,'1\n')
            path=next((folder/'CATALOGS').glob('Catshort*.DAT'))
            assert path.read_bytes()==catalogues[(variant,16)]
        report['original_standalone_matches']=['legacy','v3']
        report['fixed_density_sha256']=builder.sha(tape)
        report['fixed_density_size']=tape.stat().st_size
        assert tape.stat().st_size==28+4*16**3
        empty=scratch/'empty-v3';fixture(empty,16)
        empty_tape=scratch/'empty-fi.bin'
        empty_tape.write_bytes(struct.pack('>qqiff',16,64,1,1.,32.)+bytes(4*16**3))
        result=run([args.build_dir/'v3/replay.exe',1,16,1,'read',empty_tape],empty)
        assert 'REPLAY COMPLETE' in result.stdout
        raw=(empty/'repair-members.bin').read_bytes()
        selected,candidates,mass=struct.unpack('>qqf',raw)
        assert selected==candidates==0 and mass>0
        report['empty_catalogue_hook_passed']=True
        bad=scratch/'invalid-density';fixture(bad,16)
        bad_tape=scratch/'invalid-fi.bin';bad_tape.write_bytes(empty_tape.read_bytes()[:-4])
        result=run([args.build_dir/'v3/replay.exe',1,16,1,'read',bad_tape],bad,success=False)
        assert 'Fixed density byte count mismatch' in result.stdout
        assert not list((bad/'CATALOGS').iterdir())
        report['incomplete_density_rejected_before_finder']=True
        report['completed']=True
    finally:
        args.receipt.parent.mkdir(parents=True,exist_ok=True)
        args.receipt.write_text(json.dumps(report,indent=2)+'\n')
    if report['completed']:
        shutil.rmtree(scratch)
        print('All replay controls passed; scratch removed:',args.receipt)


if __name__=='__main__':
    main()
