"""Small independent Fourier, grid-origin, RNG and staggered-time checks."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import struct
import tempfile
import time
import numpy as np
from build import build
from generate import generate
from configure import configure

HERE = Path(__file__).resolve().parent


def receipt(path: Path) -> dict:
    values = {}
    for line in path.read_text().splitlines():
        key, value = line.split('=', 1)
        values[key.strip()] = value.strip()
    return values


def fields(path: Path, nrow: int, components: int) -> np.ndarray:
    data = np.fromfile(path, dtype='>f4').astype('f4')
    assert data.size == components*nrow**3, (path, data.size)
    return np.stack([x.reshape((nrow,)*3, order='F') for x in data.reshape(components, -1)])


def basis(n: int) -> np.ndarray:
    """Independent analytic Fourier series; FFT5 backward has no N factor."""
    phase = 2*np.pi*np.arange(n)/n
    result = np.empty((n, n))
    result[:, 0] = 1
    result[:, -1] = (-1.)**np.arange(n)
    for frequency in range(1, n//2):
        result[:, 2*frequency-1] = np.cos(frequency*phase)
        result[:, 2*frequency] = np.sin(frequency*phase)
    return result


def analytic_synthesis(packed: np.ndarray) -> np.ndarray:
    b = basis(packed.shape[0])
    result = np.einsum('ai,ijk->ajk', b, packed.astype('f8'))
    result = np.einsum('bj,ajk->abk', b, result)
    return np.einsum('ck,abk->abc', b, result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--compiler', choices=['gfortran', 'ifx'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(compiler=args.compiler, cases=[], comparisons=[], checks={},
                  test_source_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                      for name in ('checks.py','probe.f90','configure.py','build.py','generate.py','controls.f90','entry.f90')})
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_STACKSIZE='128M')
    with tempfile.TemporaryDirectory(prefix='bdm-matched-ic-check-') as temp:
        work = Path(temp)
        binpath = work / 'build'
        report['build'] = build(args.repo.resolve(), binpath, args.compiler, checked=True)
        flags = report['build']['flags']
        source = generate(args.repo.resolve()).split('! Campaign-only standalone entry:')[0]
        (binpath/'matched_module.f90').write_text(source)
        shutil.copy2(HERE/'probe.f90', binpath/'probe.f90')
        cmd = [args.compiler, *flags, '-o', 'probe.exe', 'matched_module.f90', 'probe.f90',
               'PMP2mod_tools.o', 'PMP2mod_random.o', 'PMP2mod_fft5.o']
        subprocess.run(cmd, cwd=binpath, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (work/'TableSeeds.dat').write_text('Seeds for matched IC controls\n1234567 1\n7654321 2\n')
        runs = {}

        def run(label, nrow, ngrid, alpha, threads=1, step=0.0001, master=32, seed=1, norm=False):
            path = work/label
            path.mkdir()
            control = ('&matched_ic\n'
                       f' ic_master_nrow={master}, ic_origin_ngrid=64, ic_alpha={alpha:.17g},\n'
                       f' ic_normalize_only={".true." if norm else ".false."}\n/\n')
            (path/'matched_ic.nml').write_text(control)
            completed = subprocess.run([str(binpath/'probe.exe'), str(nrow), str(ngrid), str(step)],
                                       input=f'{seed}\n', text=True, cwd=path,
                                       env=dict(env, OMP_NUM_THREADS=str(threads)), capture_output=True)
            if completed.returncode:
                raise RuntimeError(f'{label} rc={completed.returncode}: {completed.stdout}\n{completed.stderr}')
            r = receipt(path/'matched_ic_receipt.txt')
            record = dict(label=label, nrow=nrow, ngrid=ngrid, threads=threads, receipt=r,
                          control_sha256=hashlib.sha256(control.encode()).hexdigest())
            report['cases'].append(record)
            if norm:
                return r
            packed = fields(path/'packed.bin', nrow, 3)
            real = fields(path/'real.bin', nrow, 3)
            part = fields(path/'particles.bin', nrow, 6)
            assert np.isfinite(packed).all() and np.isfinite(real).all() and np.isfinite(part).all()
            assert np.count_nonzero(packed[:, -1, :, :]) == 0
            assert np.count_nonzero(packed[:, :, -1, :]) == 0
            assert np.count_nonzero(packed[:, :, :, -1]) == 0
            assert np.count_nonzero(packed[:, 0, 0, 0]) == 0
            errors = []
            for axis in range(3):
                expected = analytic_synthesis(packed[axis])
                error = np.max(np.abs(real[axis]-expected))
                assert error <= 5.e-7*max(1., np.max(np.abs(expected))), (label, error)
                errors.append(float(error))
            record['analytic_fft_max_abs_error'] = errors
            origin = float(r['origin_mpc_h'])
            xcons = float(r['xcons'])
            vcons = float(r['vcons'])
            coords = np.indices((nrow,)*3).astype('f8')*256/nrow
            physical = (part[:3].astype('f8')-1)*256/ngrid
            expected = np.mod(coords - xcons*real.astype('f8')*256/ngrid + origin, 256)
            period_error = np.abs((physical-expected+128)%256-128)
            max_position_error = float(period_error.max())
            assert max_position_error <= 5.e-5, (label, max_position_error)
            # Native output momentum is evaluated at a_v = a_x - da/2.
            a = np.float32(1./101.)
            av = np.float32(a - np.float32(step)/np.float32(2.))
            fact = np.sqrt(np.float32(.307) + np.float32(.693)*av**np.int32(3))
            expected_ratio = -float(av/a)*np.sqrt(float(av))*float(fact)
            assert np.isclose(vcons/xcons, expected_ratio, rtol=3.e-7)
            assert np.allclose(part[3:], vcons*real, rtol=2.e-7, atol=0.)
            record['physical_position_max_abs_error_mpc_h'] = max_position_error
            record['files_sha256'] = {name: hashlib.sha256((path/name).read_bytes()).hexdigest()
                                     for name in ('packed.bin','real.bin','particles.bin','matched_modes.bin')}
            result = dict(packed=packed, real=real, part=part, receipt=r, nrow=nrow, ngrid=ngrid)
            runs[label] = result
            return r

        master = run('master', 32, 64, -1)
        alpha = float(master['alpha'])
        normalize = run('normalization_only', 32, 64, -1, norm=True, threads=4)
        assert normalize['alpha'] == master['alpha'] and normalize['spectrum_sum'] == master['spectrum_sum']
        for nrow in (8,16):
            run(f'n{nrow}', nrow, 64, alpha)
        for threads in (2,4):
            run(f'master_t{threads}', 32, 64, -1, threads=threads)
        for ngrid in (32,128):
            run(f'n16_g{ngrid}',16,ngrid,alpha)
        run('half_step',32,64,alpha,step=0.00005)
        run('second_seed',32,64,-1,seed=2)
        for nrow in (8,16,32):
            run(f'production_stride_n{nrow}',nrow,64,alpha,master=1024)

        def compare_modes(left, right):
            a, b = runs[left], runs[right]
            n = min(a['nrow'], b['nrow'])
            packed_k = (np.arange(n-1)+1)//2
            px,py,pz = np.meshgrid(packed_k,packed_k,packed_k,indexing='ij')
            common = px*px+py*py+pz*pz < (n//2)**2
            assert np.array_equal(a['packed'][:, :n-1,:n-1,:n-1][:,common],
                                  b['packed'][:, :n-1,:n-1,:n-1][:,common]), (left,right)
            fa = np.fft.fftn(a['real'].astype('f8'), axes=(1,2,3))/a['nrow']**3
            fb = np.fft.fftn(b['real'].astype('f8'), axes=(1,2,3))/b['nrow']**3
            modes = np.arange(-n//2+1,n//2)
            ix,iy,iz = np.meshgrid(modes,modes,modes,indexing='ij')
            mask = (ix*ix+iy*iy+iz*iz)<(n//2)**2
            ca = fa[:,ix%a['nrow'],iy%a['nrow'],iz%a['nrow']][:,mask]
            cb = fb[:,ix%b['nrow'],iy%b['nrow'],iz%b['nrow']][:,mask]
            error = float(np.max(np.abs(ca-cb)))
            assert error < 8.e-8*max(1.,np.max(np.abs(cb))), (left,right,error)
            report['comparisons'].append(dict(left=left,right=right,packed_common_modes='bit-identical',
                                               physical_fourier_max_abs_error=error, complex_coefficients=ca.size))

        for small in ('n8','n16','n16_g32','n16_g128'):
            compare_modes(small,'master')
        for small in ('production_stride_n8','production_stride_n16'):
            compare_modes(small,'production_stride_n32')
        for label in ('master_t2','master_t4'):
            for what in ('packed','real','part'):
                assert np.array_equal(runs['master'][what],runs[label][what])
            assert runs[label]['receipt']['alpha'] == master['alpha']
        assert np.array_equal(runs['master']['packed'],runs['half_step']['packed'])
        assert np.array_equal(runs['master']['part'][:3],runs['half_step']['part'][:3])
        assert not np.array_equal(runs['master']['part'][3:],runs['half_step']['part'][3:])
        assert not np.array_equal(runs['master']['packed'],runs['second_seed']['packed'])
        for pair in (('n16','n16_g32'),('n16','n16_g128')):
            assert np.array_equal(runs[pair[0]]['real'],runs[pair[1]]['real'])
        # Exercise the real entry point, Setup/Pk parsing and on-disk PM format.
        # The fixture is deliberately small but uses native NPAGE=1024**2.
        entry_parent = work/'entry'
        entry_parent.mkdir()
        entry = entry_parent/'Run1'
        entry.mkdir()
        a0 = float(np.float32(1./101.))
        setup = ['Matched IC standalone integration fixture', 'matched-native-ic',
                 a0, .0001, .03, 256., .8228, .6777, .307, .693, .048252,
                 32, 64, 0, 4., 1., 100., 0., 0., 0., 1, 0., 0,
                 20, 0, 0, 0, 0, 0, 0, 0, 0,
                 'MG', 0, 0, 1, 0, 0, 0, 4, 1.e-8,
                 'fR', 1, 1.e-6, 'DGP', 1, 0, 1.,
                 'sym', .5, .01, 1., 'kmf', 0, 0, 2, 1., 1.,
                 'csf', 1, 1, 1., 1.]
        (entry_parent/'Setup.dat').write_text('\n'.join(map(str,setup))+'\n')
        (entry_parent/'TableSeeds.dat').write_text('Seeds\n1234567 1\n')
        (entry_parent/'PkTable.dat').write_text(
            'OmB = 0.048252\nOmC = 0.258748\nOmL = 0.693\nOmM = 0.307\ns8 = 0.8228\n'
            '1e-6 1\n1e6 1\n')
        configure(entry,None,master_nrow=32,origin_ngrid=64)
        completed = subprocess.run([str(binpath/'PMP2start.matched.exe')],input='1\n', text=True,
                                   cwd=entry,env=dict(env,OMP_NUM_THREADS='2'),capture_output=True)
        if completed.returncode:
            raise RuntimeError(f'Native entry failed: {completed.stdout}\n{completed.stderr}')
        entry_receipt = receipt(entry/'matched_ic_receipt.txt')
        assert entry_receipt['alpha'] == master['alpha']
        raw_header = (entry/'PMcrd.DAT').read_bytes()
        length = struct.unpack('>i',raw_header[:4])[0]
        assert length == len(raw_header)-8 and raw_header[:4] == raw_header[-4:]
        assert length == 529, (length, len(raw_header))
        initial_a, header_aexp0 = struct.unpack('>ff',raw_header[4+45:4+45+8])
        assert initial_a == a0 and header_aexp0 == a0
        header_nrow, header_ngrid = struct.unpack('>ii',raw_header[4+45+12*4:4+45+14*4])
        assert (header_nrow,header_ngrid) == (32,64)
        data = np.frombuffer((entry/'PMcrs0.DAT').read_bytes(), dtype='>f4')
        assert data.size == 6*1024**2
        unpadded = data.reshape(6,1024**2)[:,:32**3]
        # Fortran particle indexing has x fastest; each cube was restored with order F.
        expected = np.stack([f.ravel(order='F') for f in runs['master']['part']])
        assert np.array_equal(unpadded,expected)
        assert (entry/'pt.dat').read_bytes() == struct.pack('>ifi',4,np.float32(.0001),4)
        consumer_parent = work/'consumer'
        consumer_parent.mkdir()
        consumer = consumer_parent/'Run1'
        changed_setup = setup.copy()
        changed_setup[3] = .00005
        changed_setup[4] = .02
        changed_setup[11:13] = [16,128]
        (consumer_parent/'Setup.dat').write_text('\n'.join(map(str,changed_setup))+'\n')
        for name in ('PkTable.dat','TableSeeds.dat'):
            shutil.copy2(entry_parent/name,consumer_parent/name)
        control = configure(consumer,entry,master_nrow=32,origin_ngrid=64)
        assert control['alpha_requested'] == float(entry_receipt['alpha'])
        assert control['inputs']['astep'] != float(entry_receipt['astep'])
        assert configure(consumer,entry,master_nrow=32,origin_ngrid=64) == control
        initial_controls = (consumer/'matched_ic.nml').read_bytes()
        (consumer/'matched_ic.nml').write_text('changed controls\n')
        try:
            configure(consumer,entry,master_nrow=32,origin_ngrid=64)
        except FileExistsError:
            pass
        else:
            raise AssertionError('Existing different controls were overwritten')
        assert (consumer/'matched_ic.nml').read_text() == 'changed controls\n'
        (consumer/'matched_ic.nml').write_bytes(initial_controls)
        # Positive rejection controls: changed spectrum and altered master controls.
        with (consumer_parent/'PkTable.dat').open('a') as out:
            out.write('2e6 1\n')
        try:
            configure(consumer,entry,master_nrow=32,origin_ngrid=64)
        except ValueError:
            pass
        else:
            raise AssertionError('Changed consumer P(k) accepted')
        shutil.copy2(entry_parent/'PkTable.dat',consumer_parent/'PkTable.dat')
        with (entry/'matched_ic.nml').open('a') as out:
            out.write('! changed after normalizing\n')
        try:
            configure(consumer,entry,master_nrow=32,origin_ngrid=64)
        except ValueError:
            pass
        else:
            raise AssertionError('Changed master controls accepted')
        report['entry_smoke'] = dict(receipt=entry_receipt,header_bytes=len(raw_header),
                                    particle_bytes=(entry/'PMcrs0.DAT').stat().st_size,
                                    particle_prefix_matches_probe=True, explicit_byte_record_units=True)
        report['checks'] = dict(no_overwrite_control_publication_and_identical_resume=True,
            master_input_hash_binding_and_rejection_controls=True,
            native_entry_setup_and_power_table=True, native_big_endian_pm_output=True,
            initialized_header_epoch=True, analytic_fft_synthesis=True, common_physical_fourier_modes=True,
            zero_uniform_and_nyquist_modes=True, no_nrow_fft_amplitude_factor=True, fixed_physical_origin=True,
            three_mesh_sizes=True, same_field_across_mesh_sizes=True, master_alpha_thread_independent=True,
            byte_identical_1_2_4_thread_fields_and_particles=True, normalization_only_matches_full=True,
            half_step_identical_positions=True, half_step_correct_native_velocity_time=True,
            second_seed_changes_field=True, production_stride_1024_common_modes=True)
        report['temporary_entries_removed_on_exit'] = sum(1 for _ in work.rglob('*'))
    report['all_passed'] = True
    report['elapsed_seconds'] = time.monotonic()-start
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(all_passed=True, checks=report['checks'], cases=len(report['cases']),
                          elapsed_seconds=report['elapsed_seconds']), indent=2))


if __name__ == '__main__':
    main()
