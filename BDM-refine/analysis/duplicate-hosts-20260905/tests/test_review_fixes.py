"""Regressions for storage precision and interrupted two-file publication."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import clean_catalogue
import run_campaign
from catalogue_core import COLUMNS, FIELDS, audit, sha256
from clean_catalogue import clean, load_receipt
from test_cleaning import fixture


IDENTITY = dict(gravity='test', imodel=1, ibox=1, redshift=.25, snapnum=137)


def write_source(root, suffix, values):
    source = root/('source'+suffix)
    if suffix == '.DAT':
        header = ['synthetic', ' A = 0.8 Step = 137', 'grid', 'cosmology',
                  'radial bins', 'particle mass', 'overdensity', ' '.join(COLUMNS)]
        with source.open('w') as f:
            f.write('\n'.join(header)+'\n')
            # Preserve the exact float32 values when comparing formats.
            np.savetxt(f, values, fmt='%.17g')
    else:
        with h5py.File(source, 'w') as f:
            f.attrs['labels'] = np.array(['test', 'original'], dtype=object)
            f.create_dataset('scalar_metadata', data=99)
            for j, key in enumerate(COLUMNS):
                column = values[:, j]
                if key in ['Nhalo', 'Nparticles']:
                    column = column.astype(np.int64)
                d = f.create_dataset(key, data=column, compression='gzip')
                d.attrs['unit'] = 'original unit'
    return source


def quiet_clean(*args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return clean(*args, **kwargs)


class PrecisionTests(unittest.TestCase):
    def test_float32_diagnostics_match_the_same_values_in_float64(self):
        a = fixture()[:2].astype(np.float32)
        a[:, 0] = [1., 1.05]
        a[:, 1:3] = 1.
        cases = {}
        cases['mass difference'] = a.copy()
        cases['mass difference'][:, 6] = 1e11
        cases['mass difference'][:, 7] = [1e12, 1.0115794e12]
        cases['separation'] = a.copy()
        cases['separation'][:, :3] = [[0, 0, 0], [.12, .16, 0]]
        cases['relative speed'] = a.copy()
        cases['relative speed'][1, 3:6] = [3., np.nextafter(np.float32(4), np.float32(0)), 0.]
        cases['mass selection and bins'] = a.copy()
        cases['mass selection and bins'][:, 6] = 1e11
        cases['mass selection and bins'][:, 7] = np.float32(10**12.4)
        cases['radial velocity'] = a.copy()
        cases['radial velocity'][:, 0] = [0., 1.2]
        cases['radial velocity'][:, 3] = [123.4567, -98.7654]
        for name, values in cases.items():
            with self.subTest(quantity=name):
                data = {k: values[:, COLUMNS.index(k)] for k in FIELDS}
                reference = {k: v.astype(np.float64) for k, v in data.items()}
                expected, _, expected_drop, *_ = audit(reference, box=32., velocities=True)
                actual, _, drop, *_ = audit(data, box=32., velocities=True)
                if name in ['mass difference', 'separation', 'relative speed']:
                    self.assertEqual(expected_drop.tolist(), [False, True])
                np.testing.assert_array_equal(drop, expected_drop)
                self.assertEqual(actual, expected)
                self.assertTrue(all(v.dtype == np.float32 for v in data.values()))

    def test_float32_hdf5_and_ascii_masks_agree_and_storage_dtypes_survive(self):
        a = fixture()[:2].astype(np.float32)
        a[:, 0] = [1., 1.05]
        a[:, 6] = 1e11
        a[:, 7] = [1e12, 1.0115794e12]
        masses = a[:, 7]
        self.assertLess(np.diff(np.log10(masses.astype(np.float64)))[0], .005)
        self.assertGreater(np.diff(np.log10(masses))[0], .005)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            masks = []
            for suffix in ['.DAT', '.hdf5']:
                source = write_source(root, suffix, a)
                output = root/('clean'+suffix)
                sidecar = root/('mask'+suffix+'.hdf5')
                receipt = quiet_clean(source, output, sidecar, IDENTITY)
                self.assertEqual(receipt['removed_rows'], 1)
                with h5py.File(sidecar, 'r') as f:
                    mask = f['drop_mask'][...]
                    masks.append(mask)
                if suffix == '.hdf5':
                    with h5py.File(source, 'r') as raw, h5py.File(output, 'r') as cleaned:
                        for key in COLUMNS:
                            self.assertEqual(raw[key].dtype, cleaned[key].dtype)
                            self.assertEqual(raw[key][...][~mask].tobytes(), cleaned[key][...].tobytes())
            np.testing.assert_array_equal(*masks)


class PublicationTests(unittest.TestCase):
    def test_sidecar_write_failure_publishes_nothing_and_cleans_staging(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = write_source(root, '.DAT', fixture())
            output, sidecar = root/'clean.DAT', root/'clean.cleaning.hdf5'
            with mock.patch.object(clean_catalogue, 'json_dataset', side_effect=OSError('injected sidecar write')):
                with self.assertRaisesRegex(OSError, 'injected'):
                    quiet_clean(source, output, sidecar, IDENTITY)
            self.assertFalse(output.exists())
            self.assertFalse(sidecar.exists())
            self.assertEqual(list(root.glob('*.partial-*')), [])
            quiet_clean(source, output, sidecar, IDENTITY)
            self.assertTrue(sidecar.exists())

    def test_campaign_recovers_an_orphan_without_replacing_it(self):
        for suffix in ['.DAT', '.hdf5']:
            with self.subTest(format=suffix), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                a = fixture()
                a[1, 0] = 1023.99
                source = write_source(root, suffix, a)
                output, sidecar = root/('clean'+suffix), root/'clean.cleaning.hdf5'
                original_link = clean_catalogue.os.link

                def interrupt_receipt_publish(src, dst):
                    if Path(dst) == sidecar:
                        raise OSError('injected interruption before receipt publication')
                    original_link(src, dst)

                with mock.patch.object(clean_catalogue.os, 'link', side_effect=interrupt_receipt_publish):
                    with self.assertRaisesRegex(OSError, 'injected'):
                        quiet_clean(source, output, sidecar, IDENTITY)
                self.assertTrue(output.exists())
                self.assertFalse(sidecar.exists())
                self.assertEqual(list(root.glob('*.partial-*')), [])
                before_hash, before_inode, source_hash = sha256(output), output.stat().st_ino, sha256(source)
                manifest = root/'manifest.json'
                manifest.write_text(json.dumps([dict(source=str(source), **IDENTITY)]))
                argv = ['run_campaign', str(manifest), str(root), '--workers', '1', '--worker', '0']
                with mock.patch.object(sys, 'argv', argv), \
                     mock.patch.object(run_campaign, 'output_paths', return_value=(output, sidecar)), \
                     contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    run_campaign.main()
                    # The next resume must use the normal verified-receipt path.
                    run_campaign.main()
                receipt = load_receipt(sidecar)
                self.assertTrue(receipt['publication']['recovered_without_receipt'])
                self.assertEqual(receipt['removed_rows'], 3)
                self.assertEqual(receipt['output_sha256'], before_hash)
                self.assertEqual(sha256(output), before_hash)
                self.assertEqual(output.stat().st_ino, before_inode)
                self.assertEqual(sha256(source), source_hash)
                self.assertEqual(list(root.glob('*.partial-*')), [])

    def test_recovery_accepts_equivalent_hdf5_container_layout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = write_source(root, '.hdf5', fixture())
            output, sidecar = root/'clean.hdf5', root/'clean.cleaning.hdf5'
            quiet_clean(source, output, sidecar, IDENTITY)
            old_hash = sha256(output)
            sidecar.unlink()
            repacked = root/'repacked.hdf5'
            with h5py.File(output, 'r') as src, h5py.File(repacked, 'w') as dst:
                dst.attrs.update(src.attrs)
                for name in reversed(list(src)):
                    src.copy(name, dst)
            repacked.replace(output)
            self.assertNotEqual(sha256(output), old_hash)
            before_hash, before_inode = sha256(output), output.stat().st_ino
            receipt = quiet_clean(source, output, sidecar, IDENTITY)
            self.assertTrue(receipt['publication']['recovered_without_receipt'])
            self.assertEqual(receipt['output_sha256'], before_hash)
            self.assertEqual(output.stat().st_ino, before_inode)

    def test_recovery_rejects_unverified_ascii_without_overwriting(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = write_source(root, '.DAT', fixture())
            output, sidecar = root/'clean.DAT', root/'clean.cleaning.hdf5'
            quiet_clean(source, output, sidecar, IDENTITY)
            sidecar.unlink()
            with output.open('ab') as f:
                f.write(b'Unexpected extra row\n')
            before_hash = sha256(output)
            with self.assertRaisesRegex(ValueError, 'Existing catalogue'):
                quiet_clean(source, output, sidecar, IDENTITY)
            self.assertEqual(sha256(output), before_hash)
            self.assertFalse(sidecar.exists())
            self.assertEqual(list(root.glob('*.partial-*')), [])

    def test_recovery_checks_all_hdf5_values_dtypes_and_attributes(self):
        for damage in ['retained value', 'scalar metadata', 'file attribute', 'dataset attribute', 'dtype']:
            with self.subTest(damage=damage), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = write_source(root, '.hdf5', fixture())
                output, sidecar = root/'clean.hdf5', root/'clean.cleaning.hdf5'
                quiet_clean(source, output, sidecar, IDENTITY)
                sidecar.unlink()
                with h5py.File(output, 'r+') as f:
                    if damage == 'retained value':
                        f['x'][0] += 1
                    elif damage == 'scalar metadata':
                        f['scalar_metadata'][()] = 100
                    elif damage == 'file attribute':
                        f.attrs['original_rows'] = 999
                    elif damage == 'dataset attribute':
                        f['x'].attrs['unit'] = 'different unit'
                    else:
                        values = f['x'][...].astype(np.float32)
                        del f['x']
                        f.create_dataset('x', data=values)
                        f['x'].attrs['unit'] = 'original unit'
                before_hash = sha256(output)
                with self.assertRaisesRegex(ValueError, 'Existing catalogue'):
                    quiet_clean(source, output, sidecar, IDENTITY)
                self.assertEqual(sha256(output), before_hash)
                self.assertFalse(sidecar.exists())
                self.assertEqual(list(root.glob('*.partial-*')), [])


if __name__ == '__main__':
    unittest.main()
