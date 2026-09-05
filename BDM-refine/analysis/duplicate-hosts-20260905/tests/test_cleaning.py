"""Meaningful regression checks for row frames, periodic groups and byte identity."""
from pathlib import Path
import contextlib
import io
import json
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from catalogue_core import COLUMNS, FIELDS, audit, components, edge_flags
from clean_catalogue import clean, load_receipt


def fixture():
    a = np.ones((9, 24), dtype=np.float64)
    a[:, 0] = [0.01, 31.99, 10.30, 10.0, 10.15, 20.0, 20.05, 25., 25.05]
    a[:, 1:3] = 8.
    a[:, 3:6] = 0.
    a[6, 3] = 200.  # Close equal-mass high-relative-speed merger control.
    a[:, 6:8] = 1.e13
    a[8, 7] = 1.1e13  # Fails the aperture-mass consistency cut.
    a[:, 8] = 500.
    a[:, 11] = np.arange(1001,1010)
    a[:, 13] = 1000.
    return a


class CleaningTests(unittest.TestCase):
    def test_periodic_chain_and_merger_controls(self):
        a = fixture()
        data = {k:a[:,COLUMNS.index(k)] for k in FIELDS}
        result, parent, drop, *_ = audit(data, box=32.)
        self.assertEqual(np.flatnonzero(drop).tolist(), [1,3,4])
        self.assertEqual(parent[[2,3,4]].tolist(), [2,2,2])
        self.assertEqual(result['remaining_strict_pairs'], 0)
        self.assertEqual(result['remaining_same_bound_pairs_below_0p2'], 2)
        p, d = components(5, np.array([[4,3],[4,2],[3,2]]))
        q, e = components(5, np.array([[2,3],[2,4],[3,4]]))
        np.testing.assert_array_equal(p,q)
        np.testing.assert_array_equal(d,e)

    def test_strict_boundaries_and_unequal_count(self):
        a = fixture()[:2].copy()
        a[:,3:6]=0.
        data = {k:a[:,COLUMNS.index(k)] for k in FIELDS}
        edges = np.array([[0,1]])
        self.assertFalse(edge_flags(data, edges, np.array([.2]))['strict'][0])
        data['vx'][1]=5.
        self.assertFalse(edge_flags(data, edges, np.array([.1]))['strict'][0])
        data['vx'][1]=0.
        data['Nparticles'][1]+=1.
        self.assertFalse(edge_flags(data, edges, np.array([.1]))['strict'][0])

    def test_ascii_and_hdf5_preserve_original_frame(self):
        base = Path(__file__).resolve().parent/'cases'
        base.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=base) as directory:
            root = Path(directory)
            a = fixture()
            header = ['synthetic', ' A = 0.8 Step = 1', 'grid', 'cosmology',
                      'radial bins', 'particle mass', 'overdensity', ' '.join(COLUMNS)]
            text = root/'source.DAT'
            with open(text, 'w') as f:
                f.write('\n'.join(header)+'\n')
                np.savetxt(f,a,fmt='%.17g')
            hdf = root/'source.hdf5'
            with h5py.File(hdf, 'w') as f:
                f.attrs['box_size']=32.
                for j,k in enumerate(COLUMNS):
                    f.create_dataset(k,data=a[:,j],compression='gzip')
                f.create_dataset('scalar_metadata',data=99)
            results=[]
            for source in [text,hdf]:
                output = root/('clean'+source.suffix)
                sidecar = root/('mask'+source.suffix+'.hdf5')
                with contextlib.redirect_stdout(io.StringIO()):
                    r = clean(source,output,sidecar,
                              dict(gravity='test',imodel=0,ibox=1,redshift=.25),box=32.)
                self.assertTrue(r['bitwise_validation']['all_surviving_columns_byte_identical'])
                self.assertEqual(r['removed_rows'],3)
                self.assertEqual(load_receipt(sidecar)['source_rows'],9)
                with h5py.File(sidecar,'r') as f:
                    mask=f['drop_mask'][...]
                    groups=json.loads(f['groups_json'][...].tobytes())
                    self.assertEqual(groups[0]['member_Nhalo'],[1001,1002])
                    self.assertEqual(groups[1]['kept_Nhalo'],1003)
                results.append(mask)
                if source == text:
                    # Fixed skiprows=8 still works: the extra provenance is a comment.
                    cleaned=np.loadtxt(output,skiprows=8)
                    self.assertEqual(cleaned.tobytes(),a[~mask].tobytes())
                    original_lines=text.read_bytes().splitlines(keepends=True)[8:]
                    kept_lines=output.read_bytes().splitlines(keepends=True)[9:]
                    self.assertEqual(kept_lines,[line for i,line in enumerate(original_lines) if not mask[i]])
                else:
                    with h5py.File(output,'r') as f:
                        self.assertEqual(f['Nhalo'][...].tobytes(),a[~mask,11].tobytes())
                        self.assertEqual(f['scalar_metadata'][()],99)
                with self.assertRaises(FileExistsError):
                    clean(source,output,sidecar,{},box=32.)
            np.testing.assert_array_equal(*results)


if __name__ == '__main__':
    unittest.main()
