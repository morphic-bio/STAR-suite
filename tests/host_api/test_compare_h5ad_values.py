from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from compare_h5ad_values import compare


class H5ValuesTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.a, self.b = [Path(self.tmp.name) / name for name in ('a.h5ad', 'b.h5ad')]
        for path in (self.a, self.b):
            with h5py.File(path, 'w') as handle:
                handle.attrs['encoding-type'] = 'anndata'
                handle.create_dataset('X', data=np.array([1., float('nan'), 3.]))
                handle.create_dataset('obs/index', data=np.array(['A', 'B'], dtype=h5py.string_dtype()))

    def test_identical_including_nan(self):
        self.assertTrue(compare(self.a, self.b)['pass'])

    def test_numeric_difference(self):
        with h5py.File(self.b, 'r+') as handle:
            handle['X'][0] = 2
        self.assertFalse(compare(self.a, self.b)['pass'])

    def test_attribute_difference(self):
        with h5py.File(self.b, 'r+') as handle:
            handle.attrs['encoding-type'] = 'wrong'
        self.assertFalse(compare(self.a, self.b)['pass'])

    def test_missing_column(self):
        with h5py.File(self.b, 'r+') as handle:
            del handle['obs/index']
        self.assertFalse(compare(self.a, self.b)['pass'])

    def test_axis_order(self):
        with h5py.File(self.b, 'r+') as handle:
            handle['obs/index'][:] = ['B', 'A']
        self.assertFalse(compare(self.a, self.b)['pass'])

    def test_provenance_roots_explicit(self):
        for path, root in ((self.a, '/one'), (self.b, '/two')):
            with h5py.File(path, 'r+') as handle:
                handle.create_dataset('uns/source', data=root + '/matrix.mtx')
        self.assertFalse(compare(self.a, self.b)['pass'])
        self.assertTrue(compare(self.a, self.b, ['/one'], ['/two'])['pass'])


if __name__ == '__main__':
    unittest.main()
