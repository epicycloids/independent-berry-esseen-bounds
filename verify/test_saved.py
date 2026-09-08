"""Focused rejection tests for the public certificate interface."""
import base64
import copy
import hashlib
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import zlib

from check_saved import ROOT, CertificateError, load_inputs, tree_counts, validate_cover


def encoded(text):
    raw = text.encode()
    return base64.b64encode(zlib.compress(raw)).decode(), hashlib.sha256(raw).hexdigest()


class SavedCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cover, cls.result = load_inputs()

    def test_valid_saved_cover(self):
        self.assertEqual(validate_cover(self.cover, self.result)['counts']['bands'], 767)

    def test_incomplete_tree(self):
        with self.assertRaisesRegex(CertificateError, 'incomplete'):
            tree_counts(*encoded('DA'))

    def test_trailing_tree(self):
        with self.assertRaisesRegex(CertificateError, 'trailing'):
            tree_counts(*encoded('AA'))

    def test_invalid_instruction(self):
        with self.assertRaisesRegex(CertificateError, 'Unknown'):
            tree_counts(*encoded('?'))

    def test_tree_digest_changed(self):
        tree, _ = encoded('A')
        with self.assertRaisesRegex(CertificateError, 'digest'):
            tree_counts(tree, '0'*64)

    def test_missing_band(self):
        cover = dict(self.cover, bands=self.cover['bands'][1:])
        with self.assertRaisesRegex(CertificateError, 'Gap'):
            validate_cover(cover, self.result)

    def test_incomplete_band(self):
        cover = dict(self.cover, bands=[dict(self.cover['bands'][0], complete=False), *self.cover['bands'][1:]])
        with self.assertRaisesRegex(CertificateError, 'Incomplete'):
            validate_cover(cover, self.result)

    def test_unsupported_upper_bound(self):
        cover = dict(self.cover, bands=[dict(self.cover['bands'][0], upper=.476), *self.cover['bands'][1:]])
        with self.assertRaisesRegex(CertificateError, 'stated constant'):
            validate_cover(cover, self.result)

    def test_nonfinite_upper_bound(self):
        cover = dict(self.cover, bands=[dict(self.cover['bands'][0], upper=float('nan')), *self.cover['bands'][1:]])
        with self.assertRaisesRegex(CertificateError, 'Invalid upper'):
            validate_cover(cover, self.result)

    def test_changed_kernel(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'certificates').mkdir()
            (root/'verify/kernel').mkdir(parents=True)
            shutil.copyfile(ROOT/'result.json', root/'result.json')
            shutil.copyfile(ROOT/'certificates/cover.json', root/'certificates/cover.json')
            for name in self.cover['kernel_sha256']:
                shutil.copyfile(ROOT/'verify/kernel'/name, root/'verify/kernel'/name)
            with (root/'verify/kernel/interval_bounds.py').open('ab') as stream:
                stream.write(b'\n# deliberate test mutation\n')
            with self.assertRaisesRegex(CertificateError, 'Arithmetic source changed'):
                load_inputs(root)

    def test_geometric_replay_rejects_changed_root(self):
        sys.path.insert(0, str(ROOT/'verify/kernel'))
        from stable_cover import replay_band
        record = copy.deepcopy(self.cover['bands'][0])
        record['root'][1] /= 2
        with self.assertRaises(AssertionError):
            replay_band(record, numerical=False)


if __name__ == '__main__':
    unittest.main()
