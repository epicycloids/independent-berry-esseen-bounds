"""Reject incorrect claims attached to otherwise valid replay receipts."""
import copy
import json
import unittest

from check_saved import ROOT, CertificateError, load_inputs
from check_replays import validate_receipt


class ReplayReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cover, cls.result = load_inputs()
        cls.receipt = json.loads((ROOT/'certificates/selected-replay.json').read_text())

    def validate(self, value):
        return validate_receipt(value, 'selected', self.cover, self.result)

    def test_valid_selected_receipt(self):
        self.assertEqual(self.validate(self.receipt)['accepted'], 4944)

    def test_rejects_scope_promotion(self):
        value = dict(self.receipt, all_bands_numerically_replayed=True)
        with self.assertRaisesRegex(CertificateError, 'full numerical'):
            self.validate(value)

    def test_rejects_wrong_wrapper(self):
        value = dict(self.receipt, wrapper_sha256={})
        with self.assertRaisesRegex(CertificateError, 'wrapper digest'):
            self.validate(value)

    def test_rejects_changed_selected_value(self):
        value = copy.deepcopy(self.receipt)
        value['bands'][0]['upper'] = 0.4
        with self.assertRaisesRegex(CertificateError, 'record mismatch'):
            self.validate(value)


if __name__ == '__main__':
    unittest.main()
