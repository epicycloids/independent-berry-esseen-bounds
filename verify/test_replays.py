"""Checks of the saved numerical reevaluation of three intervals."""
from copy import deepcopy
import json
import unittest

from check_saved import ROOT,CertificateError,load_inputs
from check_replays import validate_selected


class SelectedReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cover,cls.result=load_inputs()
        cls.replay=json.loads((ROOT/"certificates/selected-replay.json").read_text())

    def test_selected_receipts_match(self):
        result=validate_selected(self.replay,self.cover,self.result)
        self.assertEqual(result["recorded_selected_accepted_leaves"],1346)
        self.assertFalse(result["fresh_numerical_evaluation_performed_by_this_check"])

    def test_duplicate_band_rejected(self):
        replay=deepcopy(self.replay);replay["bands"][1]=deepcopy(replay["bands"][0])
        with self.assertRaises(CertificateError):validate_selected(replay,self.cover,self.result)

    def test_wrong_leaf_chain_rejected(self):
        replay=deepcopy(self.replay);replay["bands"][0]["origin_leaf_chain_sha256"]="0"*64
        with self.assertRaises(CertificateError):validate_selected(replay,self.cover,self.result)

    def test_full_replay_promotion_rejected(self):
        replay=deepcopy(self.replay);replay["whole_union_numerically_replayed"]=True
        with self.assertRaises(CertificateError):validate_selected(replay,self.cover,self.result)

    def test_gaussian_node_promotion_rejected(self):
        replay=deepcopy(self.replay);replay["gaussian_node_replay_claimed"]=True
        with self.assertRaises(CertificateError):validate_selected(replay,self.cover,self.result)

    def test_wrong_cover_rejected(self):
        replay=deepcopy(self.replay);replay["cover_sha256"]="0"*64
        with self.assertRaises(CertificateError):validate_selected(replay,self.cover,self.result)


if __name__ == "__main__":unittest.main()
