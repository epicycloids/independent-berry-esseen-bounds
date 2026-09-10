"""Regression checks for certificate coverage and metadata."""
from copy import deepcopy
from fractions import Fraction
import json
import unittest

from check_saved import ROOT, CertificateError, canonical, load_inputs, mixed_record, validate_cover


class SavedCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cover, cls.result = load_inputs()
        cls.mixed = next(b["record"] for b in cls.cover["bands"] if b["type"] == "mixed")

    def test_candidate_target_is_exact(self):
        upper = Fraction(self.result["largest_recorded_upper_fraction"])
        self.assertLess(upper,Fraction("0.454"))
        self.assertEqual(upper,Fraction(self.result["largest_recorded_upper"]))

    def test_curated_source_identity_loaders(self):
        import adaptive_b_cover_engine
        import energy_reference
        import energy_small_table
        import maximum_variance_bounds
        adaptive_b_cover_engine.install()
        self.assertIs(energy_reference.reviewed_reference(),maximum_variance_bounds.reviewed_reference())
        self.assertEqual(energy_small_table.worker().REMAINDER_SHA256,
                         energy_small_table.REMAINDER_SHA256)
        table = json.loads((ROOT/"verify/kernel/energy_table_complete.json").read_text())
        self.assertEqual(table["source_sha256"],energy_reference.REVIEWED_SHA256)
        self.assertEqual(energy_small_table._compiled_endpoint(table["small_frequency_prefix"]),Fraction(1,2))

    def test_inflated_review_status_rejected(self):
        result = dict(self.result,independent_mathematical_review_completed=True)
        with self.assertRaises(CertificateError):validate_cover(self.cover,result)

    def test_unsupported_full_replay_rejected(self):
        result = dict(self.result,all_leaf_quadratures_replayed=True)
        with self.assertRaises(CertificateError):validate_cover(self.cover,result)

    def test_promoted_global_proof_rejected(self):
        result = dict(self.result,global_proof=True)
        with self.assertRaises(CertificateError):validate_cover(self.cover,result)

    def test_old_target_cannot_be_relabelled(self):
        record = deepcopy(self.mixed)
        record["origin"]["target"] = "0.454"
        with self.assertRaises(CertificateError):mixed_record(record)

    def test_frontier_must_be_original(self):
        record = deepcopy(self.mixed)
        record["origin"]["stack"].pop()
        with self.assertRaises((CertificateError,ValueError)):mixed_record(record)

    def test_closed_root_cannot_shrink(self):
        record = deepcopy(self.mixed)
        record["root"][1] *= .99
        with self.assertRaises(CertificateError):mixed_record(record)

    def test_empty_continuation_cannot_claim_completion(self):
        import stable_cover
        import hashlib
        record = deepcopy(self.mixed)
        record["continuation"].update(tree_zlib_base64=stable_cover.pack(""),tree_sha256=hashlib.sha256(b"").hexdigest(),nodes=0)
        with self.assertRaises(CertificateError):mixed_record(record)

    def test_mixed_topology_counts_replay(self):
        counts = mixed_record(self.mixed)
        self.assertEqual(counts,{k:self.mixed[k] for k in ("nodes","accepted","infeasible")})

    def test_wrong_b_coverage_policy_rejected(self):
        record = deepcopy(self.mixed)
        record["evaluator"]["parameters"]["b_max_depth"] = 7
        with self.assertRaises(CertificateError):mixed_record(record)

    def test_numerical_mode_requires_both_roles(self):
        with self.assertRaises(CertificateError):mixed_record(self.mixed,numerical=True)


if __name__ == "__main__":unittest.main()
