"""Failure controls for the current closed moment partition."""
from copy import deepcopy
from fractions import Fraction
import unittest
import zipfile
from check_saved import ROOT,CertificateError,UPPER,load_inputs,load_band,validate_band,validate_evaluators

class SavedCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cover,cls.result=load_inputs()
        with zipfile.ZipFile(ROOT/cls.cover['band_archive']) as archive:cls.band=load_band(archive,cls.cover['bands'][0])
    @classmethod
    def tearDownClass(cls):
        del cls.cover,cls.result,cls.band
    def test_exact_upper(self):
        self.assertEqual(UPPER,Fraction(self.result['largest_recorded_upper']));self.assertLess(UPPER,Fraction('0.44995'))
    def test_complete_tree(self):
        local=validate_band(self.band,self.cover['evaluators']);self.assertEqual(local['accepted'],len(self.band['leaves']))
    def test_root_cannot_shrink(self):
        band=deepcopy(self.band);band['root'][1]*=.99
        with self.assertRaises(CertificateError):validate_band(band,self.cover['evaluators'])
    def test_tree_cannot_end_early(self):
        band=deepcopy(self.band);band['tree']=band['tree'][:-1]
        with self.assertRaises(CertificateError):validate_band(band,self.cover['evaluators'])
    def test_leaf_box_cannot_change(self):
        band=deepcopy(self.band);band['leaves'][0]['box'][1]*=.99
        with self.assertRaises(CertificateError):validate_band(band,self.cover['evaluators'])
    def test_underlying_replay_scope_cannot_change(self):
        band=deepcopy(self.band);band['leaves'][0]['arithmetic_replayed']=False
        with self.assertRaises(CertificateError):validate_band(band,self.cover['evaluators'])
    def test_evaluator_depth_is_fixed(self):
        values=deepcopy(self.cover['evaluators']);key=next(k for k,v in values.items()if v['kind']=='signed');values[key]['b_max_depth']=7
        with self.assertRaises(CertificateError):validate_evaluators(values)
    def test_cap_above_global_upper_rejected(self):
        band=deepcopy(self.band);band['leaves'][0]['upper']=.45
        with self.assertRaises(CertificateError):validate_band(band,self.cover['evaluators'])
if __name__=='__main__':unittest.main()
