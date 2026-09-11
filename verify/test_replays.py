"""Check that the numerical settings retain their source-specific meaning."""
from copy import deepcopy
import unittest
from check_saved import CertificateError,load_inputs,validate_evaluators

class ReplaySettingsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.cover,cls.result=load_inputs()
    @classmethod
    def tearDownClass(cls):
        del cls.cover,cls.result
    def test_profiles_are_valid(self):validate_evaluators(self.cover['evaluators'])
    def test_target_relabelling_rejected(self):
        values=deepcopy(self.cover['evaluators']);key=next(k for k,v in values.items()if v['kind']=='signed' and v['stop_target']=='0.450')
        values[key]['stop_target']='0.44995'
        with self.assertRaises(CertificateError):validate_evaluators(values)
    def test_coupled_recipe_rejected(self):
        values=deepcopy(self.cover['evaluators']);key=next(k for k,v in values.items()if v['kind']=='signed');values[key]['coupling']=True
        with self.assertRaises(CertificateError):validate_evaluators(values)
if __name__=='__main__':unittest.main()
