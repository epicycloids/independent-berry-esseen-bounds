"""A failed interval recheck must stop the actual certificate consumer."""
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

STRUCTURE = Path(__file__).resolve().parent / 'structure'
spec = importlib.util.spec_from_file_location('be_structure_check', STRUCTURE / 'check.py')
check = importlib.util.module_from_spec(spec)
sys.path.insert(0, str(STRUCTURE))
try:
    spec.loader.exec_module(check)
finally:
    sys.path.pop(0)


class ArithmeticFailureTests(unittest.TestCase):
    def test_failed_descendant_rejected(self):
        with patch.object(check.arbitrary, 'recheck', return_value=(False, None)) as recheck:
            with self.assertRaisesRegex(ValueError, 'Independent arithmetic failed'):
                check.arbitrary_checks(arithmetic=True)
            recheck.assert_called_once()


if __name__ == '__main__':
    unittest.main()
