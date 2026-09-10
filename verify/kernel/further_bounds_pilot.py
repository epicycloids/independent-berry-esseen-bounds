from pathlib import Path
import sys
import time
HERE = Path(__file__).resolve().parent
for path in (HERE / 'joint_deletion', HERE / 'odd_quadratic', HERE / 'concave_majorant', HERE / 'concave_majorant/quadratic_variant'):
    if path.is_dir():
        sys.path.insert(0, str(path))
from odd_cutoff_pilot import CUTOFFS, CutoffBatchWeights
from odd_prefactor import odd_cells
from dual_log_bounds import dual_range
from joint_deletion_bounds import joint_range
from quadratic_prefactor import quadratic_cells as polynomial_cells
from concave_bounds import concave_cells
from quadratic_bounds import quadratic_cells as maxwell_cells
from paired_bounds import CosineWeights
from quadratic_certificate import certificate
from verify_quadratic import quick_checks
import stable_bounds
import numpy as np

def original_product(*a, **kw):
    return dual_range(*a, **kw, baseline='original')

def joint_product(*a, **kw):
    return joint_range(*a, **kw, baseline='original')

class LatestWeights(CosineWeights):
    prefactor = staticmethod(maxwell_cells)
