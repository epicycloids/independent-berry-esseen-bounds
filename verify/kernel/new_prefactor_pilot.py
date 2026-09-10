from pathlib import Path
import sys
import time
HERE = Path(__file__).resolve().parent
for p in ('capped_concavity', 'variance_resolved', 'odd_quartic'):
    if (HERE / p).is_dir():
        sys.path.insert(0, str(HERE / p))
from further_bounds_pilot import CutoffBatchWeights, CUTOFFS, maxwell_cells, original_product, LatestWeights
from capped_bounds import capped_cells
from variance_bounds import variance_cells
from quartic_bounds import quartic_cells
import stable_bounds
import numpy as np

def combined_cells(*args):
    return np.minimum.reduce([capped_cells(*args), variance_cells(*args), quartic_cells(*args)])

class CombinedWeights(LatestWeights):
    prefactor = staticmethod(combined_cells)
