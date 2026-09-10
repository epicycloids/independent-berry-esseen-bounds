from pathlib import Path
import sys
import time
HERE = Path(__file__).resolve().parent
if (HERE / 'scalar_table').is_dir():
    sys.path.insert(0, str(HERE / 'scalar_table'))
from new_prefactor_pilot import CutoffBatchWeights, CUTOFFS, maxwell_cells, original_product, LatestWeights, combined_cells
from scalar_table_bounds import table_cells, load_table
import stable_bounds
import numpy as np

def combined_table_cells(*args):
    return np.minimum(table_cells(*args), combined_cells(*args))

class TableWeights(LatestWeights):
    prefactor = staticmethod(combined_table_cells)
