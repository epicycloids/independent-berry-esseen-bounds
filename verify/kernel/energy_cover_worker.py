from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib
import importlib.util
from pathlib import Path
import time
REVIEWED_SHA256 = 'f1725128c15ba70b0e51029293e22b1a1e77c589fd907adea196a3fc41821c46'

@lru_cache(maxsize=1)
def reviewed():
    try:
        module = importlib.import_module('energy_majorant')
    except ModuleNotFoundError as error:
        if error.name != 'energy_majorant':
            raise
        path = Path(__file__).with_name('majorant.py')
        spec = importlib.util.spec_from_file_location('author_reviewed_energy_majorant', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    actual = hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
    if actual != REVIEWED_SHA256:
        raise AssertionError('The scalar kernel differs from the reviewed source')
    return module

def grid_jobs(denominator=128, price=Fraction(9, 8), descending=True):
    denominator = int(denominator)
    price = Fraction(price)
    if denominator < 2 or denominator % 2 or price not in (1, Fraction(9, 8)):
        raise ValueError('Require an even denominator and price 1 or 9/8')
    indices = range(denominator // 2, denominator)
    if descending:
        indices = reversed(tuple(indices))
    return [dict(index=i, denominator=denominator, xlo=str(Fraction(i, denominator)), xhi=str(Fraction(i + 1, denominator)), price=str(price)) for i in indices]

def run_cell(index, denominator=128, *, price=Fraction(9, 8), cpu_seconds=15.0, max_cells=400000, retain_boxes=True, precision_bits=128):
    index, denominator, price = (int(index), int(denominator), Fraction(price))
    if denominator < 2 or denominator % 2 or (not denominator // 2 <= index < denominator) or (price not in (1, Fraction(9, 8))) or (cpu_seconds <= 0):
        raise ValueError('Invalid grid cell, price, or CPU budget')
    module = reviewed()
    result = module.certify(Fraction(index, denominator), Fraction(index + 1, denominator), price=price, cpu_seconds=cpu_seconds, max_cells=int(max_cells), retain_boxes=retain_boxes, precision_bits=int(precision_bits))
    result.update(index=index, denominator=denominator, source_sha256=REVIEWED_SHA256, scalar_target='Energy-majorant gap defined by energy_majorant.gap_quotient', saved_partition=bool(retain_boxes), dispatch_performed=False)
    return result

def run_cells(indices, denominator=128, *, price=Fraction(9, 8), cpu_seconds=110.0, per_cell_seconds=15.0, max_cells_per_cell=400000, retain_boxes=True, precision_bits=128):
    reviewed()
    began = time.process_time()
    indices, rows = (tuple(indices), [])
    for index in indices:
        remaining = cpu_seconds - (time.process_time() - began)
        if remaining <= 0:
            break
        rows.append(run_cell(index, denominator, price=price, cpu_seconds=min(per_cell_seconds, remaining), max_cells=max_cells_per_cell, retain_boxes=retain_boxes, precision_bits=precision_bits))
    return dict(status='energy frequency-cell calculations', rows=rows, requested_indices=list(indices), attempted_cells=len(rows), completed_cells=sum((row['complete'] for row in rows)), unattempted_indices=list(indices[len(rows):]), denominator=int(denominator), price=str(Fraction(price)), cpu_seconds=time.process_time() - began, source_sha256=REVIEWED_SHA256, dispatch_performed=False)
