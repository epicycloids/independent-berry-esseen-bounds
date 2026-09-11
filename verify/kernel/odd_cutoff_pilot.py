from pathlib import Path
import importlib.util
import sys
import time
HERE = Path(__file__).resolve().parent
for path in (HERE.parent / '26', HERE / 'odd_refinement', HERE / 'dual_logarithm', HERE / 'smoothing_geometry', HERE / 'smoothing_geometry/adjoint', HERE / 'smoothing_geometry/adjoint/shared_cutoff'):
    if path.exists():
        sys.path.insert(0, str(path))
import numpy as np
import stable_bounds
from cutoff_transfer import CutoffBatchWeights
from dual_log_bounds import dual_range
if (HERE / 'odd_refinement/odd_bounds.py').is_file():
    spec = importlib.util.spec_from_file_location('odd_prefactor', HERE / 'odd_refinement/odd_bounds.py')
    odd_prefactor = importlib.util.module_from_spec(spec)
    sys.modules['odd_prefactor'] = odd_prefactor
    spec.loader.exec_module(odd_prefactor)
from odd_prefactor import odd_cells, OddWeights, checks
from paired_bounds import cosine_cells
CUTOFFS = (0.75, 0.78125, 0.8125, 0.84375, 0.875, 0.90625, 0.9375, 0.96875, 1.0)

class OddCutoffWeights(CutoffBatchWeights):
    prefactor = staticmethod(odd_cells)

def broad_probe(Ls, N=256, max_seconds=180.0):
    began = time.monotonic()
    old = stable_bounds.log_range_cells
    stable_bounds.log_range_cells = lambda *a, **kw: dual_range(*a, **kw, baseline='original')
    records, evaluations = ([], 0)
    stopped = False
    try:
        for L in Ls:
            cache, rows = ({}, [])
            for fraction in (0.55, 0.65, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0):
                tau = min(1.0, L) * fraction
                bucket = min(20, max(1, round(20 * tau / L)))
                if bucket not in cache:
                    cache[bucket] = OddCutoffWeights(L, N, L * bucket / 20, cutoff_fractions=CUTOFFS)
                weights = cache[bucket]
                for position in (0.0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0):
                    d = tau * tau + position * (tau ** (2 / 3) - tau * tau)
                    for skew in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
                        b = d ** 1.5 + skew * (L - tau)
                        e = 2e-12
                        value = weights.box(L - e, L + e, d - e, d + e, b - e, b + e, tau - e, tau + e)
                        evaluations += 1
                        if value is not None:
                            rows.append(dict(L=L, d=d, b=b, tau=tau, upper=value, N=N, fraction=fraction, position=position, skew=skew))
                        if time.monotonic() - began >= max_seconds:
                            stopped = True
                            break
                    if stopped:
                        break
                if stopped:
                    break
            largest = sorted(rows, key=lambda row: row['upper'], reverse=True)[:10]
            refined = []
            if not stopped:
                for row in largest[:4]:
                    cap = L * min(20, max(1, round(20 * row['tau'] / L))) / 20
                    w = OddCutoffWeights(L, 1024, cap, cutoff_fractions=CUTOFFS)
                    e = 2e-12
                    d, b, tau = (row[k] for k in ('d', 'b', 'tau'))
                    refined.append(dict(row, N=1024, upper=w.box(L - e, L + e, d - e, d + e, b - e, b + e, tau - e, tau + e)))
            records.append(dict(L=L, grid_completed=not stopped, evaluated=len(rows), largest=largest, refined=refined))
            print(dict(L=L, evaluated=len(rows), largest=largest[:1]), flush=True)
            if stopped:
                break
    finally:
        stable_bounds.log_range_cells = old
    return dict(status='Interval enclosures around sampled parameter points', rows=records, complete_grid=not stopped, evaluations=evaluations, elapsed_seconds=time.monotonic() - began)
