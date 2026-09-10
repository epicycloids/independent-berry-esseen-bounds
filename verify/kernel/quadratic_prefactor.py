from pathlib import Path
import importlib.util
import sys
HERE = Path(__file__).resolve().parent
if (HERE.parent / 'odd_refinement/odd_bounds.py').is_file():
    spec = importlib.util.spec_from_file_location('odd_prefactor', HERE.parent / 'odd_refinement/odd_bounds.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules['odd_prefactor'] = module
    spec.loader.exec_module(module)
from odd_prefactor import odd_cells
import numpy as np
from flint import arb, ctx
from interval_bounds import up, au, ia, sqpos
from paired_bounds import rational_poly, alpha_upper

def quadratic_cells(thi, Bhi, dhi, taulo, tauhi):
    tau = min(Bhi, max(0.0, taulo))
    x = up(np.asarray(thi) * au(arb(dhi).sqrt()))
    z = sqpos(ia(np.minimum(x, 2.0)))
    p = rational_poly(z, [arb(1), arb(1) / 6, -arb(1) / 40, arb(1) / 1008])[1]
    q = up(1 + up(up(x * x) / 3))
    C = up(1 + up(q * q))
    alpha = alpha_upper(x)
    D = np.maximum(up(up(8 * up(p * p)) / 3), up(up(2 * alpha) * up(np.sqrt(C))))
    excess = max(0.0, au(arb(Bhi) - arb(tau)))
    first = up(alpha * tau)
    linear = up(up(D * tau) * excess)
    quadratic = up(C * up(excess * excess))
    square = up(up(up(first * first) + linear) + quadratic)
    candidate = np.minimum(Bhi, up(np.sqrt(np.maximum(0.0, square))))
    candidate = np.where(x <= 2.0, candidate, Bhi)
    return np.minimum(odd_cells(thi, Bhi, dhi, taulo, tauhi), candidate)

def exact(t, B, d, tau):
    t, B, d, tau = map(arb, (t, B, d, tau))
    x = t * d.sqrt()
    if x > 2:
        return B
    alpha = 2 * x / 3 - x ** 3 / 15 + x ** 5 / 420
    p = 1 + x * x / 6 - x ** 4 / 40 + x ** 6 / 1008
    q = 1 + x * x / 3
    C = 1 + q * q
    D = max(8 * p * p / 3, 2 * alpha * C.sqrt())
    return min(B, (alpha ** 2 * tau ** 2 + D * tau * (B - tau) + C * (B - tau) ** 2).sqrt())

def checks():
    previous = ctx.prec
    ctx.prec = 192
    count, improved = (0, 0)
    try:
        ts = np.array([0.0, 0.01, 0.1, 0.5, 1.0, 1.5, 1.99, 2.0, 2.01, 3.0, 10.0])
        for B, d, lo, hi in [(0.5, 0.25, 0.35, 0.45), (0.5, 1.0, 0.45, 0.5), (0.5, 0.04, 0.0, 0.5), (0.5, 0.25, 0.5, 0.5), (0.5, 0.0, 0.35, 0.5), (1.0, 1.0, 0.1, 0.2), (0.6, 0.4, 0.5, 0.55)]:
            values = quadratic_cells(ts, B, d, lo, hi)
            old = odd_cells(ts, B, d, lo, hi)
            assert np.all(values <= old)
            improved += int(np.count_nonzero(values < old))
            for t, value in zip(ts, values):
                for fraction in (0.0, 0.5, 1.0):
                    for tau in (lo, (lo + hi) / 2, hi):
                        assert arb(float(value)) >= exact(t * fraction, B, d, tau)
                        count += 1
        return dict(all_passed=True, interval_comparisons=count, strict_improvements=improved)
    finally:
        ctx.prec = previous
