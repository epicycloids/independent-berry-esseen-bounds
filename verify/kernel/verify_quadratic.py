import time
import numpy as np
from flint import arb, ctx
from quadratic_bounds import quadratic_cells, exact_prefactor, exact_majorant, previous_cells
from concave_bounds import exact_prefactor as previous_exact

def quick_checks():
    began = time.monotonic()
    old_precision = ctx.prec
    ctx.prec = 180
    count, strict, general_checks = (0, 0, 0)
    try:
        for a, b, u, v in [(0.3, 2.0, 1.0, 2.0), (0.4, 3.0, 1.0, 6.0), (0.5, 8.0, 9.0, 10.0), (0.5, 0.3, 9.0, 10.0), (0.2, 1.0, 0.0, 3.0), (0.0, 2.0, 1.0, 5.0)]:
            for e in [0.0, 0.01, 0.1, 0.3, 1.0, 3.0, 100.0]:
                value = exact_majorant(a, b, u, v, e, capped=False)
                aa, bb, uu, vv, ee = map(arb, (a, b, u, v, e))
                first = (aa * aa + bb * ee + uu * ee * ee).sqrt()
                second = (bb * ee + vv * ee * ee).sqrt()
                target = max(first, second)
                assert float(value) + 2e-13 >= float(target)
                general_checks += 1
        ts = np.array([0.0, 0.0001, 0.01, 0.1, 0.5, 1.0, 1.5, 1.99, 2.0, 2.01, 3.0, 10.0])
        for B, d, tl, th in [(0.5, 0.25, 0.35, 0.45), (0.5, 1.0, 0.45, 0.5), (0.5, 0.04, 0.0, 0.5), (0.5, 0.25, 0.5, 0.5), (0.5, 0.0, 0.35, 0.5), (0.6, 0.4, 0.52, 0.56)]:
            bound = quadratic_cells(ts, B, d, tl, th)
            previous = previous_cells(ts, B, d, tl, th)
            assert np.all(bound <= previous)
            strict += int(np.count_nonzero(bound < previous))
            for t, upper in zip(ts, bound):
                for fraction in [0.0, 0.5, 1.0]:
                    for tau in [tl, (tl + th) / 2, th]:
                        reference = exact_prefactor(t * fraction, B, d, tau)
                        if not arb(float(upper)) >= reference:
                            raise AssertionError((t, B, d, tau, upper, reference))
                        count += 1
        rows = []
        for x in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]:
            for e in [0.05, 0.1, 0.2, 0.5]:
                old = previous_exact(x, 1 + arb(e), 1, 1)
                new = exact_prefactor(x, 1 + arb(e), 1, 1)
                assert float(new) <= float(old) + 2e-13
                rows.append(dict(x=x, excess=e, previous=float(old), quadratic=float(new), improvement=float(old - new)))
        return dict(all_passed=True, general_coefficient_comparisons=general_checks, outward_whole_cell_comparisons=count, strict_array_entries=strict, scalar_rows=rows, elapsed_seconds=time.monotonic() - began, status='implementation diagnostics conditional on the stronger scalar theorem')
    finally:
        ctx.prec = old_precision
