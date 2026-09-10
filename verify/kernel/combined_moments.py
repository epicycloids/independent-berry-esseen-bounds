from pathlib import Path
import sys
import time
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
for path in (HERE.parent / '26', HERE / 'logarithm', HERE / 'entropy_moments'):
    if path.exists():
        sys.path.insert(0, str(path))
import stable_bounds
from interval_bounds import K, KLO, KHI, au, up, down, sqpos, cubepos, m_upper, exp_upper_from_nonnegative_lower
from moment_fast import square_sum_lower, tangent_charge
from weighted_log_bounds import weighted_charge_lower
from entropy_charge import discrete_charge_lower
ORIGINAL = stable_bounds.log_range_cells

def joint_charge(Q, S, m, mode='full', joint_delete=True):
    Qdel = np.maximum(0.0, down(Q - up(m * m)))
    Sdel = up(S - m) if joint_delete else S
    if mode == 'continuous':
        return weighted_charge_lower(Qdel, Sdel)
    return discrete_charge_lower(Qdel, Sdel, mode=mode)

def combined_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, mode='full', joint_delete=True):
    old = ORIGINAL(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted)
    if Vlo <= 0:
        return old
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    Q = square_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo)
    m = m_upper(thi, dhi) if deleted else np.zeros_like(thi)
    t2lo, t2hi = sqpos((tlo, thi))
    t3lo, t3hi = cubepos((tlo, thi))
    khi = up(KHI * t3hi)
    tau = min(Bhi, max(0.0, taulo))
    S = np.maximum(0.0, up(up(up(Vhi / 2) * t2hi) - down(down(down(2 * KLO) * t3lo) * tau)))
    penalty = np.maximum(tangent_charge(Q, S, m), joint_charge(Q, S, m, mode, joint_delete))
    positive = up(khi * au(arb(Bhi) + arb(tauhi)))
    negative = down(down(Vlo / 2) * t2lo)
    exponent = up(up(up(positive + m) - negative) - penalty)
    return np.minimum(old, exp_upper_from_nonnegative_lower(np.maximum(0.0, -exponent)))

def install(mode='full', joint_delete=True):
    stable_bounds.log_range_cells = lambda *args, **kwargs: combined_range(*args, **kwargs, mode=mode, joint_delete=joint_delete)

def checks():
    rng = np.random.default_rng(270931)
    products = joint = 0
    for n in [1, 2, 3, 8, 25]:
        for _ in range(10):
            v = rng.dirichlet(np.ones(n))
            b = v ** 1.5 * (1 + rng.uniform(0.0, 0.6, n))
            V = float(v.sum())
            B = float(b.sum())
            d = float(v.max())
            tau = float(np.sum(v ** 1.5))
            ts = np.linspace(0.0001, 5.0, 32)
            eps = 2e-12
            q = np.maximum(0.0, ts[:, None] ** 2 * v / 2 - float(K) * ts[:, None] ** 3 * (b + v ** 1.5))
            factors = np.sqrt(np.maximum(0.0, 1 - 2 * q))
            for deleted in [False, True]:
                actual = np.max([np.prod(np.delete(factors, j, axis=1), axis=1) for j in range(n)], axis=0) if deleted else np.prod(factors, axis=1)
                for mode in ['continuous', 'full']:
                    bound = combined_range(ts, ts, V - eps, V + eps, B + eps, d + eps, tau - eps, tau + eps, deleted, mode)
                    assert np.all(actual <= bound + 2e-12), (n, deleted, mode)
                    products += len(ts)
    for n in [1, 2, 3, 8, 25]:
        for _ in range(24):
            q = rng.uniform(0.0001, 0.499, n)
            H = q.sum() - rng.uniform(0.0, 0.3)
            S = q.sum() + rng.uniform(0.0, 0.3)
            Q = np.sum(q * q) * rng.uniform(0.5, 1.0)
            m = rng.uniform(q.max(), 0.4999)
            actual = np.max([np.prod(np.sqrt(1 - 2 * np.delete(q, j))) for j in range(n)])
            for mode in ['continuous', 'full']:
                penalty = joint_charge(np.array([Q]), np.array([S]), np.array([m]), mode)
                exponent = up(up(m - H) - penalty)
                bound = exp_upper_from_nonnegative_lower(np.maximum(0.0, -exponent))[0]
                assert actual <= bound + 2e-12, (n, mode, q, H, S, Q, m, actual, bound)
                joint += 1
    return dict(finite_product_comparisons=products, joint_deleted_comparisons=joint, all_passed=True)
