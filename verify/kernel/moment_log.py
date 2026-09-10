import math
import time
import numpy as np
from flint import arb
import stable_bounds
from interval_bounds import K, KHI, al, au, up, down, ia, sqpos, cubepos, m_upper, exp_upper_from_nonnegative_lower
ORIGINAL = stable_bounds.log_range_cells
SUPPORTS = [arb(1) / 10000, arb(1) / 1000, arb(3) / 1000]
SUPPORTS += [arb(j) / 100 for j in range(1, 22)]
ALPHA = np.array([al(2 * u * (1 - u) * (1 - 2 * u) / (16 * K)) for u in SUPPORTS])
BETA = np.array([au(u * u * (1 - u) * (1 - 3 * u) / (64 * K * K)) for u in SUPPORTS])
LIMIT = np.array([al(1 - u - (2 * u * (1 - u)).sqrt()) for u in SUPPORTS])

def square_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo):
    tlo, thi = (np.asarray(tlo), np.asarray(thi))
    if Vhi <= 0 or dhi <= 0:
        return np.zeros_like(tlo)
    tau = min(Bhi, max(0.0, taulo))
    positive = down(cubepos(ia(tlo))[0][..., None] * down(ALPHA * tau))
    negative = up(sqpos(ia(thi))[1][..., None] * up(BETA * Vhi))
    support = down(positive - negative)
    argument = up(up(au(4 * K) * thi) * au(arb(dhi).sqrt()))
    valid = argument[..., None] <= LIMIT
    r2 = np.maximum(0.0, np.max(np.where(valid, support, 0.0), axis=-1))
    norm = np.maximum(0.0, down(np.sqrt(r2)))
    excess = max(0.0, au(arb(Bhi) - arb(tau)))
    k = up(KHI * cubepos(ia(thi))[1])
    q = np.maximum(0.0, down(norm - up(k * excess)))
    return np.maximum(0.0, down(q * q))

def moment_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False):
    old = ORIGINAL(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted)
    if Vlo <= 0:
        return old
    penalty = square_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo)
    m = m_upper(thi, dhi) if deleted else np.zeros_like(thi)
    if deleted:
        penalty = np.maximum(0.0, down(penalty - up(m * m)))
    k = up(KHI * cubepos(ia(thi))[1])
    positive = up(k * au(arb(Bhi) + arb(tauhi)))
    negative = down(down(Vlo / 2) * sqpos(ia(tlo))[0])
    exponent = up(up(up(positive + m) - negative) - penalty)
    extra = exp_upper_from_nonnegative_lower(np.maximum(0.0, -exponent))
    return np.minimum(old, extra)

def install():
    stable_bounds.log_range_cells = moment_range

def checks():
    rng = np.random.default_rng(270901)
    scalar = finite = 0
    for _ in range(80):
        V = float(rng.uniform(0.03, 1.0))
        B = float(rng.uniform(0.02, 1.5))
        a = float(rng.uniform(0.01, V))
        tau = float(rng.uniform(0.0, B))
        tlo = float(rng.uniform(0.01, 5.0))
        thi = tlo + float(rng.uniform(0.0, 0.03))
        lower = float(square_sum_lower(np.array([tlo]), np.array([thi]), V, B, a, tau)[0])
        for tt in [tlo, (tlo + thi) / 2, thi]:
            t = arb(tt)
            candidates = [arb(0)]
            for u in SUPPORTS:
                limit = 1 - u - (2 * u * (1 - u)).sqrt()
                if 4 * K * t * arb(a).sqrt() <= limit:
                    alpha = 2 * u * (1 - u) * (1 - 2 * u)
                    beta = u * u * (1 - u) * (1 - 3 * u)
                    r2 = alpha * t ** 3 * arb(tau) / (16 * K) - beta * t * t * arb(V) / (64 * K * K)
                    norm = max(arb(0), r2).sqrt()
                    candidates.append(max(arb(0), norm - K * t ** 3 * (arb(B) - arb(tau))) ** 2)
            ref = max(candidates)
            assert lower == 0 or arb(lower) <= ref, (lower, ref)
            scalar += 1
    for n in [1, 2, 3, 8, 25]:
        for _ in range(15):
            v = rng.dirichlet(np.ones(n))
            b = v ** 1.5 * (1 + rng.uniform(0.0, 0.5, n))
            V = float(v.sum())
            B = float(b.sum())
            a = float(v.max())
            tau = float(np.sum(v ** 1.5))
            ts = np.linspace(0.001, 5.0, 32)
            eps = 2e-12
            qlow = square_sum_lower(ts, ts, V + eps, B + eps, a + eps, max(0.0, tau - eps))
            q = np.maximum(0.0, ts[:, None] ** 2 * v / 2 - float(K) * ts[:, None] ** 3 * (b + v ** 1.5))
            assert np.all(qlow <= np.sum(q * q, axis=1) + 2e-12)
            for deleted in [False, True]:
                bound = moment_range(ts, ts, V - eps, V + eps, B + eps, a + eps, max(0.0, tau - eps), tau + eps, deleted)
                factors = np.sqrt(np.maximum(0.0, 1 - 2 * q))
                if deleted:
                    actual = np.max([np.prod(np.delete(factors, j, axis=1), axis=1) for j in range(n)], axis=0)
                else:
                    actual = np.prod(factors, axis=1)
                assert np.all(actual <= bound + 2e-12), (n, deleted)
                finite += len(ts)
    return {'scalar_comparisons': scalar, 'finite_product_comparisons': finite, 'all_passed': True}
