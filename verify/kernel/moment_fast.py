from pathlib import Path
import sys
import time
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
for path in (HERE.parent / '26', HERE / 'logarithm'):
    if path.exists():
        sys.path.insert(0, str(path))
import stable_bounds
from interval_bounds import K, KLO, KHI, al, au, up, down, ia, sqpos, cubepos, m_upper, exp_upper_from_nonnegative_lower
from weighted_log_bounds import weighted_range, tangent_table
ORIGINAL = stable_bounds.log_range_cells

def choose_support(tlo, thi, Vhi, dhi, tau):
    M = np.clip(4 * float(K) * thi * np.sqrt(dhi), 0.0, 1.0)
    limit = (1 - M) ** 2 / (2 - M + np.sqrt(1 + 2 * M - 2 * M * M))
    mean = np.divide(4 * float(K) * tlo ** 3 * tau, thi ** 2 * Vhi, out=np.zeros_like(tlo), where=thi > 0)
    return np.maximum(0.0, np.minimum(np.minimum(mean, limit), 0.21)) * (1 - 2 ** (-20))

def support_coefficients(u):
    one_lo = np.maximum(0.0, down(1 - u))
    alpha = np.maximum(0.0, down(down(down(2 * u) * one_lo) * np.maximum(0.0, down(1 - up(2 * u)))))
    alpha = np.maximum(0.0, down(alpha / au(16 * K)))
    beta = up(up(up(u * u) * up(1 - u)) * up(1 - down(3 * u)))
    beta = np.maximum(0.0, up(beta / al(64 * K * K)))
    radical = up(up(2 * u) * up(1 - u))
    limit = down(down(1 - u) - up(np.sqrt(np.maximum(0.0, radical))))
    return (alpha, beta, limit)

def square_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo, return_support=False):
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    if Vhi <= 0 or dhi <= 0:
        zero = np.zeros_like(tlo)
        return (zero, zero, np.zeros_like(tlo, dtype=bool)) if return_support else zero
    tau = min(Bhi, max(0.0, taulo))
    u = choose_support(tlo, thi, Vhi, dhi, tau)
    alpha, beta, limit = support_coefficients(u)
    t3lo, t3hi = cubepos((tlo, thi))
    t2hi = sqpos(ia(thi))[1]
    positive = down(t3lo * down(alpha * tau))
    negative = up(t2hi * up(beta * Vhi))
    argument = up(up(au(4 * K) * thi) * au(arb(dhi).sqrt()))
    valid = argument <= limit
    r2 = np.where(valid, np.maximum(0.0, down(positive - negative)), 0.0)
    norm = np.maximum(0.0, down(np.sqrt(r2)))
    excess = max(0.0, au(arb(Bhi) - arb(tau)))
    k = up(KHI * t3hi)
    q = np.maximum(0.0, down(norm - up(k * excess)))
    Q = np.maximum(0.0, down(q * q))
    return (Q, u, valid) if return_support else Q

def tangent_charge(Q, S, m):
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.divide(2 * np.maximum(0.0, Q - m * m), S - m, out=np.zeros_like(Q), where=S > m)
        index = np.rint(512 * np.clip(ratio, 0.0, 504 / 512)).astype(np.intp)
    aa, bb = tangent_table()
    alpha, beta = (aa[index], bb[index])
    loss = np.maximum(0.0, up(up(alpha * up(m * m)) - down(beta * m)))
    return np.maximum(0.0, down(down(down(alpha * Q) - up(beta * S)) - loss))

def moment_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, charge='tangent', baseline='weighted'):
    old = weighted_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted, method='tangent', passes=1) if baseline == 'weighted' else ORIGINAL(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted)
    if Vlo <= 0:
        return old
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    Q = square_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo)
    m = m_upper(thi, dhi) if deleted else np.zeros_like(thi)
    t2lo, t2hi = sqpos((tlo, thi))
    t3lo, t3hi = cubepos((tlo, thi))
    khi = up(KHI * t3hi)
    if charge == 'quadratic':
        penalty = np.maximum(0.0, down(Q - up(m * m)))
    else:
        tau = min(Bhi, max(0.0, taulo))
        S = np.maximum(0.0, up(up(up(Vhi / 2) * t2hi) - down(down(down(2 * KLO) * t3lo) * tau)))
        penalty = tangent_charge(Q, S, m)
    positive = up(khi * au(arb(Bhi) + arb(tauhi)))
    negative = down(down(Vlo / 2) * t2lo)
    exponent = up(up(up(positive + m) - negative) - penalty)
    return np.minimum(old, exp_upper_from_nonnegative_lower(np.maximum(0.0, -exponent)))

def install(charge='tangent', baseline='weighted'):
    stable_bounds.log_range_cells = lambda *a, **kw: moment_range(*a, **kw, charge=charge, baseline=baseline)

def checks():
    rng = np.random.default_rng(270922)
    scalar = products = 0
    for _ in range(100):
        V = float(rng.uniform(0.01, 1.0))
        B = float(rng.uniform(0.01, 1.5))
        d = float(rng.uniform(0.001, V))
        tau = float(rng.uniform(0.0, B))
        tl = float(rng.uniform(0.0, 5.0))
        th = tl + float(rng.uniform(0.0, 0.1))
        Q, us, valid = square_sum_lower(np.array([tl]), np.array([th]), V, B, d, tau, True)
        u = arb(float(us[0]))
        alpha = 2 * u * (1 - u) * (1 - 2 * u)
        beta = u * u * (1 - u) * (1 - 3 * u)
        for tt in [tl, (tl + th) / 2, th]:
            t = arb(tt)
            if valid[0]:
                assert 4 * K * t * arb(d).sqrt() <= 1 - u - (2 * u * (1 - u)).sqrt()
                R = alpha * t ** 3 * arb(tau) / (16 * K) - beta * t * t * arb(V) / (64 * K * K)
                reference = max(arb(0), max(arb(0), R).sqrt() - K * t ** 3 * (arb(B) - arb(tau))) ** 2
                assert Q[0] == 0 or arb(float(Q[0])) <= reference
            else:
                assert Q[0] == 0
            scalar += 1
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
            Q = square_sum_lower(ts, ts, V + eps, B + eps, d + eps, tau - eps)
            assert np.all(Q <= np.sum(q * q, axis=1) + 2e-12)
            factors = np.sqrt(np.maximum(0.0, 1 - 2 * q))
            for deleted in [False, True]:
                bound = moment_range(ts, ts, V - eps, V + eps, B + eps, d + eps, tau - eps, tau + eps, deleted)
                actual = np.max([np.prod(np.delete(factors, j, axis=1), axis=1) for j in range(n)], axis=0) if deleted else np.prod(factors, axis=1)
                assert np.all(actual <= bound + 2e-12), (n, deleted)
                products += len(ts)
    return {'scalar_selected_supports': scalar, 'finite_products': products, 'all_passed': True}
