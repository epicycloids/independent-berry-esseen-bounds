import time
import numpy as np
from flint import arb
import stable_bounds
from interval_bounds import K, KHI, al, au, up, down, ia, sqpos, cubepos, m_upper, exp_upper_from_nonnegative_lower
from gaussian_batch import GaussianBatchWeights
ORIGINAL = stable_bounds.log_range_cells
ORDER = 4
COEFFICIENT = {p: al(arb(2) ** (p - 1) / p) for p in range(2, 9)}

def power_lower(x, p):
    value = np.ones_like(x)
    for _ in range(p):
        value = np.maximum(0, down(value * x))
    return value

def power_upper(x, p):
    value = np.ones_like(x)
    for _ in range(p):
        value = up(value * x)
    return value

def higher_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, order=None):
    old = ORIGINAL(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted)
    if Vlo <= 0:
        return old
    order = ORDER if order is None else order
    t2 = sqpos(ia(tlo))[0]
    k = up(KHI * cubepos(ia(thi))[1])
    A = np.maximum(0, down(down(t2 / 2) - up(up(2 * k) * au(arb(dhi).sqrt()))))
    tau = min(Bhi, max(0, taulo))
    excess = max(0, au(arb(Bhi) - arb(tau)))
    m = m_upper(thi, dhi) if deleted else np.zeros_like(thi)
    penalty = np.zeros_like(thi)
    for p in range(2, order + 1):
        coefficient = max(0, al(arb(tau) ** (2 - arb(2) / p) / arb(Vhi) ** (2 - arb(3) / p)))
        norm = np.maximum(0, down(down(A * coefficient) - up(k * excess)))
        Q = power_lower(norm, p)
        charge = np.maximum(0, down(Q - power_upper(m, p))) if deleted else Q
        penalty = np.maximum(0, down(penalty + down(COEFFICIENT[p] * charge)))
    exponent = up(up(up(k * au(arb(Bhi) + arb(tauhi))) + m) - down(Vlo / 2 * t2))
    exponent = up(exponent - penalty)
    return np.minimum(old, exp_upper_from_nonnegative_lower(np.maximum(0, -exponent)))

def install(order=4):
    global ORDER
    assert order in [4, 8]
    ORDER = order
    stable_bounds.log_range_cells = higher_range
    import stable_cover
    if getattr(stable_cover, '_higher_log_installed', False):
        return
    original_factory = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind != 'higher':
            return original_factory(hi, N, kind, target)
        cache = {}

        def choose(tauhi):
            bucket = min(20, max(1, round(20 * tauhi / hi)))
            if bucket not in cache:
                cache[bucket] = GaussianBatchWeights(hi, N, hi * bucket / 20)
            return cache[bucket]
        return choose
    stable_cover.weight_factory = factory
    stable_cover._higher_log_installed = True

def checks():
    rng = np.random.default_rng(260938)
    count = 0
    for _ in range(40):
        V = float(rng.uniform(0.01, 1))
        B = float(rng.uniform(0.01, 2))
        d = float(rng.uniform(0.001, V))
        tl = float(rng.uniform(0, B))
        th = float(rng.uniform(tl, B))
        a = float(rng.uniform(0.001, 5))
        b = a + float(rng.uniform(0, 0.02))
        for order in [4, 8]:
            for deleted in [False, True]:
                y = float(higher_range(np.array([a]), np.array([b]), V, V, B, d, tl, th, deleted, order)[0])
                for tt in [a, (a + b) / 2, b]:
                    t = arb(tt)
                    k = K * t ** 3
                    A = t * t / 2 - 2 * k * arb(d).sqrt()
                    z = min(t * arb(d).sqrt(), 1 / (6 * K))
                    m = z * z / 2 - 2 * K * z ** 3 if deleted else arb(0)
                    for tau in [tl, (tl + th) / 2, th]:
                        s = arb(tau)
                        penalty = arb(0)
                        for p in range(2, order + 1):
                            Q = max(arb(0), A * s ** (2 - arb(2) / p) / arb(V) ** (2 - arb(3) / p) - k * (arb(B) - s)) ** p if A >= 0 else arb(0)
                            penalty += arb(2) ** (p - 1) / p * max(arb(0), Q - m ** p)
                        ref = (-arb(V) * t * t / 2 + k * (arb(B) + s) + m - penalty).exp()
                        assert y == 1 or arb(y) >= ref, (V, B, d, tl, th, tt, tau, order, deleted, y, ref)
                        count += 1
    return {'higher_log_cell_enclosures': count, 'all_passed': True}
