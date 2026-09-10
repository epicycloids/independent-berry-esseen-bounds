import math
import time
import numpy as np
from flint import arb
from interval_bounds import Weights, K, KLO, KHI, ia, sqpos, cubepos, up, down, au, al, m_upper, exp_upper_from_nonnegative_lower, f_cells, exp_cubic_cells, elementary_checks

def root_upper(x):
    hi = x.upper()
    return 0.0 if hi <= 0 else au(hi.sqrt())

def moment_coefficients(Bhi, dhi, taulo, tauhi):
    B = arb(Bhi)
    d = arb(dhi)
    tl = arb(min(Bhi, max(0, taulo)))
    th = arb(min(Bhi, max(0, tauhi)))
    J = root_upper((B * B + B * (B * B + 8 * tl * tl).sqrt()) / 2 - 2 * tl * tl)
    mid = min(th, max(tl, B / 2))
    G = au(arb(root_upper(B * (B - tl))) + arb(root_upper(2 * mid * (B - mid) / 3)))
    return (G, au(2 * d.sqrt() * th / 3), au(d * arb(J) / 3))

def stability_cells(thi, Bhi, dhi, taulo, tauhi):
    c0, c1, c2 = moment_coefficients(Bhi, dhi, taulo, tauhi)
    return np.minimum(Bhi, up(up(up(c2 * up(thi * thi)) + up(c1 * thi)) + c0))

def log_range_cells(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False):
    if Vlo <= 0:
        return np.ones_like(tlo)
    t2lo = sqpos(ia(tlo))[0]
    t3lo = cubepos(ia(tlo))[0]
    t3hi = cubepos(ia(thi))[1]
    k = up(KHI * t3hi)
    A = down(down(t2lo / 2) - up(up(2 * k) * au(arb(dhi).sqrt())))
    a = down(down(np.maximum(0, A) / au(arb(Vhi).sqrt())) + down(KLO * t3lo))
    c = up(k * Bhi)
    m = m_upper(thi, dhi) if deleted else np.zeros_like(thi)
    valid = (A >= 0) & (a > 0)
    aa = np.where(valid, a, 1.0)
    x = (down(down(k / 2) / aa), up(up(k / 2) / aa))
    lo = down(down(c + np.maximum(m, x[0])) / aa)
    hi = up(up(c + np.maximum(m, x[1])) / aa)
    tl = np.maximum(taulo, np.minimum(tauhi, lo))
    th = np.maximum(taulo, np.minimum(tauhi, hi))
    q = np.maximum(0, down(down(aa * tl) - c))
    penalty = np.maximum(0, down(down(q * q) - up(m * m)))
    scalar = np.where(valid, up(up(k * th) - penalty), up(k * tauhi))
    exponent = up(up(up(c + m) + scalar) - down(Vlo / 2 * t2lo))
    return exp_upper_from_nonnegative_lower(np.maximum(0, -exponent))

def feasible_clip(box, Lhi):
    dl, dh, bl, bh, tl, th = box
    for _ in range(2):
        dl = max(dl, max(0, al(arb(tl) ** 2)))
        dh = min(dh, 1.0, au(arb(th) ** (arb(2) / 3)))
        if dl > dh:
            return None
        dl3 = arb(dl) ** (arb(3) / 2)
        dh3 = arb(dh) ** (arb(3) / 2)
        bl = max(bl, max(0, al(dl3)))
        bh = min(bh, Lhi, au(arb(Lhi) - arb(tl) + dh3))
        if bl > bh:
            return None
        tl = max(tl, max(0, al(dl3)))
        th = min(th, 1.0, Lhi, au(arb(dh).sqrt()), au(arb(Lhi) - arb(bl) + dh3))
        if tl > th:
            return None
    return (dl, dh, bl, bh, tl, th)

class StableWeights(Weights):

    def __init__(self, Lhi, N=1024):
        super().__init__(Lhi, N, True)

    @staticmethod
    def prefactor(thi, Bhi, dhi, taulo, tauhi):
        return stability_cells(thi, Bhi, dhi, taulo, tauhi)

    def single_cells(self, t, dlo, dhi, bhi):
        d = arb(dhi)
        Ej = au(arb(bhi) + d * d.sqrt())
        return f_cells(t, dlo, Ej, True)

    def envelopes(self, Llo, Lhi, box):
        clipped = feasible_clip(box, Lhi)
        if clipped is None:
            return None
        dl, dh, bl, bh, tl, th = clipped
        dla, dha = (arb(dl), arb(dh))
        dl3 = dla * dla.sqrt()
        dh3 = dha * dha.sqrt()
        wl = max(0, al(1 - dha))
        wh = au(1 - dla)
        BRlo = max(0, al(arb(Llo) - arb(bh)))
        BRhi = max(0, au(arb(Lhi) - arb(bl)))
        trlo = max(0, al(arb(tl) - dh3))
        trhi = max(0, au(arb(th) - dl3))
        trhi = min(trhi, BRhi)
        d2 = min(dh, wh, au(arb(trhi) ** (arb(2) / 3)), au(arb(BRhi) ** (arb(2) / 3)))
        E = au(arb(Lhi) + arb(th))
        ER = au(arb(BRhi) + arb(trhi))
        Ej = au(arb(bh) + dh3)
        lamlo = max(0, al(arb(bl) / arb(Lhi)))
        lamhi = min(1, au(arb(bh) / arb(Llo)))

        def remainder(t, low, high):
            return np.minimum(f_cells(t, wl, ER), log_range_cells(low, high, wl, wh, BRhi, d2, trlo, trhi))
        R = remainder(self.t, self.tlo, self.thi)
        A = self.single_cells(self.t, dl, dh, bh)
        D = np.minimum(exp_cubic_cells(self.tlo, self.thi, wl, ER, d2), log_range_cells(self.tlo, self.thi, wl, wh, BRhi, d2, trlo, trhi, True))
        cfirst = self.prefactor(self.thi, bh, dh, max(0, al(dl3)), au(dh3))
        cremain = self.prefactor(self.thi, BRhi, d2, trlo, trhi)
        redfirst = np.minimum(1, up(cfirst / bl)) if bl > 0 else np.ones_like(self.thi)
        redrest = np.minimum(1, up(cremain / BRlo)) if BRlo > 0 else np.ones_like(self.thi)
        first = np.minimum(1, up(R * redfirst))
        other = np.minimum(1, up(up(A * D) * redrest))

        def mix(lam):
            return up(up(lam * first) + up(up(1 - lam) * other))
        B = np.minimum(1, np.maximum(mix(lamlo), mix(lamhi)))
        if getattr(self, 'direct_sum', False):
            direct = up(up(up(cfirst * R) + up(up(cremain * A) * D)) / Llo)
            B = np.minimum(B, direct)
        total = self.prefactor(self.thi, Lhi, dh, tl, th)
        total = np.minimum(1, up(total / Llo))
        uniform = log_range_cells(self.tlo, self.thi, 1.0, 1.0, Lhi, dh, tl, th, True)
        B = np.minimum(B, up(total * uniform))
        B = np.minimum(B, exp_cubic_cells(self.tlo, self.thi, 1.0, E, dh))
        hl, hh = (self.ht[0][:-1], self.ht[1][1:])
        Rhi = remainder(self.ht, hl, hh)
        F = np.minimum.reduce([f_cells(self.ht, 1.0, E), up(self.single_cells(self.ht, dl, dh, bh) * Rhi), log_range_cells(hl, hh, 1.0, 1.0, Lhi, dh, tl, th)])
        lowcf = np.minimum.reduce([f_cells(self.t, 1.0, E), up(A * R), log_range_cells(self.tlo, self.thi, 1.0, 1.0, Lhi, dh, tl, th)])
        return (B, F, lowcf)

    def box(self, Llo, Lhi, dl, dh, bl, bh, tl, th):
        result = self.envelopes(Llo, Lhi, (dl, dh, bl, bh, tl, th))
        if result is None:
            return None
        B, F, lowcf = result
        return self.integrate(B, F, Llo, lowcf)

def stable_checks():
    rng = np.random.default_rng(2626)
    scalar_count = 0
    log_count = 0
    for _ in range(80):
        B = float(rng.uniform(0.001, 1.7))
        tauhi = float(rng.uniform(0, B))
        taulo = float(rng.uniform(0, tauhi))
        d = float(rng.uniform(0.001, 1.0))
        ts = np.array([0.001, 0.1, 0.5, 1.0, 2.0, 4.0])
        out = stability_cells(ts, B, d, taulo, tauhi)
        for tau in [taulo, (taulo + tauhi) / 2, tauhi]:
            bb, tt, dd = (arb(B), arb(tau), arb(d))
            J = ((bb * bb + bb * (bb * bb + 8 * tt * tt).sqrt()) / 2 - 2 * tt * tt).sqrt()
            G = (bb * (bb - tt)).sqrt() + (2 * tt * (bb - tt) / 3).sqrt()
            for t, y in zip(ts, out):
                x = arb(float(t))
                ref = 2 * x * dd.sqrt() * tt / 3 + x * x * dd * J / 3 + G
                assert y == B or arb(float(y)) >= ref, (B, d, taulo, tauhi, tau, t, y, ref)
                scalar_count += 1
        V = float(rng.uniform(0.03, 1))
        th = float(rng.uniform(0.01, 5))
        tl = max(0.001, th - 0.02)
        for deleted in [False, True]:
            y = float(log_range_cells(np.array([tl]), np.array([th]), V, V, B, d, taulo, tauhi, deleted)[0])
            for t in [tl, (tl + th) / 2, th]:
                x = arb(t)
                k = K * x ** 3
                A = x * x / 2 - 2 * k * arb(d).sqrt()
                m = min(x * arb(d).sqrt(), 1 / (6 * K))
                m = m * m / 2 - 2 * K * m ** 3 if deleted else arb(0)
                if A >= 0:
                    a = A / arb(V).sqrt() + k
                    c = k * arb(B)
                    tau = max(arb(taulo), min(arb(tauhi), (c + max(m, k / (2 * a))) / a))
                    Q = max(arb(0), a * tau - c) ** 2
                    penalty = max(arb(0), Q - m * m)
                else:
                    tau = arb(tauhi)
                    penalty = arb(0)
                ref = (-arb(V) * x * x / 2 + k * arb(B) + m + k * tau - penalty).exp()
                assert y == 1 or arb(y) >= ref, (V, B, d, taulo, tauhi, t, y, ref)
                log_count += 1
    return {'stability_enclosures': scalar_count, 'log_range_enclosures': log_count, 'all_passed': True}

def finite_law_checks(weights_class=StableWeights, extra=False):
    rng = np.random.default_rng(260926)
    count = 0
    checks = 0
    for n in [2, 3, 4, 8, 20]:
        kinds = ['fair', 'biased', 'unequal', 'three atoms']
        if extra:
            kinds += ['rare symmetric', 'rare asymmetric']
        for kind in kinds:
            laws = []
            vs = []
            bs = []
            for _ in range(n):
                if kind == 'rare symmetric':
                    e = 3e-11
                    x = np.array([-1000.0, -1.0, 1.0, 1000.0])
                    p = np.array([e / 2, (1 - e) / 2, (1 - e) / 2, e / 2])
                elif kind == 'rare asymmetric':
                    e = 3e-11
                    x = np.array([-1.0, 1.0, 1000.0])
                    p = np.array([(1 - e) / 2, (1 - e) / 2, e])
                elif kind == 'three atoms':
                    x = np.array([-1.0, 0.0, 2.0])
                    p = rng.dirichlet([2, 4, 3])
                else:
                    p0 = 0.5 if kind == 'fair' else 0.418861 if kind == 'biased' else rng.uniform(0.01, 0.99)
                    p = np.array([p0, 1 - p0])
                    x = np.array([-1.0, 1.0])
                x = (x - np.dot(x, p)) * (1 if kind in ['fair', 'biased'] else rng.uniform(0.3, 2))
                laws.append((x, p))
                vs.append(np.dot(p, x * x))
                bs.append(np.dot(p, np.abs(x) ** 3))
            V = math.fsum(vs)
            vs = np.array(vs) / V
            bs = np.array(bs) / V ** 1.5
            laws = [(x / math.sqrt(V), p) for x, p in laws]
            L = float(bs.sum())
            tau = float(np.sum(vs ** 1.5))
            j = int(vs.argmax())
            d = float(vs[j])
            b = float(bs[j])
            eps = 2e-12
            Llo, Lhi = (L - eps, L + eps)
            w = weights_class(Lhi, 128)
            box = (d - eps, d + eps, b - eps, b + eps, tau - eps, tau + eps)
            result = w.envelopes(Llo, Lhi, box)
            assert result is not None, (kind, n, box)
            B, F, low = result

            def actual(t):
                fs = np.array([np.exp(1j * np.outer(t, x)) @ p for x, p in laws])
                stars = np.array([np.exp(0.5j * np.outer(t, x)) * np.sinc(np.outer(t, x) / (2 * math.pi)) @ (p * x * x / v) for (x, p), v in zip(laws, vs)])
                full = np.prod(fs, axis=0)
                star = np.zeros_like(full)
                for i in range(n):
                    star += vs[i] * stars[i] * np.prod(np.delete(fs, i, axis=0), axis=0)
                return (np.abs(full), np.abs(full - star))
            for t in [np.maximum(1e-14, w.tlo), w.thi]:
                ff, gg = actual(t)
                assert np.all(ff <= low + 1e-11), (kind, n, 'CF low')
                error = gg - L * t / 2 * B
                at = int(np.argmax(error))
                assert error[at] <= 1e-10, (kind, n, 'zero bias', t[at], gg[at], B[at], error[at])
                checks += len(t)
            for t in [w.ht[0][:-1], w.ht[1][1:]]:
                ff, _ = actual(t)
                assert np.all(ff <= F + 1e-11), (kind, n, 'CF high')
            count += 1
    return {'cases': count, 'zero_bias_points': checks, 'all_passed': True}
