from fractions import Fraction as F
from math import factorial, comb
import argparse
import json
import time

def add(p, ij, c):
    p[ij] = p.get(ij, F(0)) + c
    if not p[ij]:
        del p[ij]

def q_coefficients(m, n):
    k = m + n
    num = [F(0)] * (2 * k + 1)
    for i in range(m):
        for j in range(n):
            for d, c in [(2, 1), (3, 2), (4, 1)]:
                num[2 * i + 2 * j + d] += c
    v = F(4 * m * n, 9)
    num[0] -= v * F(5 - 3 * k, 2)
    num[1] -= v * (5 - 3 * k)
    num[2] -= v * F(3 * k + 1, 2)
    num[3] -= v * (3 * k + 1)
    q = [F(0)] * (2 * k - 1)
    for j in range(2 * k - 2, -1, -1):
        q[j] = num[j + 2]
        num[j + 2] -= q[j]
        num[j + 1] += 2 * q[j]
        num[j] -= q[j]
    assert not any(num)
    assert all((c > 0 for c in q))
    return q

def components(N, absolute=False):
    F1 = {}
    F0 = {}
    G = {}
    if not absolute:
        for i, c in enumerate([F(1), F(2, 3), F(1, 9)]):
            add(F1, (i, 0), c / 144)
            add(G, (i, 0), c / 36)
    for k in range(3, N + 1):
        sign = 1 if absolute else (-1) ** (k + 1)
        for m in range(1, k):
            n = k - m
            w = [F(1, factorial(2 * m) * factorial(2 * n)), F(1, factorial(2 * m) * factorial(2 * n + 1)) + F(1, factorial(2 * m + 1) * factorial(2 * n)), F(1, factorial(2 * m + 1) * factorial(2 * n + 1))]
            for j, q in enumerate(q_coefficients(m, n)):
                xp = 2 * k - 2 - j
                for i, wi in enumerate(w):
                    c = sign * q * wi
                    if xp % 2:
                        add(G, ((xp - 1) // 2 + i, (j - 1) // 2), c)
                    else:
                        a = xp // 2 + i
                        b = j // 2
                        if a == 0:
                            add(F0, (0, b - 2), c)
                        else:
                            add(F1, (a - 1, b), c)
    return (F1, F0, G)

def evaluate(p, u, v):
    return sum((c * u ** a * v ** b for (a, b), c in p.items()))

def restrict(p, u0, u1, v0, v1):
    out = {}
    for (a, b), c in p.items():
        for i in range(a + 1):
            cu = c * comb(a, i) * u0 ** (a - i) * (u1 - u0) ** i
            if not cu:
                continue
            for j in range(b + 1):
                add(out, (i, j), cu * comb(b, j) * v0 ** (b - j) * (v1 - v0) ** j)
    return out

def bernstein_lower(p, box):
    p = restrict(p, *box)
    du = max((a for a, b in p))
    dv = max((b for a, b in p))
    lo = None
    for i in range(du + 1):
        for j in range(dv + 1):
            z = sum((c * F(comb(i, a), comb(du, a)) * F(comb(j, b), comb(dv, b)) for (a, b), c in p.items() if a <= i and b <= j))
            lo = z if lo is None else min(lo, z)
    return lo

def certificate(order=14, grid=1):
    assert 3 <= order <= 24
    began = time.monotonic()
    ps = components(order)
    result = {}
    abs_lo = components(order, True)
    abs_hi = components(24, True)
    B25 = F(49, 900) * 25 ** 2 * 100 ** 25 / F(factorial(48))
    assert F(100 * 26 ** 2, 25 ** 2 * 50 * 49) < F(1, 20)
    infinite_component_tail = F(5, 19) * B25
    assert infinite_component_tail < F(1, 10 ** 9)
    for name, p in zip(['F1', 'F0', 'G'], ps):
        mins = []
        for i in range(grid):
            for j in range(grid):
                box = (F(4 * i, grid), F(4 * (i + 1), grid), F(25 * j, grid), F(25 * (j + 1), grid))
                mins.append(bernstein_lower(p, box))
        idx = ['F1', 'F0', 'G'].index(name)
        finite_tail = evaluate(abs_hi[idx], F(4), F(25)) - evaluate(abs_lo[idx], F(4), F(25))
        tail = finite_tail + infinite_component_tail
        assert tail < F(1, 8000), (name, tail)
        assert min(mins) - tail > F(1, 200), (name, min(mins), tail)
        result[name] = {'minimum_Bernstein_coefficient': str(min(mins)), 'minimum_float': float(min(mins)), 'absolute_tail_upper': str(tail), 'absolute_tail_float': float(tail), 'certified_component_lower': '1/200'}
    return {'status': 'exact rational scalar certificate passed', 'order': order, 'grid': grid, 'rectangle': {'x_squared': [0, 4], 'z_squared': [0, 25]}, 'components': result, 'infinite_component_tail_float': float(infinite_component_tail), 'elapsed_seconds': time.monotonic() - began}
