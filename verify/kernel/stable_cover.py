import base64
import hashlib
import json
import math
import time
import zlib
from flint import arb
from interval_bounds import al, au, elementary_checks
from stable_bounds import StableWeights, feasible_clip, stable_checks
from directional_bounds import DirectionalWeights, directional_checks
from refined_integral import RefinedWeights

def pack(tree):
    return base64.b64encode(zlib.compress(''.join(tree).encode(), 9)).decode()

def unpack(code):
    return list(zlib.decompress(base64.b64decode(code)).decode())

def children(box, axis):
    offset = {'D': 0, 'B': 2, 'T': 4}[axis]
    lo, hi = box[offset:offset + 2]
    mid = (lo + hi) / 2
    assert lo < mid < hi, (box, axis)
    left = list(box)
    right = list(box)
    left[offset + 1] = mid
    right[offset] = mid
    return (tuple(left), tuple(right))

def weight_factory(hi, N, kind, target=None):
    cache = {}

    def choose(tauhi):
        bucket = min(20, max(1, math.ceil(20 * tauhi / hi))) if kind == 'refined' else min(20, max(1, round(20 * tauhi / hi))) if kind == 'cosine' else 20
        if bucket not in cache:
            if kind == 'cosine':
                from paired_bounds import CosineWeights

                class Portfolio:

                    def __init__(self):
                        self.weights = [CosineWeights(hi, N, hi * bucket / 20, scale, shift) for scale, shift in [(1.0, 0.0), (1.04, 0.0), (1.0, -0.02), (1.04, -0.02)]]
                        self.T = self.weights[0].T
                        self.s = self.weights[0].s

                    def uniform(self, lo, high):
                        return self.weights[0].uniform(lo, high)

                    def box(self, *args):
                        best = float('inf')
                        for w in self.weights:
                            value = w.box(*args)
                            if value is None:
                                return None
                            best = min(best, value)
                            if target is not None and best < al(arb(target)):
                                break
                        return best
                cache[bucket] = Portfolio()
            elif kind == 'refined':
                cache[bucket] = RefinedWeights(hi, N, hi * bucket / 20)
            elif kind == 'directional':
                cache[bucket] = DirectionalWeights(hi, N)
            else:
                cache[bucket] = StableWeights(hi, N)
        return cache[bucket]
    return choose

def cover_band(lo, hi, target, N, state=None, seconds_limit=190, kind='stability'):
    began = time.monotonic()
    choose = weight_factory(hi, N, kind, target)
    w = choose(hi)
    threshold = al(arb(target))
    root = (0.0, min(1.0, au(arb(hi) ** (arb(2) / 3))), 0.0, hi, 0.0, min(1.0, hi))
    common = {'Llo': lo, 'Lhi': hi, 'target': target, 'N': N, 'T': w.T, 's': w.s, 'weight_parameter': hi, 'root': root, 'kind': kind}
    if state is None:
        vbound = au(arb('0.54093655') / arb(lo))
        if vbound < threshold:
            return dict(common, complete=True, method='variance', upper=vbound, nodes=1)
        uniform = w.uniform(lo, hi)
        if uniform < threshold:
            return dict(common, complete=True, method='uniform', upper=uniform, nodes=1)
        stack = [(root, 0)]
        tree = []
        accepted = 0
        infeasible = 0
        worst = 0.0
        worstbox = None
        chain = '0' * 64
        elapsed = 0.0
        resumptions = 0
    else:
        assert all((state[k] == common[k] for k in ['Llo', 'Lhi', 'target', 'N', 'T', 's']))
        assert state.get('kind', 'stability') == kind
        stack = state['stack']
        tree = unpack(state['tree_zlib_base64'])
        accepted = state['accepted']
        infeasible = state['infeasible']
        worst = state['upper']
        worstbox = state['worst_box']
        chain = state['leaf_chain_sha256']
        elapsed = state['elapsed_seconds']
        resumptions = state['resumptions'] + 1
    reason = None
    bad = None
    while stack:
        if time.monotonic() - began > seconds_limit:
            reason = 'bounded work interval'
            break
        if len(tree) >= 2000000:
            reason = 'node limit'
            break
        raw, depth = stack.pop()
        box = feasible_clip(raw, hi)
        if box is None:
            tree.append('X')
            infeasible += 1
            continue
        value = choose(box[-1]).box(lo, hi, *box)
        assert value is not None and math.isfinite(value), (box, value)
        if value < threshold:
            tree.append('A')
            accepted += 1
            chain = hashlib.sha256(bytes.fromhex(chain) + json.dumps([box, value], separators=(',', ':')).encode()).hexdigest()
            if value > worst:
                worst = value
                worstbox = box
            continue
        if depth >= 60:
            stack.append((raw, depth))
            reason = 'depth limit'
            bad = {'box': box, 'upper': value}
            break
        dl, dh, bl, bh, tl, th = box
        widths = [dh - dl, 2 * (bh - bl), 2 * (th - tl)]
        axis = ['D', 'B', 'T'][max(range(3), key=lambda i: widths[i])]
        left, right = children(box, axis)
        tree.append(axis)
        stack.extend([(right, depth + 1), (left, depth + 1)])
    encoded = pack(tree)
    record = dict(common, complete=not stack, method='stability', nodes=len(tree), accepted=accepted, infeasible=infeasible, upper=worst, worst_box=worstbox, tree_zlib_base64=encoded, tree_sha256=hashlib.sha256(''.join(tree).encode()).hexdigest(), leaf_chain_sha256=chain, elapsed_seconds=elapsed + time.monotonic() - began, resumptions=resumptions)
    if stack:
        record.update(stack=stack, unresolved=len(stack), reason=reason, bad=bad)
    return record

def replay_band(record, numerical=True):
    lo, hi = (record['Llo'], record['Lhi'])
    limit = al(arb(record['target']))
    choose = weight_factory(hi, record['N'], record.get('kind', 'stability'), record['target']) if numerical else None
    w = choose(hi) if numerical else None
    if numerical:
        assert w.T == record['T'] and w.s == record['s']
    assert record['complete']
    if record['method'] == 'uniform':
        if numerical:
            assert w.uniform(lo, hi) < limit
        return {'nodes': 1, 'accepted': 1, 'infeasible': 0}
    if record['method'] == 'variance':
        assert arb('0.54093655') / arb(lo) < arb(record['target'])
        return {'nodes': 1, 'accepted': 1, 'infeasible': 0}
    root = (0.0, min(1.0, au(arb(hi) ** (arb(2) / 3))), 0.0, hi, 0.0, min(1.0, hi))
    assert list(root) == list(record['root'])
    tree = unpack(record['tree_zlib_base64'])
    raw = ''.join(tree).encode()
    assert hashlib.sha256(raw).hexdigest() == record['tree_sha256']
    stream = iter(tree)
    stack = [root]
    counts = {'nodes': 0, 'accepted': 0, 'infeasible': 0}
    chain = '0' * 64
    worst = 0.0
    while stack:
        code = next(stream)
        counts['nodes'] += 1
        box = feasible_clip(stack.pop(), hi)
        if code == 'X':
            assert box is None
            counts['infeasible'] += 1
            continue
        assert box is not None
        if code == 'A':
            counts['accepted'] += 1
            if numerical:
                value = choose(box[-1]).box(lo, hi, *box)
                assert value is not None and value < limit
                chain = hashlib.sha256(bytes.fromhex(chain) + json.dumps([box, value], separators=(',', ':')).encode()).hexdigest()
                worst = max(worst, value)
            continue
        assert code in 'DBT'
        left, right = children(box, code)
        stack.extend([right, left])
    assert next(stream, None) is None
    assert all((record[k] == counts[k] for k in counts))
    if numerical:
        assert chain == record['leaf_chain_sha256'] and worst == record['upper']
    return counts
