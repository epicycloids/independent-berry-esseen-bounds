import hashlib
import json
import math
import time
if __package__:
    from .full_b_paths import HERE
else:
    from full_b_paths import HERE
from full_b_consumer import check_helper_source
from flint import arb
from interval_bounds import al, au
import stable_cover as core
KIND = 'full_b_fast_maximum_fixed085_dual5_cutoff_v1'
SPLIT_AXES = 'DT'
B_POLICY = 'full_feasible_each_evaluation'
KINDS = {KIND: (1,), KIND + '_b2': (1, 2), KIND + '_b4': (1, 2, 4), KIND + '_b8': (1, 2, 4, 8)}
KINDS.update({kind + '_jointcf': parts for kind, parts in tuple(KINDS.items())})
COUPLING = {kind: kind.endswith('_jointcf') for kind in KINDS}

def split_axis(box):
    dl, dh, _, _, tl, th = box
    available = []
    for axis, lo, hi, weight in (('D', dl, dh, 1.0), ('T', tl, th, 2.0)):
        if lo < (lo + hi) / 2 < hi:
            available.append((weight * (hi - lo), axis))
    if not available:
        return None
    return max(available, key=lambda item: item[0])[1]

def _resume_matches(state, common):
    for key in ('Llo', 'Lhi', 'target', 'N', 'T', 's', 'weight_parameter', 'kind', 'split_axes', 'b_policy', 'cf_helper_sha256'):
        if state.get(key) != common[key]:
            raise ValueError('Resume changes ' + key)
    if list(state.get('root', ())) != list(common['root']):
        raise ValueError('Resume changes the exact root')
    tree = core.unpack(state['tree_zlib_base64'])
    if any((code not in 'ADTX' for code in tree)):
        raise ValueError('This source kind never admits a B split')
    if hashlib.sha256(''.join(tree).encode()).hexdigest() != state['tree_sha256']:
        raise ValueError('Saved tree checksum differs')
    if len(tree) != state['nodes']:
        raise ValueError('Saved node count differs')
    if state.get('complete') or not state.get('stack') or state['unresolved'] != len(state['stack']):
        raise ValueError('Resume requires its nonempty unfinished stack')
    if state.get('method') != 'stability':
        raise ValueError('Only an unfinished stability tree can resume')
    frontier = [(tuple(common['root']), 0)]
    counts = dict(accepted=0, infeasible=0)
    for code in tree:
        if not frontier:
            raise ValueError('Saved tree extends beyond its closed root')
        raw, depth = frontier.pop()
        box = core.feasible_clip(raw, common['Lhi'])
        if code == 'X':
            if box is not None:
                raise ValueError('Saved infeasible leaf is feasible')
            counts['infeasible'] += 1
        elif box is None:
            raise ValueError('Saved feasible node is empty')
        elif code == 'A':
            counts['accepted'] += 1
        else:
            if depth >= 60 or code != split_axis(box):
                raise ValueError('Saved split differs from the D/tau policy')
            left, right = core.children(box, code)
            frontier.extend(((right, depth + 1), (left, depth + 1)))
    saved = [(tuple(raw), depth) for raw, depth in state['stack']]
    if saved != frontier or any((state[k] != v for k, v in counts.items())):
        raise ValueError('Saved stack or leaf counts differ from closed-tree replay')
    return tree

def full_b_cover_band(lo, hi, target, N, state=None, seconds_limit=190, kind=KIND):
    if kind not in KINDS:
        raise ValueError('The full-b split policy requires its distinct source kind')
    began = time.monotonic()
    choose = core.weight_factory(hi, N, kind, target)
    w = choose(hi)
    threshold = al(arb(target))
    root = (0.0, min(1.0, au(arb(hi) ** (arb(2) / 3))), 0.0, hi, 0.0, min(1.0, hi))
    common = dict(Llo=lo, Lhi=hi, target=target, N=N, T=w.T, s=w.s, weight_parameter=hi, root=root, kind=kind, split_axes=SPLIT_AXES, b_policy=B_POLICY, cf_helper_sha256=check_helper_source())
    if state is None:
        vbound = au(arb('0.54093655') / arb(lo))
        if vbound < threshold:
            return dict(common, complete=True, method='variance', upper=vbound, nodes=1)
        uniform = w.uniform(lo, hi)
        if uniform < threshold:
            return dict(common, complete=True, method='uniform', upper=uniform, nodes=1)
        stack = [(root, 0)]
        tree = []
        accepted = infeasible = 0
        worst = 0.0
        worstbox = None
        chain = '0' * 64
        elapsed = 0.0
        resumptions = 0
    else:
        tree = _resume_matches(state, common)
        stack = list(state['stack'])
        accepted, infeasible, worst = (state['accepted'], state['infeasible'], state['upper'])
        worstbox, chain = (state['worst_box'], state['leaf_chain_sha256'])
        elapsed, resumptions = (state['elapsed_seconds'], state['resumptions'] + 1)
    reason = bad = None
    while stack:
        if time.monotonic() - began > seconds_limit:
            reason = 'bounded work interval'
            break
        if len(tree) >= 2000000:
            reason = 'node limit'
            break
        raw, depth = stack.pop()
        box = core.feasible_clip(raw, hi)
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
                worst, worstbox = (value, box)
            continue
        axis = split_axis(box)
        if depth >= 60 or axis is None:
            stack.append((raw, depth))
            reason = 'depth limit' if depth >= 60 else 'no representable D/tau split'
            bad = dict(box=box, upper=value)
            break
        left, right = core.children(box, axis)
        tree.append(axis)
        stack.extend([(right, depth + 1), (left, depth + 1)])
    record = dict(common, complete=not stack, method='stability', nodes=len(tree), accepted=accepted, infeasible=infeasible, upper=worst, worst_box=worstbox, tree_zlib_base64=core.pack(tree), tree_sha256=hashlib.sha256(''.join(tree).encode()).hexdigest(), leaf_chain_sha256=chain, elapsed_seconds=elapsed + time.monotonic() - began, resumptions=resumptions)
    if stack:
        record.update(stack=stack, unresolved=len(stack), reason=reason, bad=bad)
    return record
