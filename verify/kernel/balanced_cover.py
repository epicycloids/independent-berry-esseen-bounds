import hashlib, json, math, time
from flint import arb
from interval_bounds import al, au
import stable_cover as core

def balanced_cover_band(lo, hi, target, N, state=None, seconds_limit=190, kind='stability'):
    began = time.monotonic()
    choose = core.weight_factory(hi, N, kind, target)
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
        tree = core.unpack(state['tree_zlib_base64'])
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
                worst = value
                worstbox = box
            continue
        if depth >= 60:
            stack.append((raw, depth))
            reason = 'depth limit'
            bad = {'box': box, 'upper': value}
            break
        dl, dh, bl, bh, tl, th = box
        widths = [dh - dl, 0.25 * (bh - bl), 2 * (th - tl)]
        axis = ['D', 'B', 'T'][max(range(3), key=lambda i: widths[i])]
        left, right = core.children(box, axis)
        tree.append(axis)
        stack.extend([(right, depth + 1), (left, depth + 1)])
    encoded = core.pack(tree)
    record = dict(common, complete=not stack, method='stability', nodes=len(tree), accepted=accepted, infeasible=infeasible, upper=worst, worst_box=worstbox, tree_zlib_base64=encoded, tree_sha256=hashlib.sha256(''.join(tree).encode()).hexdigest(), leaf_chain_sha256=chain, elapsed_seconds=elapsed + time.monotonic() - began, resumptions=resumptions)
    if stack:
        record.update(stack=stack, unresolved=len(stack), reason=reason, bad=bad)
    return record
