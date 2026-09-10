import base64
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import re
import time
import zlib
if __package__:
    from .full_b_paths import HERE
else:
    from full_b_paths import HERE
from flint import arb, ctx
import interval_bounds
import stable_bounds
import stable_cover as core
KIND = 'dt_verified_prefix_continuation_v1'
SCHEMA = 'be.round27.mixed-source-prefix.v1'
MAX_NODES = 2000000
GEOMETRY_PRECISION = 100
CORE_FIELDS = ('Llo', 'Lhi', 'target', 'N', 'T', 's', 'weight_parameter', 'root', 'kind', 'split_axes', 'b_policy', 'cf_helper_sha256', 'method', 'complete', 'nodes', 'accepted', 'infeasible', 'upper', 'worst_box', 'tree_zlib_base64', 'tree_sha256', 'leaf_chain_sha256', 'stack', 'unresolved')

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def _require(condition, message):
    if not condition:
        raise ValueError(message)

def _hash(value):
    _require(isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None, 'Expected a lowercase SHA256 digest')

def _number(value):
    _require(type(value) in (int, float) and math.isfinite(value), 'Expected a finite number')
    return value

def _target(value):
    _require(isinstance(value, str), 'Target must be an exact decimal/rational string')
    answer = Q(value)
    _require(answer > 0, 'Target must be positive')
    return answer

def geometry_sources():
    return {name: hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() for name, module in (('interval_bounds.py', interval_bounds), ('stable_bounds.py', stable_bounds), ('stable_cover.py', core))}

def _sources(mapping):
    _require(isinstance(mapping, dict) and mapping, 'A complete source mapping is required')
    for name, value in mapping.items():
        _require(isinstance(name, str) and name, 'Invalid source name')
        _hash(value)
    for name, value in geometry_sources().items():
        _require(mapping.get(name) == value, 'Incompatible frozen geometry source: ' + name)

def _canonical_root(hi):
    return (0.0, min(1.0, interval_bounds.au(arb(hi) ** (arb(2) / 3))), 0.0, hi, 0.0, min(1.0, hi))

def _unpack(record):
    n = record['nodes']
    _require(type(n) is int and 0 <= n <= MAX_NODES, 'Invalid tree length')
    compressed = base64.b64decode(record['tree_zlib_base64'], validate=True)
    decoder = zlib.decompressobj()
    raw = decoder.decompress(compressed, MAX_NODES + 1)
    _require(decoder.eof and (not decoder.unused_data) and (not decoder.unconsumed_tail), 'Incomplete, trailing, or excessive compressed tree')
    _require(len(raw) == n, 'Saved node count differs from the tree')
    _require(hashlib.sha256(raw).hexdigest() == record['tree_sha256'], 'Tree checksum differs')
    tree = raw.decode('ascii')
    _require(all((c in 'ADTX' for c in tree)), 'Only closed D/tau trees are supported')
    return tree

def _stack(rows):
    answer = []
    for raw, depth in rows:
        _require(len(raw) == 6 and type(depth) is int and (0 <= depth <= 60), 'Invalid frontier entry')
        box = tuple((_number(x) for x in raw))
        _require(all((box[i] <= box[i + 1] for i in (0, 2, 4))), 'Unordered frontier box')
        answer.append((box, depth))
    return answer

def _chain(chain, box, value):
    return hashlib.sha256(bytes.fromhex(chain) + json.dumps([box, value], separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def _evaluate(callback, descriptor, lo, hi, box):
    previous = ctx.prec
    try:
        return _number(callback(deepcopy(descriptor), lo, hi, box))
    finally:
        ctx.prec = previous

def _walk(tree, frontier, hi, *, evaluate=None, descriptor=None, lo=None, target=None):
    stack = list(frontier)
    accepted = infeasible = 0
    chain = '0' * 64
    worst = 0.0
    worstbox = None
    for code in tree:
        _require(bool(stack), 'Tree extends past the closed root/frontier')
        raw, depth = stack.pop()
        box = core.feasible_clip(raw, hi)
        if code == 'X':
            _require(box is None, 'A recorded infeasible leaf is feasible')
            infeasible += 1
            continue
        _require(box is not None, 'A recorded feasible node is empty')
        if code == 'A':
            accepted += 1
            if evaluate is not None:
                value = _evaluate(evaluate, descriptor, lo, hi, box)
                _require(Q(value) < target, 'Numerical leaf does not satisfy its own target')
                chain = _chain(chain, box, value)
                if value > worst:
                    worst, worstbox = (value, box)
            continue
        _require(depth < 60, 'Split exceeds the inherited depth limit')
        offset = 0 if code == 'D' else 4
        l, r = box[offset:offset + 2]
        _require(l < (l + r) / 2 < r, 'Unrepresentable closed split')
        left, right = core.children(box, code)
        stack.extend(((right, depth + 1), (left, depth + 1)))
    return dict(stack=stack, accepted=accepted, infeasible=infeasible, chain=chain, upper=worst, worst_box=worstbox)

def _descriptor(record, sources, parameters):
    return {**{k: deepcopy(record[k]) for k in ('kind', 'target', 'N', 'T', 's', 'weight_parameter')}, 'source_sha256': deepcopy(sources), 'parameters': deepcopy(parameters)}

def _origin(record, provenance, lo, hi, target, root):
    digest(record)
    digest(provenance)
    for key in ('record_sha256', 'manifest_sha256', 'checkpoint_sha256'):
        _hash(provenance[key])
    _sources(provenance['source_sha256'])
    _require(provenance.get('geometry_precision_bits') == GEOMETRY_PRECISION, 'Origin must attest the frozen 100-bit geometry convention')
    _require(digest(record) == provenance['record_sha256'], 'Original record hash mismatch')
    _require(isinstance(provenance.get('evaluator_parameters'), dict), 'Origin evaluator parameters must be explicitly recorded')
    _require(record['Llo'] == lo and record['Lhi'] == hi, 'Reuse cannot change the L interval')
    _require(list(record['root']) == list(root) == list(_canonical_root(hi)), 'Reuse cannot reduce or change the canonical root')
    _require(record['kind'] != KIND and record['method'] == 'stability' and (record['complete'] is False) and (record['split_axes'] == 'DT'), 'Require an original unfinished D/tau stability record')
    old_target = _target(record['target'])
    _evaluator(_descriptor(record, provenance['source_sha256'], provenance['evaluator_parameters']), record['target'], hi)
    upper = _number(record['upper'])
    _require(0 <= upper and Q(upper) < old_target and (Q(upper) < target), 'The old accepted upper must satisfy both old and requested targets')
    _hash(record['leaf_chain_sha256'])
    tree = _unpack(record)
    replay = _walk(tree, [(tuple(root), 0)], hi)
    saved = _stack(record['stack'])
    _require(saved and saved == replay['stack'] and (record['unresolved'] == len(saved)), 'Origin frontier differs from exact root replay')
    for key in ('accepted', 'infeasible'):
        _require(type(record[key]) is int and record[key] == replay[key], 'Origin ' + key + ' mismatch')
    if not record['accepted']:
        _require(upper == 0 and record['worst_box'] is None and (record['leaf_chain_sha256'] == '0' * 64), 'Empty accepted prefix has nonempty numerical evidence')
    compact = {k: deepcopy(record[k]) for k in CORE_FIELDS if k in record}
    return (dict(record_sha256=provenance['record_sha256'], core_sha256=digest(compact), record=compact, provenance=deepcopy(provenance)), tree, replay['stack'])

def _evaluator(evaluator, target, hi):
    digest(evaluator)
    _require(isinstance(evaluator.get('kind'), str) and evaluator['kind'] and (evaluator['kind'] != KIND), 'A distinct numerical evaluator kind is required')
    _require(evaluator.get('target') == target, 'New evaluator target mismatch')
    _require(type(evaluator.get('N')) is int and evaluator['N'] > 0, 'Invalid evaluator mesh')
    _require(_number(evaluator['T']) > 0 and 0 < _number(evaluator['s']) < 1, 'Invalid smoothing parameters')
    _require(evaluator['weight_parameter'] == hi, 'Evaluator weight parameter must retain Lhi')
    _require(isinstance(evaluator.get('parameters'), dict), 'Explicit evaluator parameters required')
    _sources(evaluator['source_sha256'])

def _empty_suffix():
    return dict(nodes=0, accepted=0, infeasible=0, upper=0.0, worst_box=None, tree_zlib_base64=core.pack(''), tree_sha256=hashlib.sha256(b'').hexdigest(), leaf_chain_sha256='0' * 64)

def _priority(evaluator):
    value = evaluator.get('dt_priority', (1.0, 1.0))
    _require(isinstance(value, (list, tuple)) and len(value) == 2, 'D/tau priority must contain exactly two positive finite weights')
    weights = tuple((_number(item) for item in value))
    _require(all((item > 0 for item in weights)), 'D/tau priority weights must be positive')
    return weights

def prepare(origin_record, provenance, evaluator, *, Llo, Lhi, target, root):
    previous = ctx.prec
    try:
        ctx.prec = GEOMETRY_PRECISION
        _require(0 < _number(Llo) < _number(Lhi), 'Require an unchanged positive L band')
        limit = _target(target)
        _evaluator(evaluator, target, Lhi)
        origin, _, frontier = _origin(origin_record, provenance, Llo, Lhi, limit, root)
        return dict(schema=SCHEMA, kind=KIND, method='mixed_source_prefix', implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), Llo=Llo, Lhi=Lhi, target=target, root=list(root), geometry_precision_bits=GEOMETRY_PRECISION, dt_priority=list(_priority(evaluator)), origin=origin, evaluator=deepcopy(evaluator), continuation=_empty_suffix(), complete=False, upper=origin_record['upper'], nodes=origin_record['nodes'], accepted=origin_record['accepted'], infeasible=origin_record['infeasible'], stack=deepcopy(frontier), unresolved=len(frontier), reason='prepared frontier', elapsed_seconds=0.0, resumptions=0, leaf_chain_semantics='Independent zero-seeded origin and continuation chains; never a new-evaluator chain for the prefix')
    finally:
        ctx.prec = previous

def verify(record, origin_record, provenance, evaluator, *, numerical=False, evaluate_origin=None, evaluate_new=None):
    previous = ctx.prec
    try:
        ctx.prec = GEOMETRY_PRECISION
        digest(record)
        _require(record['schema'] == SCHEMA and record['kind'] == KIND and (record['method'] == 'mixed_source_prefix'), 'Not a mixed-source prefix record')
        _require(record['geometry_precision_bits'] == GEOMETRY_PRECISION, 'Geometry precision changed')
        _require(record['implementation_sha256'] == hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'Prefix wrapper implementation changed')
        _require(type(record['resumptions']) is int and record['resumptions'] >= 0 and (_number(record['elapsed_seconds']) >= 0), 'Invalid continuation accounting')
        lo, hi = (record['Llo'], record['Lhi'])
        _require(0 < _number(lo) < _number(hi), 'Invalid L band')
        limit = _target(record['target'])
        _evaluator(evaluator, record['target'], hi)
        _require(digest(record['evaluator']) == digest(evaluator), 'Continuation evaluator/configuration changed')
        _require(record['dt_priority'] == list(_priority(evaluator)), 'Continuation D/tau priority changed')
        origin, tree, frontier = _origin(origin_record, provenance, lo, hi, limit, record['root'])
        _require(digest(record['origin']) == digest(origin), 'Origin witness differs from its pinned record')
        suffix = record['continuation']
        newtree = _unpack(suffix)
        _require(len(tree) + len(newtree) <= MAX_NODES, 'Combined tree exceeds the node limit')
        if numerical:
            _require(callable(evaluate_origin) and callable(evaluate_new), 'Mixed numerical replay requires both source-specific callbacks')
            olddesc = _descriptor(origin_record, provenance['source_sha256'], provenance['evaluator_parameters'])
            oldrun = _walk(tree, [(tuple(record['root']), 0)], hi, evaluate=evaluate_origin, descriptor=olddesc, lo=lo, target=_target(origin_record['target']))
            _require(oldrun['chain'] == origin_record['leaf_chain_sha256'] and oldrun['upper'] == origin_record['upper'], 'Original numerical leaf chain/upper changed')
            _require(digest(oldrun['worst_box']) == digest(origin_record['worst_box']), 'Original numerical worst box changed')
        run = _walk(newtree, frontier, hi, evaluate=evaluate_new if numerical else None, descriptor=evaluator, lo=lo, target=limit)
        for key in ('accepted', 'infeasible'):
            _require(type(suffix[key]) is int and suffix[key] == run[key], 'Continuation ' + key + ' mismatch')
            _require(record[key] == origin_record[key] + suffix[key], 'Combined ' + key + ' mismatch')
        upper = _number(suffix['upper'])
        _hash(suffix['leaf_chain_sha256'])
        _require(0 <= upper and Q(upper) < limit, 'Continuation accepted upper exceeds target')
        if not suffix['accepted']:
            _require(upper == 0 and suffix['worst_box'] is None and (suffix['leaf_chain_sha256'] == '0' * 64), 'Empty continuation has nonempty numerical evidence')
        if numerical:
            _require(run['chain'] == suffix['leaf_chain_sha256'] and run['upper'] == upper, 'Continuation numerical leaf chain/upper changed')
            _require(digest(run['worst_box']) == digest(suffix['worst_box']), 'Continuation numerical worst box changed')
        _require(record['upper'] == max(origin_record['upper'], upper), 'Combined upper mismatch')
        _require(record['nodes'] == origin_record['nodes'] + suffix['nodes'], 'Combined node count mismatch')
        _require(_stack(record['stack']) == run['stack'] and record['unresolved'] == len(run['stack']), 'Continuation frontier differs from exact replay')
        _require(type(record['complete']) is bool and record['complete'] == (not run['stack']), 'Completion flag disagrees with exact closed topology')
        return dict(valid_topology=True, numerically_replayed=numerical, complete=record['complete'], origin_nodes=len(tree), continuation_nodes=len(newtree), unresolved=len(run['stack']))
    finally:
        ctx.prec = previous

def advance(record, origin_record, provenance, evaluator, evaluate, *, seconds_limit=180, max_new_nodes=200000):
    _require(callable(evaluate), 'An explicit new evaluator callback is required')
    _require(type(max_new_nodes) is int and 0 <= max_new_nodes <= MAX_NODES, 'Invalid node allowance')
    _require(0 <= _number(seconds_limit) <= 190, 'Require a work limit of at most 190 seconds')
    verify(record, origin_record, provenance, evaluator)
    out = deepcopy(record)
    if out['complete']:
        return out
    began = time.monotonic()
    previous = ctx.prec
    try:
        ctx.prec = GEOMETRY_PRECISION
        stack = _stack(out['stack'])
        suffix = out['continuation']
        tree = list(_unpack(suffix))
        start = len(tree)
        lo, hi = (out['Llo'], out['Lhi'])
        target = _target(out['target'])
        reason = 'bounded work interval'
        while stack:
            if out['origin']['record']['nodes'] + len(tree) >= MAX_NODES:
                reason = 'node limit'
                break
            if len(tree) - start >= max_new_nodes or time.monotonic() - began >= seconds_limit:
                break
            raw, depth = stack.pop()
            box = core.feasible_clip(raw, hi)
            if box is None:
                tree.append('X')
                suffix['infeasible'] += 1
                continue
            value = _evaluate(evaluate, evaluator, lo, hi, box)
            if Q(value) < target:
                tree.append('A')
                suffix['accepted'] += 1
                suffix['leaf_chain_sha256'] = _chain(suffix['leaf_chain_sha256'], box, value)
                if value > suffix['upper']:
                    suffix['upper'], suffix['worst_box'] = (value, box)
                continue
            wd, wt = out['dt_priority']
            choices = [(weight * (box[i + 1] - box[i]), axis) for i, axis, weight in ((0, 'D', wd), (4, 'T', wt)) if box[i] < (box[i] + box[i + 1]) / 2 < box[i + 1]]
            if depth >= 60 or not choices:
                stack.append((raw, depth))
                reason = 'depth limit' if depth >= 60 else 'no representable D/tau split'
                break
            axis = max(choices, key=lambda row: row[0])[1]
            left, right = core.children(box, axis)
            tree.append(axis)
            stack.extend(((right, depth + 1), (left, depth + 1)))
        raw = ''.join(tree).encode()
        suffix.update(nodes=len(tree), tree_zlib_base64=core.pack(tree), tree_sha256=hashlib.sha256(raw).hexdigest())
        out.update(complete=not stack, stack=stack, unresolved=len(stack), reason=None if not stack else reason, upper=max(origin_record['upper'], suffix['upper']), nodes=origin_record['nodes'] + suffix['nodes'], accepted=origin_record['accepted'] + suffix['accepted'], infeasible=origin_record['infeasible'] + suffix['infeasible'], elapsed_seconds=out['elapsed_seconds'] + time.monotonic() - began, resumptions=out['resumptions'] + 1)
        verify(out, origin_record, provenance, evaluator)
        return out
    finally:
        ctx.prec = previous
