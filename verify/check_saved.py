"""Check saved proof records without recomputing accepted-leaf quadrature."""
import base64
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import zlib

from scalars import scalar_bounds

ROOT = Path(__file__).resolve().parent.parent
COMMON = {'Llo','Lhi','N','T','complete','kind','method','nodes','root','s',
          'target','upper','weight_parameter'}
TREE = {'accepted','infeasible','leaf_chain_sha256','tree_sha256',
        'tree_zlib_base64','worst_box'}


class CertificateError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise CertificateError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def finite_number(value):
    return type(value) in (int, float) and math.isfinite(value)


def tree_counts(encoded, expected_hash):
    try:
        raw = zlib.decompress(base64.b64decode(encoded, validate=True))
    except (ValueError, zlib.error) as error:
        raise CertificateError('Invalid compressed tree') from error
    require(sha256(raw) == expected_hash, 'Tree digest mismatch')
    slots = 1
    counts = {'nodes': 0, 'accepted': 0, 'infeasible': 0}
    for code in raw:
        require(slots > 0, 'Tree has trailing nodes')
        slots -= 1
        counts['nodes'] += 1
        if code == ord('A'):
            counts['accepted'] += 1
        elif code == ord('X'):
            counts['infeasible'] += 1
        elif code in b'DBT':
            slots += 2
        else:
            raise CertificateError('Unknown tree instruction')
    require(slots == 0, 'Tree is incomplete')
    return counts


def validate_cover(cover, result):
    require(cover['schema'] == 'independent-be-cover-v1', 'Unsupported cover schema')
    require(result['schema'] == 'independent-be-result-v1', 'Unsupported result schema')
    require(result['status'] == 'proof-candidate', 'Unexpected mathematical status')
    require(result['independent_mathematical_review_completed'] is False, 'Unsupported review promotion')
    require(result['sharp_conjecture_resolved'] is False, 'Unsupported sharpness claim')
    require(cover['kind'] == 'gaussian', 'Wrong evaluator kind')
    require(cover['target'] == result['comparison_target'] == '0.475', 'Target mismatch')
    claimed = Fraction(result['upper_bound'])
    require(claimed == Fraction('0.474999998'), 'Wrong candidate constant')
    scalar_metadata = {'small_L_coefficient': '0.3413', 'small_L_source_limit': '0.01',
        'small_L_upper_decimal': '0.471750509536', 'variance_upper_decimal': '0.540936541549',
        'large_L_upper_decimal': '0.470379601347', 'kappa': '0.099162'}
    require(all(result[k] == v for k, v in scalar_metadata.items()), 'Scalar metadata mismatch')
    require(result['all_leaf_quadratures_replayed'] is False, 'Unsupported numerical replay promotion')
    domain = cover['domain']
    for side in ('lower', 'upper'):
        value = domain[side]
        require(finite_number(value) and value > 0, 'Invalid domain endpoint')
        require(float.fromhex(domain[side+'_hex']) == value, 'Dyadic endpoint mismatch')
        require(Fraction(domain[side+'_fraction']) == Fraction(value), 'Rational endpoint mismatch')
    require(domain['lower'] == .006 and domain['upper'] == 1.15, 'Unexpected covered domain')
    previous = domain['lower']
    counts = {'bands': 0, 'nodes': 0, 'accepted': 0, 'infeasible': 0}
    methods = Counter()
    worst = 0.0
    require(bool(cover['bands']), 'Empty cover')
    for band in cover['bands']:
        method = band['method']
        require(method in ('stability', 'uniform', 'variance'), 'Unknown acceptance method')
        require(set(band) == COMMON | (TREE if method == 'stability' else set()), 'Band fields do not match schema')
        require(band['complete'] is True, 'Incomplete band')
        require(band['kind'] == cover['kind'] and band['target'] == cover['target'], 'Band evaluator mismatch')
        require(band['N'] == 512, 'Unexpected quadrature size')
        for field in ('Llo','Lhi','T','s','upper','weight_parameter'):
            require(finite_number(band[field]), f'Invalid {field}')
        require(band['Llo'] == previous and band['Lhi'] > previous, 'Gap, overlap or reversed band')
        require(band['Lhi'] <= domain['upper'], 'Band beyond domain')
        require(band['T'] > 0 and 0 < band['s'] < 1, 'Invalid smoothing parameters')
        require(band['weight_parameter'] == band['Lhi'], 'Weight parameter mismatch')
        require(0 <= band['upper'] and Fraction(band['upper']) < claimed, 'Upper bound fails stated constant')
        require(len(band['root']) == 6 and all(finite_number(x) for x in band['root']), 'Invalid root box')
        if method == 'stability':
            local = tree_counts(band['tree_zlib_base64'], band['tree_sha256'])
            require(all(type(band[k]) is int and band[k] == v for k, v in local.items()), 'Tree counts disagree')
            require(len(band['leaf_chain_sha256']) == 64, 'Invalid accepted-leaf digest')
            require(len(bytes.fromhex(band['leaf_chain_sha256'])) == 32, 'Invalid accepted-leaf digest')
        else:
            require(band['nodes'] == 1, 'Invalid single-node band')
            local = {'nodes': 1, 'accepted': 1, 'infeasible': 0}
        if method == 'variance':
            require(Fraction('0.54093655')/Fraction(band['Llo']) < Fraction(cover['target']), 'Invalid variance-only band')
        for key, value in local.items():
            counts[key] += value
        counts['bands'] += 1
        methods[method] += 1
        worst = max(worst, band['upper'])
        previous = band['Lhi']
    require(previous == domain['upper'], 'Missing final band')
    require(counts == result['counts'], 'Aggregate counts disagree')
    require(dict(methods) == result['band_methods'], 'Method counts disagree')
    require(worst == result['largest_recorded_upper'], 'Maximum disagrees')
    selected = result['selected_replay_band_indices']
    require(selected == [min(range(len(cover['bands'])), key=lambda i: abs(cover['bands'][i]['Llo']-p))
                         for p in (.006,.25,.48)], 'Selected numerical bands disagree')
    return {'counts': counts, 'band_methods': dict(methods), 'largest_recorded_upper': worst}


def load_inputs(root=ROOT):
    require(__debug__, 'Run without Python optimization: the frozen kernel uses assertions')
    result = json.loads((root/'result.json').read_text())
    raw = (root/'certificates/cover.json').read_bytes()
    require(sha256(raw) == result['cover_sha256'], 'Cover file digest mismatch')
    cover = json.loads(raw)
    require(len(cover['kernel_sha256']) == 9, 'Incomplete arithmetic kernel')
    for name, expected in cover['kernel_sha256'].items():
        require(Path(name).name == name and name.endswith('.py'), 'Invalid kernel filename')
        require(sha256((root/'verify/kernel'/name).read_bytes()) == expected, f'Arithmetic source changed: {name}')
    return cover, result


def check(root=ROOT):
    cover, result = load_inputs(root)
    summary = validate_cover(cover, result)
    summary['scalar_enclosures'] = scalar_bounds(cover['domain']['lower'], cover['domain']['upper'])
    summary.update(status='saved-record checks passed', cover_sha256=result['cover_sha256'],
        tree_geometry_replayed=False, accepted_leaf_quadrature_replayed=False,
        independent_mathematical_review_completed=False)
    return summary


if __name__ == '__main__':
    try:
        print(json.dumps(check(), indent=2))
    except (CertificateError, AssertionError, ArithmeticError, KeyError, ValueError) as error:
        print(f'Certificate check failed: {error}', file=sys.stderr)
        raise SystemExit(1)
