"""Check the supplied Bernoulli certificates and selected arithmetic formulas."""
import argparse
import base64
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import zlib

from flint import ctx

import benchmarks
import common_gap
import two_arbitrary_gaps as arbitrary
import two_equal_gaps as equal

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'certificates/structure'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(name):
    return json.loads((DATA / (name+'.json')).read_text())


def source_identity(record, filename):
    actual = hashlib.sha256(Path(__file__).with_name(filename).read_bytes()).hexdigest()
    require(record['source_sha256'] == actual, 'Arithmetic source changed: '+filename)


def unpack(record):
    artifact = record['lossless_certificate']
    raw = zlib.decompress(base64.b64decode(artifact['data_b64'], validate=True))
    require(len(raw) == artifact['json_bytes'] and
            hashlib.sha256(raw).hexdigest() == artifact['sha256_json'],
            'Lossless partition differs from its certificate')
    return json.loads(raw)


def common_checks(arithmetic=False):
    result, payload = load('common_gap'), load('common_gap-input')
    source_identity(result, 'common_gap.py')
    require(payload['price'] == '199/500' and result['all_covers_certified'],
            'Wrong price or incomplete common-gap certificate')
    require(result['parameters'] == {k:v for k,v in payload.items() if k != 'source_bundle'},
            'Common-gap parameters differ')
    ctx.prec = payload['arb_precision']
    modules = {k:common_gap.load_original(k,v,payload['price'])[0]
               for k,v in payload['source_bundle'].items()}
    checks = common_gap.fixed_controls(modules['singleton'], modules['two_summand'], payload['price'])
    require(checks['passed'] and result['controls']['passed'], 'Common-gap controls failed')
    for key, nodes, leaves in [('singleton',550563,275283), ('two_summand',13641,6822)]:
        row = result[key]
        require(row['certified'] and row['pending'] == 0 and row['processed'] == nodes
                and row['accepted'] == leaves, 'Common-gap counts or completion differ')
    if arithmetic:
        fresh = common_gap.run(payload)
        require(fresh['all_covers_certified'], 'Common-gap arithmetic did not finish')
        for key in ['singleton','two_summand']:
            require(fresh[key]['processed'] == result[key]['processed'] and
                    fresh[key]['accepted'] == result[key]['accepted'],
                    'Deterministic common-gap subdivision changed')
    return {'recorded_complete':True, 'arithmetic_recomputed':arithmetic}


def equal_checks(arithmetic=False):
    result, payload = load('two_equal_gaps'), load('two_equal_gaps-input')
    source_identity(result, 'two_equal_gaps.py')
    ctx.prec = payload.get('arb_precision',192)
    require(result['all_covers_certified'], 'Incomplete two-coordinate certificate')
    checks = equal.controls()
    require(all(v for k,v in checks.items() if isinstance(v,bool)), 'Two-coordinate controls failed')
    for key in ['diagonal','full']:
        require(result[key]['certified'] and result[key]['pending'] == 0,
                'Incomplete equal-gap domain')
    if arithmetic:
        fresh = equal.run(payload)
        require(fresh['all_covers_certified'], 'Equal-gap arithmetic did not finish')
        for key in ['diagonal','full']:
            require(fresh[key]['processed'] == result[key]['processed'] and
                    fresh[key]['accepted'] == result[key]['accepted'],
                    'Deterministic equal-gap subdivision changed')
    return {'recorded_complete':True, 'arithmetic_recomputed':arithmetic}


def arbitrary_checks(arithmetic=False):
    record, payload = load('two_arbitrary_gaps-replay'), load('two_arbitrary_gaps-input')
    source_identity(record, 'two_arbitrary_gaps.py')
    ctx.prec = record['precision_bits']
    source, pending = arbitrary.load_compound(payload)
    require(len(source) == 37287 and not pending, 'Incomplete three-cube source partition')
    replay = unpack(record)
    require(not replay['unresolved'] and replay['unstarted_source_index_range'] == [37287,37287]
            and record['complete_independent_arithmetic_cover']
            and record['verified_source_leaf_count'] == len(source), 'Incomplete recorded arithmetic replay')
    require(replay['reason_codes'] == arbitrary.REASONS and
            replay['base_certificate_sha256'] == arbitrary.BASE_CERTIFICATE_PIN and
            replay['continuation_certificate_sha256'] == arbitrary.CONTINUATION_CERTIFICATE_PIN,
            'Replay mathematical inputs differ')
    groups = defaultdict(list)
    trace = hashlib.sha256()
    for index,path,code in replay['descendants']:
        require(type(index) is int and 0 <= index < len(source) and
                type(code) is int and 0 <= code < len(arbitrary.REASONS), 'Invalid replay leaf')
        arbitrary.validate_path(source[index]['owner'],path)
        require(arbitrary.REASONS[code] == source[index]['reason'],
                'Replay changed an original acceptance condition')
        groups[index].append(path)
        trace.update(arbitrary.canonical([index,path,code])+b'\n')
        if arithmetic:
            outer = arbitrary.decode(source[index]['path'])
            inner = arbitrary.decode(path,outer)
            valid, witness = arbitrary.recheck(inner,source[index]['owner'],arbitrary.REASONS[code])
            require(valid,
                    'Independent arithmetic failed on a recorded descendant')
    require(trace.hexdigest() == record['independent_trace_sha256'], 'Replay trace differs')
    require(set(groups) == set(range(len(source))) and
            len(replay['descendants']) == record['descendant_leaf_count'], 'Missing replay source or descendant')
    for paths in groups.values():
        arbitrary.forest(paths)
    controls = arbitrary.fixed_controls()
    require(controls['all_passed'], 'Independent formula or rejection control failed')
    return {'source_leaves':len(source), 'descendant_leaves':len(replay['descendants']),
            'exact_partitions_checked':True, 'arithmetic_recomputed':arithmetic}


def benchmark_checks():
    result = load('benchmarks')
    source_identity(result,'benchmarks.py')
    fresh = benchmarks.run()
    require({k:fresh[k] for k in result} == result, 'Exact binomial benchmark differs')
    rows = {r['m']:r for r in fresh['benchmarks']}
    require(Fraction(rows[3]['ratio_lower']) > Fraction(199,500) and
            Fraction(rows[5]['ratio_lower']) > Fraction(2,5), 'Benchmark fails selection threshold')
    return {'exact_rational_benchmarks_checked':True}


def main():
    require(__debug__, 'Run without Python optimization')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arithmetic', action='store_true',
                        help='Recompute inequalities using the supplied subdivisions and cover algorithms.')
    parser.add_argument('--case', choices=['common-gap','two-equal-gaps','two-arbitrary-gaps'],
                        help='Select one certificate family; rational benchmarks are always checked.')
    args = parser.parse_args()
    result = {'benchmarks':benchmark_checks()}
    for name, check in [('common-gap',common_checks), ('two-equal-gaps',equal_checks),
                        ('two-arbitrary-gaps',arbitrary_checks)]:
        if args.case is None or args.case == name:
            result[name] = check(args.arithmetic)
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
