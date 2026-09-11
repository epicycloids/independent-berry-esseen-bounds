"""Keep the manuscript's numerical summary synchronized with result.json."""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def render(result):
    values = {
        'UpperBound': result['upper_bound'],
        'ExactUpper': '\\frac{'+result['largest_recorded_upper_fraction'].replace('/', '}{')+'}',
        'ReplayedCount': f"{result['arithmetic_replayed_leaves']:,}",
        'SharedBoundCount': f"{result['containing_band_bound_leaves']:,}",
        'SmallUpper': result['small_L_upper_decimal'],
        'VarianceUpper': result['variance_upper_decimal'],
        'LargeUpper': result['large_L_upper_decimal'],
        'BandCount': f"{result['counts']['bands']:,}",
        'NodeCount': f"{result['counts']['nodes']:,}",
        'AcceptedCount': f"{result['counts']['accepted']:,}",
        'InfeasibleCount': f"{result['counts']['infeasible']:,}",
        'RecordedUpper': str(result['largest_recorded_upper']),
    }
    return '% Generated from result.json by verify/paper_values.py.\n' + ''.join(
        '\\newcommand{\\'+name+'}{'+value+'}\n' for name, value in values.items())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='replace the generated TeX values')
    args = parser.parse_args()
    expected = render(json.loads((ROOT/'result.json').read_text()))
    path = ROOT/'paper/values.tex'
    if args.write:
        path.write_text(expected)
        print('Wrote paper/values.tex')
    elif not path.exists() or path.read_text() != expected:
        raise SystemExit('Paper values differ from result.json; review before regenerating')
    else:
        print('Paper values match result.json')
