"""Geometric or numerical replay with the preserved arithmetic kernel."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time

import flint
import numpy as np
import scipy
from check_saved import ROOT, load_inputs, validate_cover
from scalars import scalar_bounds


def run(mode='topology', band_index=None):
    began = time.monotonic()
    cover, result = load_inputs()
    validate_cover(cover, result)
    sys.path.insert(0, str(ROOT/'verify/kernel'))
    from gaussian_batch import install_factory
    from stable_cover import replay_band
    install_factory()
    if mode == 'topology':
        indices = range(len(cover['bands']))
        numerical = False
    elif mode == 'selected':
        indices = result['selected_replay_band_indices']
        numerical = True
    elif mode == 'band' and type(band_index) is int and 0 <= band_index < len(cover['bands']):
        indices = [band_index]
        numerical = True
    else:
        raise ValueError('Choose topology, selected, or band with a valid zero-based index')
    counts = Counter({'bands': 0, 'nodes': 0, 'accepted': 0, 'infeasible': 0})
    records = []
    for index in indices:
        band = cover['bands'][index]
        checked = replay_band(band, numerical=numerical)
        counts.update(checked)
        counts['bands'] += 1
        records.append({'band_index': index, 'Llo': band['Llo'], 'Lhi': band['Lhi'],
                        'counts': checked, 'upper': band['upper'],
                        'leaf_chain_sha256': band.get('leaf_chain_sha256')})
    if mode == 'topology':
        assert dict(counts) == result['counts']
    enclosures = scalar_bounds(cover['domain']['lower'], cover['domain']['upper'])
    wrapper = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
               for p in [Path(__file__), ROOT/'verify/check_saved.py', ROOT/'verify/scalars.py']}
    return {'status': 'replay passed', 'mode': mode, 'counts': dict(counts),
        'bands': records if numerical else None,
        'cover_sha256': result['cover_sha256'], 'kernel_sha256': cover['kernel_sha256'],
        'wrapper_sha256': wrapper, 'scalar_enclosures': enclosures,
        'runtime': {'python': sys.version.split()[0], 'numpy': np.__version__,
                    'python-flint': flint.__version__, 'scipy': scipy.__version__,
                    'kernel_arb_precision_bits': flint.ctx.prec},
        'tree_geometry_replayed': True, 'accepted_leaf_quadrature_replayed': numerical,
        'all_bands_numerically_replayed': False, 'independent_mathematical_review_completed': False,
        'elapsed_seconds': time.monotonic()-began}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['topology','selected','band'], default='topology')
    parser.add_argument('--band-index', type=int)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = json.dumps(run(args.mode, args.band_index), indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
    else:
        print(output, end='')
