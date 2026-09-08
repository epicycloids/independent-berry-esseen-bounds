"""Validate saved replay receipts and their exact computational scope."""
import hashlib
import json
from pathlib import Path

from check_saved import ROOT, load_inputs, require
from scalars import scalar_bounds


def validate_receipt(receipt, mode, cover, result, root=ROOT):
    require(receipt['status'] == 'replay passed' and receipt['mode'] == mode, 'Replay status mismatch')
    require(receipt['cover_sha256'] == result['cover_sha256'], 'Replay input digest mismatch')
    require(receipt['kernel_sha256'] == cover['kernel_sha256'], 'Replay arithmetic digest mismatch')
    wrappers = {name: hashlib.sha256((root/'verify'/name).read_bytes()).hexdigest()
                for name in ('replay.py', 'check_saved.py', 'scalars.py')}
    require(receipt['wrapper_sha256'] == wrappers, 'Replay wrapper digest mismatch')
    require(receipt['tree_geometry_replayed'] is True, 'Missing geometric replay')
    require(receipt['accepted_leaf_quadrature_replayed'] is (mode == 'selected'), 'Replay scope mismatch')
    require(receipt['all_bands_numerically_replayed'] is False, 'Unsupported full numerical replay claim')
    require(receipt['independent_mathematical_review_completed'] is False, 'Unsupported independent review claim')
    expected_runtime = dict(result['recorded_runtime'], kernel_arb_precision_bits=100)
    require(receipt['runtime'] == expected_runtime, 'Replay runtime mismatch')
    expected_scalars = scalar_bounds(cover['domain']['lower'], cover['domain']['upper'])
    require(receipt['scalar_enclosures'] == expected_scalars, 'Replay scalar enclosures differ')
    if mode == 'topology':
        require(receipt['counts'] == result['counts'], 'Geometric replay counts mismatch')
        require(receipt['bands'] is None, 'Unexpected numerical band records')
    else:
        indices = result['selected_replay_band_indices']
        require([r['band_index'] for r in receipt['bands']] == indices, 'Numerical selection mismatch')
        total = {'bands': len(indices), 'nodes': 0, 'accepted': 0, 'infeasible': 0}
        for row, index in zip(receipt['bands'], indices):
            band = cover['bands'][index]
            require(band['method'] == 'stability', 'Selected band has no numerical chain')
            require(all(row[k] == band[k] for k in ('Llo','Lhi','upper','leaf_chain_sha256')),
                    'Numerical replay record mismatch')
            counts = {k: band[k] for k in ('nodes','accepted','infeasible')}
            require(row['counts'] == counts, 'Numerical replay counts mismatch')
            for key, value in counts.items():
                total[key] += value
        require(receipt['counts'] == total, 'Numerical replay total mismatch')
    return receipt['counts']


def check(root=ROOT):
    cover, result = load_inputs(root)
    counts = {}
    for mode in ('topology', 'selected'):
        receipt = json.loads((root/f'certificates/{mode}-replay.json').read_text())
        counts[mode] = validate_receipt(receipt, mode, cover, result, root)
    return {'status': 'saved replay receipts match sources, inputs, and stated scope',
            'counts': counts, 'quadrature_rerun_by_this_check': False,
            'independent_mathematical_review_completed': False}


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
