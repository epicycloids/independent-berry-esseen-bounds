import sys
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace
import full_b_consumer as consumer
import full_b_cover as cover
from flint import arb

def reject(call):
    try:
        call()
    except (ValueError, AssertionError):
        return
    raise AssertionError('Invalid resume was accepted')

def run():
    began = time.process_time()
    digest = consumer.check_helper_source()
    counts = {}
    box = (0.2, 0.4, 0.3, 0.31, 0.4, 0.5)
    full = consumer.full_b_box(0.65, box)
    assert full is not None and full[2] <= box[2] and (full[3] >= box[3])
    assert full == consumer.full_b_box(0.65, (*box[:2], 0.0, 0.65, *box[4:]))
    for count in (1, 2, 4, 8):
        pieces = consumer.closed_b_slices(full, count)
        assert pieces[0][2] == full[2] and pieces[-1][3] == full[3]
        assert all((a[3] == b[2] for a, b in zip(pieces, pieces[1:])))
        assert all((p[:2] + p[4:] == full[:2] + full[4:] for p in pieces))
    point = (*full[:2], 0.25, 0.25, *full[4:])
    assert consumer.closed_b_slices(point, 4) == (point,)
    counts['closed_b_partitions'] = 5
    maximum = consumer.coupled.quadratic_exponential_max([(0.0, 1.0, 0.0)], [(0.0, 1.0)], 0.0, 1.0)
    exact = (arb(1) / 2).sqrt() * (-arb(1) / 2).exp()
    assert arb(maximum) >= exact and maximum < 0.428882
    counts['finite_root_integration'] = 1
    seen = []

    class Source:

        def envelopes(self, Llo, Lhi, box):
            seen.append(('source', tuple(box)))
            return (object(), 0.7, 0.7)

        def integrate(self, B, F, Llo, low):
            return F

    class Rectangles:

        def __init__(self, weight):
            self.weight = weight

        def score(self, Llo, Lhi, box, envelopes):
            seen.append(('reference', tuple(box)))
            return 0.45

    class Weight(consumer.FullBMaximumMixin, Source):
        target = None
        rectangles = None
        rectangle_class = Rectangles
        b_partitions = (1, 2, 4)
    old_tighten = consumer.tighten_envelopes

    def tighten(weight, Llo, Lhi, box, envelopes, mode):
        seen.append(('coupled', tuple(box)))
        return (envelopes[0], 0.6, 0.6)
    try:
        consumer.tighten_envelopes = tighten
        weight = Weight()
        assert weight.box(0.6, 0.65, *box) == 0.45
        assert not any((tag == 'coupled' for tag, b in seen))
        seen.clear()
        for count in weight.b_partitions:
            for piece in consumer.closed_b_slices(full, count):
                assert weight._score_box(0.6, 0.65, piece, full, coupling=True) == 0.45
        sources = [b for tag, b in seen if tag == 'source']
        coupled = [b for tag, b in seen if tag == 'coupled']
        refs = [b for tag, b in seen if tag == 'reference']
        assert sources == coupled and len(sources) == 7
        assert refs == [full] * 14
        assert len({b[2:4] for b in sources}) == 7
        counts['source_coupling_reference_agreement'] = 7

        class Scores(Weight):

            def _score_box(self, Llo, Lhi, piece, reference_box, **kwargs):
                n = len(self.calls)
                self.calls.append(piece)
                return (0.7, 0.4, 0.6, 0.1, 0.2, 0.3, 0.5)[n]
        scored = Scores()
        scored.calls = []
        assert scored.box(0.6, 0.65, *box) == 0.5
        counts['max_over_slices_min_over_portfolios'] = 1

        class Abort(Scores):
            target = 0.45
            b_partitions = (1, 2, 4, 8)

            def _score_box(self, Llo, Lhi, piece, reference_box, **kwargs):
                self.calls.append(piece)
                return 0.7 if len(self.calls) == 1 else 0.5
        aborted = Abort()
        aborted.calls = []
        assert aborted.box(0.6, 0.65, *box) == 0.7 and len(aborted.calls) == 4
        counts['abandoned_partitions_retain_full_bound'] = 1
        assert len(refs) == 14

        class Lazy(Weight):
            coupling = True
            target = 0.46
        seen.clear()
        assert Lazy().box(0.6, 0.65, *box) == 0.45
        assert not any((tag == 'coupled' for tag, b in seen))
        counts['cheap_acceptance_skips_optional_coupling'] = 1
    finally:
        consumer.tighten_envelopes = old_tighten
    core = cover.core
    old_factory, old_time = (core.weight_factory, cover.time)

    class Clock:
        value = -1

        def monotonic(self):
            self.value += 1
            return self.value

    class Stub:
        T, s = (1.0, 1.0)

        def uniform(self, lo, hi):
            return 1.0

        def box(self, lo, hi, *b):
            return 0.35 if b[1] - b[0] < 0.15 and b[5] - b[4] < 0.15 else 0.6

    def factory(*args, **kwargs):
        return lambda tauhi: Stub()
    try:
        core.weight_factory = factory
        cover.time = Clock()
        partial = cover.full_b_cover_band(0.6, 0.61, '0.4', 4, seconds_limit=2.5)
        assert not partial['complete'] and partial['nodes'] == 2
        initial = deepcopy(partial)
        resumed = cover.full_b_cover_band(0.6, 0.61, '0.4', 4, partial, seconds_limit=1000)
        assert partial == initial and resumed['complete']
        whole = cover.full_b_cover_band(0.6, 0.61, '0.4', 4, seconds_limit=1000)
        for key in ('tree_sha256', 'leaf_chain_sha256', 'upper', 'nodes', 'accepted', 'infeasible'):
            assert resumed[key] == whole[key], key
        assert set(core.unpack(whole['tree_zlib_base64'])) <= set('ADTX')
        assert core.replay_band(whole, numerical=True)['nodes'] == whole['nodes']
        bad = deepcopy(partial)
        bad['stack'][0] = (bad['stack'][0][0], 99)
        reject(lambda: cover.full_b_cover_band(0.6, 0.61, '0.4', 4, bad))
        reject(lambda: cover.full_b_cover_band(0.6, 0.61, '0.4', 4, partial, kind=cover.KIND + '_b2'))
        bad = deepcopy(partial)
        bad['tree_zlib_base64'] = core.pack('B')
        bad['tree_sha256'] = hashlib.sha256(b'B').hexdigest()
        bad['nodes'] = 1
        reject(lambda: cover.full_b_cover_band(0.6, 0.61, '0.4', 4, bad))
        counts['closed_tree_nodes'] = whole['nodes']
        counts['resume_rejections'] = 3
        counts['resume_and_frozen_replay'] = True
    finally:
        core.weight_factory, cover.time = (old_factory, old_time)
    return dict(passed=True, checks=counts, helper_sha256=digest, cpu_seconds=time.process_time() - began, scope='Synthetic interface, closed-partition, and resume checks')
