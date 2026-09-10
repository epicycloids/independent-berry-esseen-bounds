import hashlib
import importlib.util
from pathlib import Path
import sys
if __package__:
    from .full_b_paths import HERE
else:
    from full_b_paths import HERE
from flint import arb
try:
    import beta_coupled_cf_bounds as coupled
except ModuleNotFoundError as error:
    if error.name != 'beta_coupled_cf_bounds':
        raise
    spec = importlib.util.spec_from_file_location('beta_coupled_cf_bounds', Path(__file__).with_name('bounds.py'))
    coupled = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = coupled
    spec.loader.exec_module(coupled)
from interval_bounds import au
from stable_bounds import feasible_clip
HELPER_SHA256 = 'd5ac41819eae9856947f43155a0bd716d5d28f1a76752f4dc9d52cde9d80de3a'
tighten_envelopes = coupled.tighten_envelopes

def check_helper_source():
    digest = hashlib.sha256(Path(coupled.__file__).read_bytes()).hexdigest()
    if digest != HELPER_SHA256:
        raise ValueError('The coupled CF helper differs from the frozen reviewed source')
    return digest

def full_b_box(Lhi, box):
    if len(box) != 6:
        raise ValueError('A six-coordinate containing box is required')
    dl, dh, _, _, tl, th = box
    return feasible_clip((dl, dh, 0.0, Lhi, tl, th), Lhi)

def closed_b_slices(box, count):
    if count not in (1, 2, 4, 8):
        raise ValueError('Only one, two, four, or eight closed b slices are supported')
    pieces = [tuple(box)]
    for _ in range({1: 0, 2: 1, 4: 2, 8: 3}[count]):
        refined = []
        for piece in pieces:
            mid = (piece[2] + piece[3]) / 2
            if not piece[2] < mid < piece[3]:
                refined.append(piece)
                continue
            left, right = (list(piece), list(piece))
            left[3] = right[2] = mid
            refined.extend((tuple(left), tuple(right)))
        pieces = refined
    return tuple(pieces)

class FullBMaximumMixin:
    cf_mode = 'dual_tyurin'
    rectangle_class = None
    b_partitions = (1,)
    coupling = False

    def envelopes(self, Llo, Lhi, box):
        full = full_b_box(Lhi, box)
        if full is None:
            return None
        old = super().envelopes(Llo, Lhi, full)
        return tighten_envelopes(self, Llo, Lhi, full, old, mode=self.cf_mode) if self.coupling else old

    def box(self, Llo, Lhi, *box):
        full = full_b_box(Lhi, box)
        if full is None:
            return None
        best = float('inf')
        for coupling in (False, True) if self.coupling else (False,):
            for count in self.b_partitions:
                scores = []
                abandoned = False
                for piece in closed_b_slices(full, count):
                    value = self._score_box(Llo, Lhi, piece, full, coupling=coupling)
                    if value is not None:
                        scores.append(value)
                        if count > 1 and self.target is not None and (value >= self.target) and (best < float('inf')):
                            abandoned = True
                            break
                if not scores or abandoned:
                    continue
                best = min(best, max(scores))
                if self.target is not None and best < self.target:
                    return best
        if best == float('inf'):
            raise ArithmeticError('No feasible score in a nonempty full-b box')
        return best

    def _score_box(self, Llo, Lhi, box, reference_box, *, coupling=False):
        old_envelopes = super().envelopes(Llo, Lhi, box)
        if old_envelopes is None:
            return None
        old = self.integrate(*old_envelopes[:2], Llo, old_envelopes[2])
        outside = au(1 / (5 * arb(Llo)))
        if self.target is not None and old < self.target:
            return old
        best = old
        eligible = outside < old and (self.target is None or outside < self.target)
        if eligible:
            if self.rectangles is None:
                if self.rectangle_class is None:
                    raise TypeError('A reviewed prepared rectangle consumer is required')
                self.rectangles = self.rectangle_class(self)
            best = min(best, self.rectangles.score(Llo, Lhi, reference_box, old_envelopes))
        if self.target is not None and best < self.target:
            return best
        if not coupling:
            return best
        envelopes = tighten_envelopes(self, Llo, Lhi, box, old_envelopes, mode=self.cf_mode)
        best = min(best, self.integrate(*envelopes[:2], Llo, envelopes[2]))
        if self.target is not None and best < self.target or outside >= best or (not eligible):
            return best
        return min(best, self.rectangles.score(Llo, Lhi, reference_box, envelopes))
