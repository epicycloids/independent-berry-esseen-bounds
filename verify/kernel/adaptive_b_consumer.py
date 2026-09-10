import heapq
import math
from full_b_consumer import FullBMaximumMixin, full_b_box

def partition_upper(full, score, target, *, max_depth=6, max_evaluations=127):
    if type(max_depth) is not int or not 0 <= max_depth <= 6:
        raise ValueError('Require an integer depth from zero through six')
    if type(max_evaluations) is not int or not 1 <= max_evaluations <= 127:
        raise ValueError('Require an integer evaluation budget from one through 127')
    initial = score(tuple(full))
    if initial is None:
        return (None, dict(evaluations=1, leaves=0, reason='empty'))
    if not math.isfinite(initial):
        raise ArithmeticError('Nonfinite initial upper')
    heap = [(-initial, 0, tuple(full), 0)]
    evaluations, serial = (1, 1)
    reason = 'target'
    while heap and (target is None or -heap[0][0] >= target):
        negative, _, box, depth = heap[0]
        midpoint = (box[2] + box[3]) / 2
        if depth >= max_depth or not box[2] < midpoint < box[3]:
            reason = 'b depth or representation limit'
            break
        if evaluations + 2 > max_evaluations:
            reason = 'evaluation limit'
            break
        left, right = (list(box), list(box))
        left[3] = right[2] = midpoint
        children = []
        for child in (tuple(left), tuple(right)):
            value = score(child)
            evaluations += 1
            if value is None:
                continue
            if not math.isfinite(value):
                raise ArithmeticError('Nonfinite child upper')
            children.append((-min(-negative, value), serial, child, depth + 1))
            serial += 1
        heapq.heappop(heap)
        for entry in children:
            heapq.heappush(heap, entry)
    if not heap:
        raise ArithmeticError('Every child of a nonempty containing box was excluded')
    return (-heap[0][0], dict(evaluations=evaluations, leaves=len(heap), reason=reason))

def depth_for_box(full):
    dl, dh, bl, bh, tl, th = full
    spatial_width = max(dh - dl, 2 * (th - tl))
    if spatial_width <= 0:
        return 6
    ratio = (bh - bl) / spatial_width
    if ratio <= 0:
        return 0
    return min(6, max(1, 1 + math.ceil(math.log2(ratio))))

class AdaptiveBMaximumMixin(FullBMaximumMixin):
    coupling = False

    def box(self, Llo, Lhi, *box):
        full = full_b_box(Lhi, box)
        if full is None:
            return None
        answer, details = partition_upper(full, lambda piece: self._score_box(Llo, Lhi, piece, full, coupling=False), self.target, max_depth=depth_for_box(full), max_evaluations=127)
        self.last_partition = details
        return answer
