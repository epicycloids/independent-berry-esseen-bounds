import math
from adaptive_b_consumer import partition_upper, depth_for_box

def run():
    full = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    observed = []

    def width_score(box):
        observed.append(box)
        return 0.4 + (box[3] - box[2])
    value, details = partition_upper(full, width_score, 0.45)
    assert value == 0.43125 and details['evaluations'] == 63 and (details['leaves'] == 32)
    assert observed[1][3] == observed[2][2] == 0.5
    for budget in (1, 2, 3, 7, 13, 32):
        value, details = partition_upper(full, width_score, 0.45, max_evaluations=budget)
        assert 0.43125 <= value <= 1.4
        assert details['evaluations'] <= budget
        assert details['evaluations'] == 2 * details['leaves'] - 1
    value, details = partition_upper(full, lambda box: 0.4 if box == full else 1.0, 0.3, max_depth=2)
    assert value == 0.4
    value, details = partition_upper(full, lambda box: None if box[2] >= 0.5 else box[3], 0.6)
    assert value == 0.5 and details['leaves'] == 1
    calls = []
    adjacent = (0.0, 0.0, 1.0, math.nextafter(1.0, math.inf), 0.0, 0.0)
    value, details = partition_upper(adjacent, lambda box: calls.append(box) or 1.0, 0.5)
    assert value == 1.0 and len(calls) == 1
    assert depth_for_box(full) == 1 and depth_for_box(adjacent) == 6
    assert partition_upper(full, lambda box: 0.2, 0.3)[1]['evaluations'] == 1
    rejected = 0
    for kwargs in ({'max_depth': 7}, {'max_depth': True}, {'max_evaluations': 0}, {'max_evaluations': 128}, {'max_evaluations': 3.0}):
        try:
            partition_upper(full, width_score, 0.5, **kwargs)
        except ValueError:
            rejected += 1
    assert rejected == 5
    try:
        partition_upper(full, lambda box: math.nan, 0.5)
    except ArithmeticError:
        rejected += 1
    assert rejected == 6
    return dict(complete_partitions=True, parent_minimum=True, excluded_child=True, representability=True, invalid_inputs_rejected=rejected)
