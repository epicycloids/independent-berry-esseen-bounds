from fractions import Fraction
from functools import lru_cache
import json
import math
from pathlib import Path
import numpy as np
from flint import arb
from interval_bounds import up, au
from quadratic_bounds import quadratic_cells as previous_cells

def _upper(value):
    value = Fraction(value)
    result = float(value)
    if Fraction(result) < value:
        result = math.nextafter(result, math.inf)
    return result

@lru_cache(maxsize=1)
def load_table():
    data = json.loads((Path(__file__).resolve().parent / 'scalar_table.json').read_text())
    assert data['schema'] == 'joint-complex-scalar-table-v1'
    kappa = Fraction(data['inradius_lower'])
    assert kappa == Fraction(997, 1000)
    prices = tuple((Fraction(x) for x in data['prices']))
    assert all((x > 0 for x in prices))
    cursor = Fraction(0)
    endpoints, constants = ([], [])
    current = [None] * len(prices)
    for row in data['rows']:
        lo, hi = map(Fraction, (row['xlo'], row['xhi']))
        assert lo == cursor and hi > lo
        cursor = hi
        assert len(row['constants']) == len(prices)
        values = list(map(Fraction, row['constants']))
        assert all((value + price >= 0 for value, price in zip(values, prices)))
        current = [value if old is None else max(value, old) for value, old in zip(values, current)]
        endpoint = float(hi)
        assert Fraction(endpoint) == hi
        endpoints.append(endpoint)
        constants.append([_upper(value / kappa) for value in current])
    assert cursor > 0
    result = (np.asarray(endpoints), np.asarray(constants), np.asarray([_upper(price / kappa) for price in prices]))
    for values in result:
        values.flags.writeable = False
    return result

def table_cells(thi, Bhi, dhi, taulo, tauhi):
    previous = previous_cells(thi, Bhi, dhi, taulo, tauhi)
    if Bhi <= 0:
        return previous
    endpoints, constants, prices = load_table()
    x = up(np.asarray(thi) * au(arb(dhi).sqrt()))
    index = np.searchsorted(endpoints, x, side='left')
    eligible = index < len(endpoints)
    selected = constants[np.minimum(index, len(endpoints) - 1)]
    tau_lower = min(Bhi, max(0.0, taulo))
    tau_upper = min(Bhi, max(tau_lower, tauhi))
    tau = np.where(selected <= 0, tau_lower, tau_upper)
    lines = up(up(selected * tau) + up(prices * Bhi))
    candidate = np.min(lines, axis=-1)
    return np.minimum(previous, np.where(eligible, candidate, Bhi))
