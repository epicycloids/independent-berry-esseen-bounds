from functools import lru_cache
import math
from pathlib import Path
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
for directory in (HERE.parent / 'scalar_table', HERE.parent / 'concave_majorant' / 'quadratic_variant', HERE.parent.parent / '26'):
    if directory.is_dir() and str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
from assemble_disk_table import Q, PRICES, CELLS, CELL_COUNT, DENOMINATOR, DISK_SHA256, TABLE_SCHEMA, require, rationals, parse_json
from interval_bounds import up, au
from scalar_table_bounds import table_cells as previous_cells

def _upper(value):
    value = Q(value)
    result = float(value)
    if Q(result) < value:
        result = math.nextafter(result, math.inf)
    return result

def table_arrays(data):
    require(data['schema'] == TABLE_SCHEMA and data['center_moments'] is True and (data['cumulative'] is True), 'Require the direct centered cumulative disk table.')
    require(type(data['denominator']) is int and data['denominator'] == DENOMINATOR and (type(data['cell_count']) is int) and (data['cell_count'] == CELL_COUNT), 'Unexpected disk grid.')
    require(rationals(data['prices'], len(PRICES), 'table prices') == PRICES, 'Unexpected disk prices.')
    require(data['provenance']['disk_certificate_sha256'] == DISK_SHA256 and data['provenance']['numerical_replay_performed_by_assembler'] is False, 'Unexpected certificate generation or replay scope.')
    require(len(data['rows']) == CELL_COUNT, 'The disk table is incomplete.')
    endpoints, constants, cumulative = ([], [], [None] * len(PRICES))
    for cell, (row, interval) in enumerate(zip(data['rows'], CELLS)):
        require(type(row['cell']) is int and row['cell'] == cell and (rationals((row['xlo'], row['xhi']), 2, 'table interval') == interval), 'A disk row has a wrong index, gap, overlap, or endpoint.')
        local = rationals(row['cell_constants'], len(PRICES), 'cell constants')
        require(all((value + price >= 0 for value, price in zip(local, PRICES))), 'A disk bound is negative at rho=1.')
        cumulative = [value if old is None else max(value, old) for value, old in zip(local, cumulative)]
        require(rationals(row['constants'], len(PRICES), 'cumulative constants') == tuple(cumulative), 'Stored cumulative constants do not equal exact prefix maxima.')
        endpoint = float(interval[1])
        require(Q(endpoint) == interval[1], 'A grid endpoint is not an exact dyadic.')
        endpoints.append(endpoint)
        constants.append([_upper(value) for value in cumulative])
    result = (np.asarray(endpoints), np.asarray(constants), np.asarray([_upper(price) for price in PRICES]))
    for values in result:
        require(np.all(np.isfinite(values)), 'A table coefficient overflows binary64.')
        values.flags.writeable = False
    return result

@lru_cache(maxsize=1)
def load_table():
    return table_arrays(parse_json((HERE / 'disk_table.json').read_bytes()))

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
