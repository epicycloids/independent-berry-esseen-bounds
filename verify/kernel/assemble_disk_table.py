import argparse
import base64
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import zlib
Q = Fraction
PRICES = tuple(map(Q, ('5/4', '3/2', '2', '3', '4')))
DENOMINATOR = 64
CELL_COUNT = 128
CELLS = tuple(((Q(j, DENOMINATOR), Q(j + 1, DENOMINATOR)) for j in range(CELL_COUNT)))
DISK_SHA256 = 'b8c68e4fcb25228c0b5b3fc46820a4ac78ff6f03a7caf7a0cf030b6c941f7a81'
CERTIFICATE_SCHEMA = 'zero-bias-disk-certificate-v1'
TABLE_SCHEMA = 'centered-direct-disk-table-v1'
COEFFICIENT_ORDER = ['beta_centered', 'gamma_centered', 'c0', 'c2', 'Lambda']
REQUIRED_SOURCES = {'disk_certificate.py', 'entire_target.py', 'fixed_price_candidates.py', 'centered_disk_pilot.py', 'modal_centered_disk_pilot.py'}

class AssemblyError(ValueError):
    pass

def require(condition, message):
    if not condition:
        raise AssemblyError(message)

def exact(value):
    if isinstance(value, bool):
        raise AssemblyError('A boolean is not a rational endpoint or coefficient.')
    if isinstance(value, (int, str, Q)):
        return Q(value)
    if isinstance(value, float) and math.isfinite(value):
        return Q.from_float(value)
    raise AssemblyError('Require an integer, finite binary64 float, or rational string.')

def rationals(values, length, label):
    require(isinstance(values, (list, tuple)) and len(values) == length, f'{label} must contain {length} entries.')
    return tuple(map(exact, values))

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def valid_hash(value):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None

def hash_map(value, label):
    require(isinstance(value, dict) and value, f'{label} must be a nonempty map.')
    require(all((isinstance(k, str) and valid_hash(v) for k, v in value.items())), f'{label} has an invalid source identity.')
    return value

def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f'Duplicate JSON object key {key!r}.')
        result[key] = value
    return result

def parse_json(raw):

    def invalid(value):
        raise AssemblyError(f'Nonfinite JSON number {value}.')
    return json.loads(raw, object_pairs_hook=unique_object, parse_constant=invalid)

def decode_certificate(row):
    compressed = base64.b64decode(row['certificate_zlib_base64'], validate=True)
    decoder = zlib.decompressobj()
    raw = decoder.decompress(compressed) + decoder.flush()
    require(decoder.eof and (not decoder.unused_data) and (not decoder.unconsumed_tail), 'Certificate compression is truncated or contains trailing data.')
    sha = digest(raw)
    require(sha == row['certificate_sha256'], 'Certificate payload SHA-256 mismatch.')
    return (parse_json(raw), sha, digest(compressed))

def convert_affine(intercept, source_price, target_price):
    intercept, source_price, target_price = map(exact, (intercept, source_price, target_price))
    require(0 < source_price <= target_price, 'Ineligible affine source price.')
    return intercept + source_price - target_price

def check_tail(coefficients, correction, interval, radius, recorded):
    beta, gamma, c0, c2, lam = coefficients
    x0, x1 = interval
    if lam > 1:
        cubic, quadratic = (lam - 1, c2 - abs(gamma))
        linear = -(x1 * x1 / 2 + abs(beta))
        constant = c0 + correction + abs(gamma)
    else:
        require(x0 > 0, 'A zero-touching centered certificate requires source price > 1.')
        bound_gamma = max(abs(gamma - 2 / x0), abs(gamma - 2 / x1))
        cubic, quadratic = (lam, c2 - bound_gamma)
        linear = -(4 / (x0 * x0) + abs(beta + 2))
        constant = c0 + correction + bound_gamma - 2 / x0
    expected = (constant + linear * radius + quadratic * radius ** 2 + cubic * radius ** 3, linear + 2 * quadratic * radius + 3 * cubic * radius ** 2, quadratic + 3 * cubic * radius, cubic)
    require(rationals(recorded, 4, 'tail coefficients') == expected and min(expected) >= 0, 'Exact corrected tail coefficients disagree or fail nonnegativity.')

def check_partition(certificate, interval, radius, extra):
    nodes = certificate['nodes']
    require(isinstance(nodes, list) and nodes, 'Missing compact partition payload.')
    require(type(certificate['node_count']) is int and certificate['node_count'] == len(nodes), 'Certificate node count disagrees with its payload.')
    seen, leaves = (set(), 0)
    pending = [(0, (*interval, Q(0), radius))]
    while pending:
        index, rectangle = pending.pop()
        require(type(index) is int and 0 <= index < len(nodes) and (index not in seen), 'Partition has an invalid, repeated, or cyclic node.')
        seen.add(index)
        node = nodes[index]
        actual = (*rationals(node['x'], 2, 'node x interval'), *rationals(node['y'], 2, 'node y interval'))
        require(actual == rectangle, 'Partition rectangle disagrees with its parent.')
        children = node['children']
        if children is None:
            lower_p, lower_h = (exact(node['lower_p']), exact(node['lower_h']))
            require(lower_p + extra >= 0 and lower_h + 2 * extra * lower_p + extra ** 2 >= 0, 'Saved leaf bounds do not imply corrected positivity.')
            leaves += 1
            continue
        require(isinstance(children, list) and len(children) == 2, 'Every internal partition node needs two children.')
        x0, x1, y0, y1 = rectangle
        cut, axis = (exact(node['split_at']), node['split_axis'])
        if axis == 'x' and x0 < cut < x1:
            rectangles = ((x0, cut, y0, y1), (cut, x1, y0, y1))
        elif axis == 'y' and y0 < cut < y1:
            rectangles = ((x0, x1, y0, cut), (x0, x1, cut, y1))
        else:
            raise AssemblyError('A split does not partition its parent rectangle.')
        pending.extend(zip(children, rectangles))
    require(len(seen) == len(nodes), 'The certificate has unreachable partition nodes.')
    require(type(certificate['leaf_count']) is int and certificate['leaf_count'] == leaves, 'Certificate leaf count disagrees with the partition.')
    return (len(nodes), leaves)

def check_quality(certificate, row, replay, interval, radius, correction):
    names = {'correction_tolerance', 'required_correction_lower', 'correction_gap_upper', 'correction_tolerance_met', 'correction_lower_witness'}
    require(names <= certificate.keys(), 'Incomplete correction-quality metadata.')
    tolerance, lower, gap = (exact(certificate[k]) for k in ('correction_tolerance', 'required_correction_lower', 'correction_gap_upper'))
    require(tolerance >= 0 and 0 <= lower <= correction and (gap == correction - lower), 'Correction-quality identities are inconsistent.')
    flag = certificate['correction_tolerance_met']
    require(type(flag) is bool and flag == (gap <= tolerance), 'Invalid correction-tolerance flag.')
    witness = certificate['correction_lower_witness']
    wx, wy, wl = (exact(witness[k]) for k in ('x', 'y', 'lower'))
    require(interval[0] <= wx <= interval[1] and 0 <= wy <= radius and (wl == lower), 'Correction witness is outside the compact domain or has inconsistent metadata.')
    for summary in (row, replay):
        require(exact(summary['correction_gap_upper']) == gap and type(summary['correction_tolerance_met']) is bool and (summary['correction_tolerance_met'] == flag), 'Correction-quality summary differs from the hashed payload.')
    require(exact(replay['required_correction_lower']) == lower, 'Replay correction-lower summary differs from the payload.')
    require(certificate['stop_reason'] != 'correction_tolerance' or flag, 'Correction-tolerance stopping requires its stated tolerance to hold.')

@dataclass(frozen=True)
class Candidate:
    coefficients: tuple
    correction: Q
    intercept: Q
    sha256: str

def validate_row(row, cell):
    lo, hi = CELLS[cell]
    require(rationals((row['xlo'], row['xhi']), 2, 'row interval') == (lo, hi), "Row interval differs from the receipt's cell.")
    price = exact(row['price'])
    require(price in PRICES, 'Unexpected common price.')
    require(type(row['complete']) is bool, 'Row completeness must be an actual boolean.')
    if not row['complete']:
        return (price, None)
    certificate, sha, compressed_sha = decode_certificate(row)
    require(certificate['schema'] == CERTIFICATE_SCHEMA and certificate['certified'] is True, 'Unknown or uncertified disk payload.')
    require(certificate['center_moments'] is True and certificate['coefficient_order'] == COEFFICIENT_ORDER, 'Explicit centered mode and centered coefficient names are required.')
    require(rationals(certificate['x_interval'], 2, 'certificate interval') == (lo, hi), 'Certificate interval differs from the row.')
    coefficients = rationals(certificate['coefficients'], 5, 'disk coefficients')
    beta, gamma, c0, c2, lam = coefficients
    require(0 < lam <= price and (lo > 0 or lam > 1), 'Ineligible disk source price.')
    proposal = row['proposal']
    centered = proposal['centered']
    require(rationals([centered[k] for k in ('beta_centered', 'gamma_centered', 'c0', 'c2', 'price')], 5, 'proposal coefficients') == coefficients, 'Centered proposal coefficients differ from the certified coefficients.')
    require(exact(proposal['requested_price']) == price and proposal['source_price_eligible'] is True and (exact(proposal['price_gap']) == price - lam), 'Proposal price identity or eligibility is inconsistent.')
    midpoint = (lo + hi) / 2
    require(exact(proposal['x']) == exact(centered['x']) == midpoint, 'Proposal frequency is not the recorded cell midpoint.')
    original = proposal['candidate']
    require(exact(original['x']) == midpoint and tuple((exact(original[k]) for k in ('c0', 'c2', 'price'))) == (c0, c2, lam), 'Original and centered proposal identities disagree.')
    correction, base, requested = (exact(certificate[k]) for k in ('correction', 'base_correction', 'requested_correction'))
    require(0 <= requested <= base <= correction, 'Invalid correction or requested/base relation.')
    extra = correction - base
    require(exact(certificate['additional_correction']) == extra, 'Additional correction identity is inconsistent.')
    require(rationals(certificate['corrected_coefficients'], 5, 'corrected coefficients') == (beta, gamma, c0 + correction, c2, lam), 'Corrected coefficients are inconsistent.')
    intercept = c0 + c2 + correction
    require(exact(certificate['affine_intercept']) == intercept and exact(certificate['price']) == lam, 'Certified affine identity is inconsistent.')
    converted = convert_affine(intercept, lam, price)
    require(exact(row['converted_intercept']) == converted, 'Converted affine intercept is inconsistent.')
    require(converted + price >= 0, 'Affine bound is negative at rho=1.')
    goal = certificate['correction_goal_met']
    require(type(goal) is bool and goal == (correction <= requested), 'Invalid correction-goal flag.')
    require(exact(row['correction']) == correction and row['stop_reason'] == certificate['stop_reason'], 'Row correction or stopping identity differs from the payload.')
    radius = exact(certificate['tail_radius'])
    require(radius >= 1, 'Invalid tail radius.')
    check_tail(coefficients, correction, (lo, hi), radius, certificate['tail_coefficients'])
    node_count, leaf_count = check_partition(certificate, (lo, hi), radius, extra)
    replay = row['replay']
    require(replay['replayed'] is True and replay['center_moments'] is True, 'A successful centered worker replay is required.')
    require(type(row['node_count']) is int and row['node_count'] == node_count and (type(replay['node_count']) is int) and (replay['node_count'] == node_count) and (type(replay['leaf_count']) is int) and (replay['leaf_count'] == leaf_count), 'Worker replay/row counts differ from the hashed partition.')
    require(type(certificate['precision_bits']) is int and certificate['precision_bits'] >= 64 and (type(replay['precision_bits']) is int) and (replay['precision_bits'] >= certificate['precision_bits']), 'Invalid replay precision.')
    require(exact(replay['correction']) == correction and exact(replay['affine_intercept']) == intercept and (exact(replay['price']) == lam) and (exact(replay['tail_radius']) == radius), 'Worker replay identities differ from the hashed payload.')
    check_quality(certificate, row, replay, (lo, hi), radius, correction)
    return (price, (Candidate(coefficients, correction, converted, sha), compressed_sha))

class TableAssembler:

    def __init__(self):
        self.rows = {}
        self.sources = {}
        self.certificate_origins = {}
        self.accepted_occurrences = 0
        self.ignored_incomplete_rows = 0

    def add_document(self, data, receipt, *, name, file_sha256, receipt_name, receipt_sha256):
        try:
            require(valid_hash(file_sha256) and valid_hash(receipt_sha256), 'Invalid input file hash.')
            path, receipt_path = (PurePosixPath(name), PurePosixPath(receipt_name))
            require(not path.is_absolute() and '..' not in path.parts and (receipt_path == path.with_suffix('.call.json')), 'Invalid relative result/receipt pair.')
            match = re.fullmatch('(?:centered_disk_pilot|disk_cells)_([0-9a-f]{12})\\.json', path.name)
            require(match is not None and receipt['reservation_id'] == match.group(1), 'Receipt reservation does not match its result filename.')
            require(isinstance(receipt['function_call_id'], str) and receipt['function_call_id'].startswith('fc-'), 'Receipt lacks a worker function-call identity.')
            observed = hash_map(data['source_sha256'], 'worker sources')
            expected = hash_map(receipt['expected_sources'], 'expected sources')
            require(observed == expected and REQUIRED_SOURCES <= observed.keys(), 'Receipt/worker sources differ or omit a required module.')
            cell = data['cell']
            require(type(cell) is int and 0 <= cell < CELL_COUNT and (type(receipt['cell']) is int) and (receipt['cell'] == cell), 'Invalid or mismatched receipt cell.')
            require(type(data['denominator']) is int and data['denominator'] == DENOMINATOR and (type(receipt['denominator']) is int) and (receipt['denominator'] == DENOMINATOR), 'Invalid or mismatched grid denominator.')
            record = dict(sha256=file_sha256, receipt=receipt_name, receipt_sha256=receipt_sha256, worker_source_sha256=observed, receipt_metadata={k: v for k, v in receipt.items() if k != 'expected_sources'}, source_identity_checked=True, eligible=observed['disk_certificate.py'] == DISK_SHA256)
            if not record['eligible']:
                record['excluded_reason'] = 'disk certificate source is outside the frozen version'
                require(name not in self.sources or self.sources[name] == record, 'Conflicting duplicate filename.')
                self.sources[name] = record
                return
            require(isinstance(data['rows'], list), 'Worker rows must be a list.')
            checked = [validate_row(row, cell) for row in data['rows']]
            prices = [price for price, _ in checked]
            require(len(prices) == len(set(prices)), 'A worker document repeats a target price.')
            complete = set(prices) == set(PRICES) and all((candidate is not None for _, candidate in checked))
            require(type(data['complete']) is bool and data['complete'] == complete, 'Outer completeness disagrees with its row identities.')
            record.update(complete=complete, accepted_rows=sum((value is not None for _, value in checked)), incomplete_rows=sum((value is None for _, value in checked)))
            require(name not in self.sources or self.sources[name] == record, 'Conflicting duplicate filename.')
            if name in self.sources:
                return
            self.sources[name] = record
            for row_index, (price, evidence) in enumerate(checked):
                if evidence is None:
                    self.ignored_incomplete_rows += 1
                    continue
                candidate, compressed_sha = evidence
                alternatives = self.rows.setdefault((cell, price), {})
                require(candidate.sha256 not in alternatives or alternatives[candidate.sha256] == candidate, 'Conflicting certificate identity.')
                alternatives[candidate.sha256] = candidate
                origin = dict(file=name, row=row_index, compressed_sha256=compressed_sha)
                self.certificate_origins.setdefault(candidate.sha256, []).append(origin)
                self.accepted_occurrences += 1
        except (KeyError, TypeError, ValueError, ArithmeticError, zlib.error) as error:
            raise AssemblyError(f'{name}: {error}') from error

    def add_file(self, path, *, source_root):
        path, source_root = (Path(path).resolve(), Path(source_root).resolve())
        receipt_path = path.with_suffix('.call.json')
        raw, receipt_raw = (path.read_bytes(), receipt_path.read_bytes())
        self.add_document(parse_json(raw), parse_json(receipt_raw), name=path.relative_to(source_root).as_posix(), file_sha256=digest(raw), receipt_name=receipt_path.relative_to(source_root).as_posix(), receipt_sha256=digest(receipt_raw))

    def missing(self):
        return [dict(cell=cell, xlo=str(CELLS[cell][0]), xhi=str(CELLS[cell][1]), price=str(price)) for cell in range(CELL_COUNT) for price in PRICES if (cell, price) not in self.rows]

    def summary(self):
        return dict(files_checked=len(self.sources), eligible_files=sum((r['eligible'] for r in self.sources.values())), accepted_occurrences=self.accepted_occurrences, complete_cell_prices=len(self.rows), missing_cell_prices=CELL_COUNT * len(PRICES) - len(self.rows), unique_payloads=len(self.certificate_origins), ignored_incomplete_rows=self.ignored_incomplete_rows, numerical_replay_performed_by_assembler=False)

    def finish(self):
        missing = self.missing()
        require(not missing, f'Incomplete disk table: {len(missing)} missing cell/price entries; first {missing[:10]}.')
        cumulative = [None] * len(PRICES)
        output = []
        for cell, (lo, hi) in enumerate(CELLS):
            constants, origins = ([], [])
            for j, price in enumerate(PRICES):
                alternatives = self.rows[cell, price]
                chosen = min(alternatives.values(), key=lambda item: (item.intercept, item.sha256))
                constants.append(chosen.intercept)
                cumulative[j] = chosen.intercept if cumulative[j] is None else max(cumulative[j], chosen.intercept)
                origins.append(dict(price=str(price), certificate_sha256=chosen.sha256, alternative_certificates=sorted(alternatives)))
            output.append(dict(cell=cell, xlo=str(lo), xhi=str(hi), cell_constants=list(map(str, constants)), constants=list(map(str, cumulative)), origins=origins))
        return dict(schema=TABLE_SCHEMA, center_moments=True, denominator=DENOMINATOR, cell_count=CELL_COUNT, prices=list(map(str, PRICES)), cumulative=True, rows=output, provenance=dict(disk_certificate_sha256=DISK_SHA256, assembly_source_sha256=digest(Path(__file__).read_bytes()), files=dict(sorted(self.sources.items())), certificate_origins=dict(sorted(self.certificate_origins.items())), validation='Exact identities, hashes, topology, and saved rational consequences checked; numerical replay is recorded worker evidence.', numerical_replay_performed_by_assembler=False, source_identity_implies_audit_status=False, **{k: v for k, v in self.summary().items() if k != 'numerical_replay_performed_by_assembler'}))

def discover_files(round27):
    round27 = Path(round27)
    return sorted((path for directory, prefix in ((round27, 'centered_disk_pilot_'), (round27 / 'disk_table' / 'raw', 'disk_cells_')) for path in directory.glob(prefix + '*.json') if not path.name.endswith('.call.json')))

def assemble_files(paths, *, source_root=None, progress=None):
    paths = sorted(set((Path(path).resolve() for path in paths)))
    require(paths, 'No worker result files supplied.')
    if source_root is None:
        source_root = Path(os.path.commonpath([str(path.parent) for path in paths]))
    assembler = TableAssembler()
    for path in paths:
        assembler.add_file(path, source_root=source_root)
        if progress:
            progress(dict(file=str(path), **assembler.summary()))
    return assembler.finish()
