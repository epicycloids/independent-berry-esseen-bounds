from fractions import Fraction
import math
from flint import arb, ctx
import packed_energy_reference as frozen

def rational(value):
    return value if isinstance(value, Fraction) else Fraction(value)

def exact(value):
    value = rational(value)
    return arb(value.numerator) / arb(value.denominator)

def upper(value):
    if not value.is_finite():
        raise ArithmeticError('A finite outward enclosure is required')
    return 0.0 if value == 0 else math.nextafter(float(value.upper()), math.inf)

class PreparedPrefix:

    def __init__(self, table, base):
        self.precision = ctx.prec
        self.table = table
        compiled = table['small_frequency_prefix']
        zero = frozen.small_prefix_pair(arb(0), compiled, precision_bits=self.precision)
        if zero is None:
            raise ValueError('The small-frequency prefix must include zero')
        self.small_right = exact(rational(compiled['right']))
        self.price = zero[1]
        self.rows = tuple(((base._arb(rational(row['right'])), (row['prefix_Q_upper'], row['prefix_lambda_upper'])) for row in table['rows']))

    def pair(self, argument):
        if ctx.prec != self.precision:
            raise ValueError('Prepared prefix precision differs from active precision')
        if not argument.is_finite():
            raise ValueError('The frequency enclosure must be finite')
        if argument >= 0 and argument <= self.small_right:
            x = arb(argument.upper())
            q = arb(8) * x.sinc() ** 4 / (3 * x.cos() ** 2)
            return (upper(q), self.price)
        if not (argument >= 0 and argument <= 1):
            return None
        for right, pair in self.rows:
            if argument <= right:
                return pair
        raise ValueError('Compiled complete energy table is missing a required prefix')

class PreparedMoment:

    def __init__(self, L, tau):
        self.precision = ctx.prec
        self.Llo, self.Lhi = map(rational, L)
        a, b = map(rational, tau)
        self.a, self.b = (max(Fraction(0), a), min(self.Lhi, b))
        if not (0 < self.Llo <= self.Lhi and self.a <= self.b):
            raise ValueError('Invalid coefficients or moment box')
        self.cache = {}

    def value(self, pair):
        if ctx.prec != self.precision:
            raise ValueError('Prepared moment precision differs from active precision')
        if pair in self.cache:
            return self.cache[pair]
        Q, price = map(rational, pair)
        if Q < 0 or price < 0:
            raise ValueError('Invalid coefficients or moment box')
        candidates = [self.a, self.b]
        if Q > price:
            vertex = (Q - 2 * price) * self.Lhi / (2 * (Q - price))
            if self.a <= vertex <= self.b:
                candidates.append(vertex)
        squares = []
        price_arb = exact(price)
        for t in candidates:
            if t == self.Lhi:
                squares.append(0.0)
            else:
                e = exact(self.Lhi - t)
                squares.append(upper(e * (exact(Q * t) + price_arb * e)))
        result = upper(arb(max(squares)).sqrt())
        self.cache[pair] = result
        return result

class PreparedPacking:
    PRICE = Fraction(9, 8)
    PADDED_ARGUMENT = Fraction(10001, 10000)
    TAU_DECREASING_THRESHOLD = Fraction(15, 23)

    def __init__(self, L, d, tau):
        self.precision = ctx.prec
        Llo, Lhi = map(rational, L)
        dlo, dhi = map(rational, d)
        a, b = map(rational, tau)
        a, b = (max(Fraction(0), a), min(b, Lhi))
        if not (0 < Llo <= Lhi and 0 <= dlo <= dhi <= 1 and (dhi > 0) and (a <= b)):
            raise ValueError('Invalid frequency or moment box')
        packing_globals = frozen.packing_source_upper.__globals__
        packing_globals['derivative_cap_certificate']()
        dhi_arb = exact(dhi)
        self.sqrt_dhi = dhi_arb.sqrt()
        cap = Fraction(upper(dhi_arb * self.sqrt_dhi))
        total = a if a >= self.TAU_DECREASING_THRESHOLD * Lhi else b
        count = total // cap
        residual = total - count * cap
        assert 0 <= residual < cap
        self.cap_root = exact(cap).root(3)
        self.full_weight = exact(count * cap)
        self.residual = bool(residual)
        self.residual_weight = exact(residual) if residual else None
        self.residual_root = exact(residual).root(3) if residual else None
        self.excess = Lhi - a
        self.excess_arb = exact(self.excess)
        self.price_excess = exact(self.PRICE * self.excess)
        self.padded_limit = exact(self.PADDED_ARGUMENT)

    @staticmethod
    def q_ball(x):
        return arb(8) * x.sinc() ** 4 / (3 * x.cos() ** 2)

    def source_at(self, frequency, frequency_square, deletion):
        if ctx.prec != self.precision:
            raise ValueError('Prepared packing precision differs from active precision')
        actual_argument = frequency * self.sqrt_dhi
        if not actual_argument <= 1:
            return None
        padded_argument = frequency * self.cap_root
        if not padded_argument <= self.padded_limit:
            return None
        packed = self.full_weight * self.q_ball(padded_argument)
        if self.residual:
            packed += self.residual_weight * self.q_ball(frequency * self.residual_root)
        square = self.excess_arb * (packed + self.price_excess)
        moment = 0.0 if self.excess == 0 else upper(square.sqrt())
        return upper(frequency_square * deletion * arb(moment) / 2)

class PreparedGeometry:

    def __init__(self, nodes, base):
        self.precision = ctx.prec
        self.parsed = base._nodes(nodes)
        self.rows = tuple(((base._arb(rh), base._arb(rh) ** 2, base._arb(rh - ll), base._arb(ll)) for (ll, lh), (rl, rh) in zip(self.parsed, self.parsed[1:])))

class PreparedBox:

    def __init__(self, L, d, tau, base):
        self.precision = ctx.prec
        self.key = (tuple(map(rational, L)), tuple(map(rational, d)), tuple(map(rational, tau)), self.precision)
        self.L, self.d, self.tau = self.key[:3]
        self.dhi = self.d[1]
        self.sqrt_dhi = base._arb(self.dhi).sqrt()
        self.moment = PreparedMoment(self.L, self.tau)
        self.packing = None
        self.base = base

    def packing_source(self, frequency, frequency_square, deletion):
        if not 0 < self.dhi <= 1:
            return None
        if self.packing is None:
            self.packing = PreparedPacking(self.L, self.d, self.tau)
        return self.packing.source_at(frequency, frequency_square, deletion)

def fast_packed_energy_reference_bounds(nodes, deletion_caps, *, energy_table, L, d, tau, max_argument=Fraction(3, 2), precision_bits=128, prepared_prefix=None, prepared_geometry=None, prepared_box=None):
    if energy_table.get('source_sha256') != frozen.REVIEWED_SHA256:
        raise ValueError('Require a compiled envelope for the reviewed transformed target')
    if rational(energy_table['small_frequency_prefix']['right']) != Fraction(1, 2):
        raise ValueError('Require the full small-frequency prefix')
    base = frozen.reviewed_reference()
    nodes, deletion_caps = (tuple(nodes), tuple(deletion_caps))
    old = base.reference_bounds(nodes, deletion_caps, L=L, d=d, tau=tau, max_argument=max_argument, precision_bits=precision_bits)
    taulo = Fraction(old['taulo_used'])
    radii, endpoints, cells = ([], [0.0], [])
    source_changes = 0
    previous_precision = ctx.prec
    try:
        ctx.prec = max(previous_precision, 128, int(precision_bits))
        if prepared_geometry is None:
            prepared_geometry = PreparedGeometry(nodes, base)
        elif prepared_geometry.precision != ctx.prec or prepared_geometry.parsed != base._nodes(nodes):
            raise ValueError('Prepared geometry differs from grid or precision')
        expected_box = (tuple(map(rational, L)), tuple(map(rational, d)), tuple(map(rational, tau)), ctx.prec)
        if prepared_box is None:
            prepared_box = PreparedBox(L, d, tau, base)
        elif prepared_box.key != expected_box:
            raise ValueError('Prepared box differs from moment inputs or precision')
        if prepared_prefix is not None and (prepared_prefix.precision != ctx.prec or prepared_prefix.table is not energy_table):
            raise ValueError('Prepared prefix differs from table or precision')
        u = base._arb(taulo)
        for j, (geometry, deletion) in enumerate(zip(prepared_geometry.rows, deletion_caps)):
            original = old['cells'][j]
            if not original['eligible']:
                radii.append(math.inf)
                endpoints.append(math.inf)
                cells.append(dict(original))
                continue
            if prepared_prefix is None:
                prepared_prefix = PreparedPrefix(energy_table, base)
            frequency, frequency_square, width, left = geometry
            deletion_arb = base._arb(deletion)
            old_source = original['source_upper']
            source, energy_source = (old_source, None)
            pair = prepared_prefix.pair(frequency * prepared_box.sqrt_dhi)
            if pair is not None:
                moment = prepared_box.moment.value(pair)
                energy_source = base._upper(frequency_square * deletion_arb * arb(moment) / 2)
                source = min(old_source, energy_source)
            packed_source = prepared_box.packing_source(frequency, frequency_square, deletion_arb)
            if packed_source is not None:
                source = min(source, packed_source)
            source_changes += source < old_source
            k = left if taulo == 0 else (left * u).tan() / u
            exp_decay = (-k * width).exp()
            psi = width if k == 0 else -(-k * width).expm1() / k
            previous = endpoints[-1]
            whole = max(previous, base._upper(exp_decay * arb(previous) + arb(source) * psi))
            whole = min(whole, old['E_upper'][j])
            endpoint = min(whole, old['endpoint_upper'][j + 1], base._upper(arb(original['propagation_ratio_upper']) * arb(previous) + arb(source) * psi))
            radii.append(whole)
            endpoints.append(endpoint)
            row = dict(original)
            row.update(source_upper=source, old_source_upper=old_source, energy_source_upper=energy_source, packed_source_upper=packed_source, source_branch='energy' if source < old_source else 'old', E_upper=whole, endpoint_upper=endpoint, old_E_upper=old['E_upper'][j])
            cells.append(row)
    finally:
        ctx.prec = previous_precision
    result = dict(old)
    result.update(status='source-only packed energy/reference wrapper; no signed smoothing certificate', E_upper=radii, endpoint_upper=endpoints, cells=cells, old_E_upper=old['E_upper'], old_endpoint_upper=old['endpoint_upper'], energy_improved_source_cells=int(source_changes), energy_improved_radius_cells=sum((eligible and new < prior for eligible, new, prior in zip(old['eligible'], radii, old['E_upper']))), M_upper_role='Inherited old moment factor; the selected energy source varies by cell', reference_source_sha256=frozen.REFERENCE_SHA256, energy_source_sha256=frozen.REVIEWED_SHA256, energy_table_replayed_per_box=False, small_frequency_source_sha256=energy_table['small_frequency_prefix']['source_sha256'])
    return result
