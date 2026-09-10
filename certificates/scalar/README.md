# Scalar certificate verification

These files record the coefficient bounds used in
[the paper](../../paper/main.tex).

From the archive root, run:

```sh
uv run --frozen python verify/check_scalar_certificates.py
```

The default check regenerates the three exact polynomials and their
Bernstein coefficient bounds, verifies every tangent lower bound for the
logarithmic correction with Arb, checks all saved partitions and exact
rational inequalities derived from their leaf bounds, verifies the analytic
tails, and reconstructs the scalar tables used by the evaluator.
**It uses the saved transcendental enclosures at the leaves without
recomputing them.**

To recompute the saved leaf bounds, run:

```sh
uv run --frozen python verify/check_scalar_certificates.py --replay directional --replay disk --replay energy
```

This recomputes every leaf in the requested families using the supplied Arb
formulas. Adding `--limit N` restricts each affine family to its first N
certificates and the energy family to its first N frequency-cell records.
The output states whether each family's numerical reevaluation is complete.
The directional certificates contain over a million leaves, so their full
reevaluation takes substantially longer than the default checks. All
calculations use one numerical thread and process one affine certificate at
a time. For checks of the moment cover and smoothing calculations, see the
[verification guide](../../verify/README.md).

## Files and mathematical roles

| Input | Content |
|---|---|
| `odd-quadratic.json` | Three exact polynomials and bounds on their Bernstein coefficients, with analytic series and spatial tails for the quadratic bound on the imaginary part of the zero-bias error |
| `direct-log.json` | 2,006 tangent lower bounds at five values of $\lambda$, with exact dyadic coefficients, tangent points and domains |
| `affine-index.json`, `directional.zip` | 61 frequency cells and five slopes, each with 25 directions in the closed upper half-plane: 7,625 certificates with 1,558,869 leaves |
| `affine-index.json`, `disk.zip` | 128 frequency cells and five common slopes: 640 certificates for modulus bounds with 75,046 leaves |
| `energy.json` | 128 frequency cells for the cosine comparison: 25,175 rectangles on frequencies $[0,1/2]$ and 150,453 on $[1/2,1]$, with analytic tail bounds |
| `source/odd_polynomial.py` | Exact rational generator used to regenerate every polynomial coefficient and Bernstein minimum |
| `source/direct_log.py` | Construction of tangent lower bounds at 100-bit precision for the five values of $\lambda$ |
| `source/entire_target.py`, `source/directional.py`, `source/disk.py` | Analytic expression for the normalized zero-bias error, including its value at zero frequency, and Arb formulas for tail checks and leaf reevaluation |
| `source/energy_majorant.py`, `source/energy_small.py` | Quotients used to check the reference-error inequality, including their continuous values at zero frequency, remainder bounds and exact tail polynomials |

The directional and disk certificates bound the same normalized zero-bias
error, through real projections and the complex modulus, respectively. The
energy certificates bound the forcing term in the differential equation
for comparison with the cosine product. At zero frequency, the checker uses
the continuous quotient after removing the factors that vanish there.

The logarithmic bounds use $\lambda\in\{3/4,1,3/2,2,3\}$. The checker
verifies the tangent and endpoint conditions at 192-bit Arb precision;
the evaluator constructs these bounds at 100-bit precision. The affine
zero-bias bounds use the common slopes $5/4,3/2,2,3,4$.

## Exact compact record format

Every rational string denotes its exact value. A finite JSON float, where
present in a runtime table, denotes its exact binary64 value.

Each affine ZIP entry retains its coefficients, correction, closed frequency
interval, precision, tail radius, tail coefficients and complete binary
partition. A split node is `[axis, cut, left_index, right_index]`, with axis
zero for frequency and one for the spatial variable. A directional leaf is
`[lower]`; a disk leaf is `[lower_P, lower_H]`. The root rectangle is the
frequency interval times `[-R,R]` for a directional certificate, or `[0,R]`
for a disk. Child rectangles are reconstructed from exact cuts. Both closed
children cover their parent, including the cut; every node must be reachable
exactly once.

For a directional certificate, each saved lower bound plus the final
correction must be nonnegative. For a disk, the two lower bounds apply at
the saved base correction. If the final additional correction is `e`, the
exact tests are `lower_P + e >= 0` and
`lower_H + 2*e*lower_P + e^2 >= 0`. Both forms have complete analytic tails.
The disk certificate's slope is converted to a common slope using
$\rho=\mathbb E|W|^3\ge1$ for a mean-zero, variance-one variable $W$.
Directional bounds are converted to modulus bounds using a polygon whose
inradius is at least 997/1000; this inclusion is checked exactly.

An energy rectangle is `[xlo, xhi, ylo, yhi, lower]`. The rectangles must
cover the full spatial interval at every frequency boundary and between
consecutive boundaries. The spatial interval is $[0,48]$ for frequencies in
$[0,1/2]$ and $[0,16]$ for $[1/2,1]$. Exact positive tail polynomials complete
the real line. Numerical reevaluation includes the continuous values at
removable singularities.

The affine index identifies the certificates used for each table row.
The checker reconstructs maxima over directions, conversions to common
disk slopes, coefficients for the energy bound, and cumulative maxima over
frequency cells. They must agree exactly with the corresponding rows in
`verify/kernel/scalar_table.json`, `disk_table.json` and
`energy_table_complete.json`.

The checker reports the number of numerically reevaluated leaves or
rectangles separately for each family.
