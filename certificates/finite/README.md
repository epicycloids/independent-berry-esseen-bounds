# Computational supplement for the finite results

The [finite-results paper](../../paper/finite-arrays/finite-arrays.pdf)
uses five primary certificate calculations. The paper states their
inequalities, parameter domains, and analytic complements.

| Calculation | Use in the paper | Method and recorded primary result |
| --- | --- | --- |
| Curvature inequality for a binomial sum plus one Bernoulli variable | Common-diameter theorem | Three parametrized domains; 448,835 processed boxes and 224,419 accepted boxes at 192-bit Arb precision. |
| Scalar inequality used in the maximal-atom estimate | Atom and short-interval estimates | 129,024 rational boxes covering the order and skewness parameters. The recorded lower bound for the gap exceeds $0.0194747706542$. |
| Inequality for products of Bernoulli odds | Support-endpoint estimate | Exact Bernstein coefficients in $\mathbb Q(\sqrt5)$: 1,281 coefficients in the first table, with one specified zero, and 70 strictly positive coefficients in the second. |
| Interior threshold $(0,1)$ for two Bernoulli summands | Two-summand theorem | 2,231 processed boxes and no unresolved boxes. |
| Interior threshold $(1,0)$ for two Bernoulli summands | Two-summand theorem | 3,509 processed boxes and no unresolved boxes. |

The common-diameter proof uses the curvature inequality and the two-summand
theorem. The latter uses the support-endpoint estimate and both interior
calculations. The scalar atom inequality is used in the separate atom and
short-interval estimates.

[results.json](results.json) records the primary calculations and their
internal checks. Separate evaluators checked the scalar atom inequality at
256 and 384 bits, the odds coefficients exactly, and the first interior
two-summand inequality at 128, 192, and 256 bits. The curvature calculation
was rerun with its primary program. The second
interior inequality was rerun at 256 bits, with separate exact formula and
boundary checks. The record identifies each method and its result.

## Default checks

From the repository root:

```sh
uv run --frozen python verify/finite/check.py
```

This checks program hashes, constants, selected interval enclosures,
analytic boundary estimates, exact four-atom probabilities, reflection,
and degenerate cases. It also tests that omitted pieces of the scalar
partition, exhausted subdivision budgets, and fixed perturbations that
make inequalities false are rejected.

These checks are included in `./verify.sh --no-latex`. Add `--verbose`
to the finite-checker command to see each result.

## Complete calculations

To rerun all five primary certificate programs in sequence:

```sh
uv run --frozen python verify/finite/check.py --full
```

To include the supplied separate implementations, higher-precision
checks, and additional failure controls:

```sh
uv run --frozen python verify/finite/check.py --full --independent
```

Complete calculations can take substantially longer than the default
checks. Each interval program covers its stated domain and succeeds
only when every remaining box has a strictly positive bound. Exhausting
a depth or box limit is a failure. The algebraic programs reconstruct
the Bernstein coefficients and determine their signs by exact arithmetic.

The sources are grouped by their mathematical use:

- [Binomial-plus-one-Bernoulli curvature](../../verify/finite/common_diameter/singleton_certificate.py).
- [Scalar maximal-atom inequality](../../verify/finite/scalar_atom/compact_certificate.py)
  and its [separate evaluator](../../verify/finite/scalar_atom/independent_certificate.py).
- [Product-odds inequality](../../verify/finite/support_odds/exact_certificate.py)
  and its [separate algebraic calculation](../../verify/finite/support_odds/independent_certificate.py).
- [First interior threshold](../../verify/finite/two_summands/certificate.py),
  [second interior threshold](../../verify/finite/two_summands/row10_certificate.py),
  and [exact boundary and reflection checks](../../verify/finite/two_summands/independent_frame.py).

The two interior thresholds are indexed by the Bernoulli outcomes
$(0,1)$ and $(1,0)$ in the programs. When their locations coincide, their
probabilities must be combined before evaluating the distribution
function. The manuscript and the exact checks include that case.

The recorded primary runs used Python 3.14, python-flint 0.8.0, and
SymPy 1.14.0. The locked environment uses Python 3.12 with the same
python-flint and SymPy versions. The separate interval implementations
check formulas and domain coverage; both use the Arb arithmetic library.

For the scalar maximal-atom inequality, the analytic treatment of orders
above 128 ends with a Bessel-function inequality at one endpoint. The
function `bessel_tail_control` in the separate scalar evaluator checks
that inequality as part of the proof. It runs in the default checks and
before every complete calculation.

Arb is described by
[Johansson (2017)](https://doi.org/10.1109/TC.2017.2690633).
The exact algebraic calculations use
[SymPy](https://doi.org/10.7717/peerj-cs.103).
