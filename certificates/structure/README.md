# Bernoulli extremizers: computational proofs

The [Bernoulli manuscript](../../paper/bernoulli-structure/main.tex) uses
three interval calculations. The probability and diameter variables range
over closed domains; zero masses and coincident support points are handled
by the analytic boundary arguments in the paper.

| Calculation | Supplied evidence |
| --- | --- |
| Curvature with coefficient $199/500$ | Complete deterministic subdivisions: 550,563 processed boxes and 275,283 accepted boxes for the one-exception inequality; 13,641 processed boxes and 6,822 accepted boxes for the two-coordinate inequalities. |
| Two nonconstant Bernoulli coordinates of equal diameter | Complete interval covers of the diagonal and the full probability square, using the ratio threshold $37/100$ and the gradient and curvature conditions in the paper. |
| Two Bernoulli coordinates of arbitrary diameters | Exact forests covering the three threshold cubes: 37,287 source leaves and 37,421 leaves in a separately formulated arithmetic check. Every original acceptance condition is checked, with no unresolved region. |

The $199/500$ inequalities also hold at the larger coefficient $2/5$
by monotonicity. [benchmarks.json](benchmarks.json) supplies exact lower
bounds at explicit three- and five-coordinate binomial laws. These bounds
provide the competitiveness thresholds used in the proof.

The records and mathematical formulas are in this directory and
[verify/structure/](../../verify/structure/check.py). The arbitrary-diameter
forests are stored losslessly in their JSON records. Their labels identify
the specific necessary inequality that excludes each box; they are
defined in the verifier and derived in the manuscript.

From the repository root, run:

```sh
uv run --frozen python verify/structure/check.py
```

This checks source hashes, completion and counts, exact rational benchmarks,
formula controls, and the entire stored arbitrary-diameter partition,
including each arithmetic-check subpartition. It does not reevaluate every
interval inequality. The common-diameter records supply counts and a
deterministic subdivision algorithm; their full subdivisions are reconstructed
by the arithmetic command below.

```sh
uv run --frozen python verify/structure/check.py --arithmetic
```

This recomputes the two common-diameter covers and checks every supplied
arbitrary-diameter descendant with independently arranged derivative formulas.
It runs on one process in the locked Python environment. Reaching a time,
node or depth limit before completion is a verification failure. Use
`--case common-gap`, `--case two-equal-gaps`, or
`--case two-arbitrary-gaps` to select one family. Both implementations use
Arb through python-flint; independence of formulas does not imply use of a
different arithmetic library.

These computations support binary-law statements. They do not impose
binary support on arbitrary independent laws or settle the unrestricted
two-summand problem.
