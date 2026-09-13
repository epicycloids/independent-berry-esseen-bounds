# Checking the upper-bound calculation

From the repository root, run:

```sh
./verify.sh --no-latex
```

This checks the file hashes, the closed partition of the moment parameters,
the supplied Gaussian interpolation bounds, and the comparison of every
stored upper bound with

$$
U_* = \frac{8104460571768603}{18014398509481984}<0.44988794.
$$

It also checks the endpoint estimates, the
[scalar certificates](../certificates/scalar/README.md), the default
[finite-array certificates](../certificates/finite/README.md), and the
extremizer paper's
[exact rational Gaussian certificate](../paper/extremizers/gaussian_window_certificate.py).
Running `./verify.sh` adds the PDF builds. The locked environment uses
Python 3.12.

## The supplied certificate

After normalizing total variance to one, the calculation uses the
third-moment sum $L$ and three further parameters: the largest summand
variance $d$, the third absolute moment $b$ of that summand, and
$\tau=\sum_i v_i^{3/2}$, where $v_i$ are the individual variances.
For each interval of $L$, a binary tree subdivides a containing region
in $(d,b,\tau)$ into closed boxes. Each final box either has an upper
bound for $\Delta/L$ or is excluded by the moment constraints.

The [cover index](../certificates/cover.json) lists 1,670 closed intervals
of $L$, with exact endpoints, record hashes, and evaluator settings.
Their trees and final boxes are in `certificates/cover.zip`. The index
also identifies the Gaussian interpolation partitions and smoothing
weights in the numbered `certificates/gaussian-NNN.zip` archives. The
[result metadata](../result.json) gives the total counts and the exact
maximum.

A final box records either its own evaluated upper bound or a bound valid
on its entire containing interval. The record identifies the choice.
The certificate contains 143,618 individual values and 92,041 bounds valid
on their entire containing intervals. Within 8,671 final boxes, a further
complete proof covers the feasible third-moment interval or partitions the
largest-variance interval. These proofs retain containing parent bounds,
exact closed children, and the moment inequalities used to exclude empty
children. Their numerical nodes specify fixed evaluator inputs. They may
also use a complete signed threshold partition, a rational split between
two characteristic-function disks, or a forcing bound obtained by convex maximization over the moment
weights, as derived in the paper. Every signed partition includes its outside-tail bound.

The counts of nodes and final boxes in `result.json` describe the outer
moment trees. Nested proof records do not increase those counts. The
`complete_owner_proofs` field counts these nested moment proofs. The
default checks reconstruct both levels of coverage and check the stored
Gaussian and scalar certificates. The commands below recompute the
fixed numerical inputs, including each selected nested proof.

The [paper](../paper/independent-berry-esseen.pdf) gives the enclosure
formulas and arithmetic assumptions. The
[scalar guide](../certificates/scalar/README.md) describes the coefficient
tables used by those formulas.

## Numerical reevaluation

To recompute the first interval:

```sh
uv run --frozen python verify/replay.py --band-index 0
```

An interval's zero-based index is given in the cover index. To recompute
only its first ten final boxes:

```sh
uv run --frozen python verify/replay.py --band-index 0 --first-leaf 0 --max-leaves 10
```

The output reports the interval and the range of boxes evaluated. For
an individual stored value, the new calculation must agree exactly.
For a containing-interval bound, the new value must be no larger than
that bound. Each evaluator retains the precision, smoothing parameters,
and early-stopping level that affect its numerical output. A nested proof
retains the completed subdivisions and evaluated parent bounds, so its
reevaluation does not choose a new subdivision or stop at a new target.

Numerical reevaluation uses one thread. Large intervals can take
substantially longer than the default checks. The evaluators use
[NumPy](https://numpy.org/doc/), [SciPy](https://docs.scipy.org/doc/scipy/),
and [python-flint](https://python-flint.readthedocs.io/) for Arb arithmetic.
The paper states the bounds used for binary64 rounding and integration;
see Johansson's papers on [Arb](https://doi.org/10.1109/TC.2017.2690633)
and [ball arithmetic integration](https://doi.org/10.1007/978-3-319-96418-8_30).
