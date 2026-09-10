# Checking the upper-bound calculation

From the repository root, run:

```sh
./verify.sh --no-latex
```

This checks file hashes, reconstructs the closed moment partition, checks
the supplied Gaussian interpolation partitions, and compares recorded
upper bounds with the exact target `227/500`. It also checks the analytic
endpoint estimates with Arb, the [scalar certificates](../certificates/scalar/README.md),
and the saved reevaluation records for three selected intervals. The same
command includes the default checks in the
[finite-results supplement](../certificates/finite/README.md).

These checks use the stored integral bounds at the final parameter boxes.
The numerical reevaluation commands below recompute them. Running
`./verify.sh` also builds all three PDFs; `./verify.sh --help` lists the
options. The locked verification environment uses Python 3.12 and is
managed by [uv](https://docs.astral.sh/uv/).

## Calculation and recorded checks

After normalizing total variance to one, the calculation uses three further
moment parameters: the largest summand variance $d$, that summand's third
absolute moment $b$, and $\tau=\sum_i v_i^{3/2}$, where $v_i$ are the summand
variances. For each interval of $L$, it subdivides a containing region in
these parameters into closed boxes. A final box is an `accepted leaf` when
its recorded upper bound for $\Delta/L$ meets the stated target; boxes that
cannot contain realizable moments are excluded.

The cover contains 1,683 closed intervals of $L$, called `bands` in the
files, and 200,239 accepted leaves. The 100 records labeled `mixed` use
two subdivision trees: the second covers every unresolved box left by the
first. Each tree has its own recorded evaluator and acceptance target.

The [Gaussian partition file](../certificates/gaussian-partitions.json)
supplies node enclosures, interpolation intervals, and cutoff corrections
for the mixed records. Other intervals regenerate their Gaussian arrays
during numerical reevaluation.

The [paper's account of the computation](../paper/independent-berry-esseen.pdf)
states the arithmetic assumptions and the enclosure formulas. The
[result metadata](../result.json) records the counts, maximum, and
evaluation methods. The supplied numerical reevaluation records cover
1,346 accepted leaves in three intervals, traversing 2,689 nodes.

## Numerical reevaluation

To recompute the three selected intervals:

```sh
uv run --frozen python verify/replay.py --selected
```

To recompute one interval, supply its zero-based index in
[the cover](../certificates/cover.json). For example, for the first interval:

```sh
uv run --frozen python verify/replay.py --band-index 0
```

Each tree must reproduce its recorded leaf checksum, maximum, and
maximizing box. Mixed records retain separate acceptance targets of
0.453 and 0.454; other records use their own stored targets and evaluators.
The calculation uses one numerical thread and a separate process for
each interval. Large intervals can take substantial time. These commands
are separate from the default checks.

The evaluators use [NumPy](https://numpy.org/doc/),
[SciPy](https://docs.scipy.org/doc/scipy/), and
[python-flint](https://python-flint.readthedocs.io/), which provides bindings
to Arb. The paper states the assumptions for outward rounding, binary64
arithmetic, and Arb integration. See Johansson's papers on
[Arb](https://doi.org/10.1109/TC.2017.2690633) and
[ball arithmetic integration](https://doi.org/10.1007/978-3-319-96418-8_30).
