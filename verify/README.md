# Checking the supplied calculation

Run `./verify.sh --no-latex` at the repository root. The command reconstructs
the closed moment partition, checks source hashes and supplied Gaussian
interpolation partitions, compares all recorded upper bounds with the exact
target `227/500`, and checks the analytic endpoint estimates with Arb. It
also checks the [scalar certificates](../certificates/scalar/README.md) and
the saved numerical reevaluation of three selected intervals. Run
`./verify.sh` to include the PDF build.

The cover contains 1,683 closed intervals of L, called bands in the files,
and 200,239 accepted leaves. The 100 records labeled `mixed` use two trees:
the second covers every unresolved box left by the first. Each tree is
checked with its recorded evaluator and acceptance target.

The [Gaussian partition file](../certificates/gaussian-partitions.json)
supplies node enclosures, interpolation intervals and cutoff corrections
for the mixed records. Other intervals regenerate their Gaussian arrays
during numerical reevaluation.

The saved numerical reevaluation reproduced 1,346 accepted leaves in three
intervals, across 2,689 nodes. The default command checks those saved records
without repeating their quadratures. A full numerical reevaluation of
accepted leaves and Gaussian quadrature nodes remains outstanding, as does
independent mathematical review. The [result metadata](../result.json)
retains `global_proof: false`. The [paper](../paper/main.tex) states the
analytic bounds and arithmetic assumptions.

To recompute the three selected intervals, run:

```sh
uv run --frozen python verify/replay.py --selected
```

To recompute one interval, use `--band-index` with its zero-based index in
[the cover](../certificates/cover.json). Each tree must reproduce its recorded
leaf checksum, maximum and maximizing box. Mixed records retain separate
acceptance targets of 0.453 and 0.454; other records use their own stored
targets and evaluators. The calculation uses one numerical thread and a
separate process for each interval. Large intervals can take substantial
time. These calculations are optional and are not part of the default check.

The evaluators use [NumPy](https://numpy.org/doc/),
[SciPy](https://docs.scipy.org/doc/scipy/) and
[python-flint](https://python-flint.readthedocs.io/), which provides bindings
to Arb; see Johansson's papers on
[Arb](https://doi.org/10.1109/TC.2017.2690633) and
[ball arithmetic integration](https://doi.org/10.1007/978-3-319-96418-8_30).
