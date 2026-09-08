# The independent Berry–Esseen constant

**Research draft · 7 September 2026 · version 0.1.0.dev0**

Let $C_{\mathrm{ind}}$ be the smallest universal constant in the classical
Berry–Esseen inequality for finite independent centered real arrays with
positive total variance and finite third absolute moments. This draft
presents a proof candidate for the upper bound

$$
C_{\mathrm E}=\frac{\sqrt{10}+3}{6\sqrt{2\pi}}
=0.40973218370239634299\ldots
\le C_{\mathrm{ind}}<0.474999998<0.475,
$$

where the lower bound is due to Esseen. The proposed upper bound combines
analytic estimates with a completed interval calculation. **Independent
mathematical review remains outstanding.** The conjecture
$C_{\mathrm{ind}}=C_{\mathrm E}$ is unresolved.

## Reading

The paper, [*An upper bound on the Berry–Esseen constant for independent
summands*](paper/independent-berry-esseen.pdf), gives the argument and describes
the computation. Its [LaTeX source](paper/main.tex) and
[bibliography](paper/refs.bib) are included.

| Guide | Contents |
| --- | --- |
| [Overview](docs/overview.md) | Definitions and the main steps of the argument |
| [Bounds and their scope](docs/bounds.md) | Published comparisons and the bounds on the three moment ranges |
| [Remaining questions](docs/questions.md) | Sources of loss and possible refinements |
| [Verification guide](verify/README.md) | Saved records, arithmetic assumptions, checks, and replay commands |

## Computation and replay

The interval calculation covers the middle of the
[three moment ranges](docs/bounds.md#the-three-moment-ranges). It divides
that range into 767 closed intervals, called moment bands. The recorded
partition has 5,577,029 nodes, 2,750,409 accepted leaves, and 38,489
infeasible leaves. The largest recorded upper endpoint is
0.47499999772796003. All accepted-leaf bounds were evaluated in the
original completed calculation.

The supplied replay records cover the full geometric partition and
numerically recompute three bands containing 4,944 accepted leaves. The
remaining accepted-leaf quadratures have not been replayed in full.

[result.json](result.json) summarizes the status and values.
[cover.json](certificates/cover.json) contains the completed bands
and partition tree. The nine [arithmetic source files](verify/kernel/README.md)
are preserved byte-for-byte and identified by SHA-256 digests.

## Check or build

Install [uv](https://docs.astral.sh/uv/) and, to build the PDF,
[Tectonic](https://tectonic-typesetting.github.io/en-US/). From this directory:

```sh
./verify.sh --no-latex   # saved records, scalar enclosures, and rejection tests
./build.sh              # rebuild the paper into dist/
```

`./verify.sh` without an option performs both steps. These commands do not
dispatch cloud computation. The [verification guide](verify/README.md)
documents optional local and Modal replay.

## Citation and reuse

Author: Logan Bell. Cite this as a research draft and include its version;
[CITATION.cff](CITATION.cff) supplies the metadata. Explanatory text and the
manuscript use CC BY 4.0; code, configuration, and certificate data use MIT.
See [LICENSE](LICENSE).
