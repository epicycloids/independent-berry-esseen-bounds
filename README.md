# The Berry–Esseen constant for independent summands

**Research draft · 10 September 2026 · version 0.2.0.dev0**

This draft proposes the upper bound

$$
C_{\mathrm E}=\frac{\sqrt{10}+3}{6\sqrt{2\pi}}
=0.40973218370239634299\ldots
\le C_{\mathrm{ind}}<0.454
$$

where $C_{\mathrm{ind}}$ is the Berry–Esseen constant for finite families of
independent centered real random variables with finite third absolute
moments. Esseen proved the lower bound. The argument for the upper bound
combines analytic inequalities with an interval calculation. **Independent
mathematical review remains outstanding, and the interval calculation has
not been recomputed in full.** The conjecture
$C_{\mathrm{ind}}=C_{\mathrm E}$ is open.

Read the paper, [*An upper bound on the Berry–Esseen constant for independent
summands*](paper/independent-berry-esseen.pdf), for the definitions and argument.
Its [LaTeX source](paper/main.tex) and [bibliography](paper/refs.bib) are included.

| Guide | Contents |
| --- | --- |
| [Overview](docs/overview.md) | How the moment estimates lead to a uniform bound |
| [Bounds](docs/bounds.md) | Normalization, comparison with a published bound, and estimates on three ranges |
| [Verification](verify/README.md) | Certificate contents, completed checks, arithmetic assumptions, and commands |

The interval calculation covers 1,683 closed intervals of the normalized
third-moment sum. Its largest recorded upper bound is approximately
0.45399999743026986. The complete partition has been checked, and three
selected intervals have been numerically reevaluated.
[result.json](result.json) gives the numerical summary and verification status.

## Check or build

Install [uv](https://docs.astral.sh/uv/) and, to build the PDF,
[Tectonic](https://tectonic-typesetting.github.io/en-US/). From this directory:

```sh
./verify.sh --no-latex
./build.sh
```

The first command checks the stored partitions, scalar certificates and
numerical records. The second rebuilds the paper into `dist/`.
`./verify.sh` without an option performs both steps. See the
[verification guide](verify/README.md) for the scope of each check.

## Citation and reuse

Author: Logan Bell. Cite this as a research draft and include its version;
[CITATION.cff](CITATION.cff) supplies the metadata. The manuscript and
explanatory text use CC BY 4.0; code, configuration, and certificate data
use MIT. See [LICENSE](LICENSE).
