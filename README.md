# The Berry–Esseen constant for independent summands

**10 September 2026 · collection version 0.3.0.dev0**

For a finite sum of independent centered real random variables with finite
third absolute moments and positive total variance, let $\Delta$ be the
Kolmogorov distance of the standardized sum from the standard normal law,
and let $\ell$ be the normalized sum of third absolute moments. The least
constant in $\Delta\le C\ell$ is denoted by $C_{\mathrm{ind}}$. The open
conjecture is

$$
C_{\mathrm{ind}}=C_{\mathrm E}
=\frac{3+\sqrt{10}}{6\sqrt{2\pi}}
=0.4097321837023963\ldots.
$$

[Esseen (1956)](https://doi.org/10.1080/03461238.1956.10414946) proved
the lower bound. This collection studies an unrestricted upper bound and
the conjectured inequality for finite Bernoulli sums and for general
independent arrays with sufficiently small $\ell$.

| Paper | Main results |
| --- | --- |
| [An upper bound on the Berry–Esseen constant for independent summands](paper/independent-berry-esseen.pdf) | The bound $C_{\mathrm{ind}}<0.454$, using analytic estimates and interval computation. |
| [Sharp Berry–Esseen bounds for Bernoulli sums and reductions for general arrays](paper/finite-arrays/finite-arrays.pdf) | The $C_{\mathrm E}$ inequality for centered two-point summands with a common support diameter and arbitrary probabilities, and for two summands each supported on at most two points with arbitrary diameters; concentration estimates and three-point reductions. |
| [The sharp Berry–Esseen bound for independent arrays with small Lyapunov ratio](paper/small-lyapunov/small-lyapunov.pdf) | A universal $\ell_0>0$ such that $\Delta\le C_{\mathrm E}\ell$ whenever $0<\ell\le\ell_0$, with no restriction on the summand laws. The proof gives no numerical value of $\ell_0$. |

The small-Lyapunov paper adapts the approach of
[He and Cheng (2026)](https://arxiv.org/abs/2609.06358v1) for sufficiently
large i.i.d. samples. It uses the common-diameter theorem from the
finite-results paper, including that theorem's certificate calculations.
The unrestricted conjecture remains open.

The [overview](docs/overview.md) gives reading routes and summarizes the
results and remaining questions. The [bounds note](docs/bounds.md) gives
the normalization and numerical bounds. Each paper has a bibliography and
accompanying TeX source. The [upper-bound verification guide](verify/README.md)
and [finite-results supplement](certificates/finite/README.md) describe the
computations.

## Check or build

Install [uv](https://docs.astral.sh/uv/) and, to build the PDFs,
[Tectonic](https://tectonic-typesetting.github.io/en-US/). The verification
environment uses Python 3.12. From this directory:

```sh
./verify.sh --no-latex
./build.sh
```

The first command checks stored records and runs the default certificate
checks. The verification guides give commands for repeating the interval
calculations.
The second command rebuilds all three papers under `dist/`. Running
`./verify.sh` performs the default checks and the PDF builds. Both scripts
accept `--help`.

## Citation and reuse

Author: Logan Bell. Cite a manuscript by author, title, and date
(10 September 2026); include collection version `0.3.0.dev0` to identify
the files used. [CITATION.cff](CITATION.cff) provides collection metadata.
The manuscripts and explanatory text use CC BY 4.0; code, configuration,
and certificate data use MIT. See [LICENSE](LICENSE).
