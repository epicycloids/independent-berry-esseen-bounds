# The Berry–Esseen constant for independent summands

**6 October 2026 · collection version 0.7.0**

For independent centered real random variables with finite third absolute
moments, put

```math
V=\sum_i\mathbb E X_i^2>0,\qquad
L=\frac{\sum_i\mathbb E|X_i|^3}{V^{3/2}},\qquad
\Delta=\sup_x\left|\mathbb P\left(\frac{\sum_iX_i}{\sqrt V}\le x\right)-\Phi(x)\right|.
```

Here $\Phi$ is the standard normal distribution function. The least constant
in $\Delta\le C L$ for every finite family is $C_{\mathrm{ind}}$. This
collection establishes

```math
C_{\mathrm E}=\frac{3+\sqrt{10}}{6\sqrt{2\pi}}
=0.409732\ldots\ \le C_{\mathrm{ind}}<0.44988.
```

[Esseen (1956)](https://doi.org/10.1080/03461238.1956.10414946) proved
the lower bound. The conjecture $C_{\mathrm{ind}}=C_{\mathrm E}$ remains open.

| Paper | Purpose and main results |
| --- | --- |
| [An upper bound on the Berry–Esseen constant for independent summands](paper/independent-berry-esseen.pdf) | The unrestricted bound $C_{\mathrm{ind}}\lt 0.44988$, proved by analytic estimates and a finite interval certificate. |
| [Berry–Esseen bounds for two summands and Bernoulli arrays](paper/finite-arrays/finite-arrays.pdf) | $C_2\lt 49/120$ for two arbitrary summands; the sharp $C_{\mathrm E}$ bound for common-diameter Bernoulli sums, three Bernoulli summands of arbitrary diameters, and specified larger arrays. |
| [The sharp Berry–Esseen bound for independent arrays with small Lyapunov ratio](paper/small-lyapunov/small-lyapunov.pdf) | $\Delta\le C_{\mathrm E}L$ whenever $L\le\exp(-\exp(97))$, for arbitrary independent summands. |
| [Extremizers for the independent Berry–Esseen inequality](paper/extremizers/extremizers.pdf) | Four-atom restrictions on every original marginal of a full maximizer, three-atom representatives preserving individual moments, and restrictions on hypothetical unrestricted extremizers above $C_{\mathrm E}$. |
| [Zero-bias transport and extremal Bernoulli arrays](paper/transport/transport.pdf) | Equality in the additive zero-bias transport bound, a lattice-span inequality, and the sharp asymptotic ratio when a fixed finite block is repeated. |
| [Bernoulli extremizers for the independent Berry–Esseen ratio](paper/bernoulli-structure/bernoulli-structure.pdf) | Classification of common-diameter maximizers, the exact two-Bernoulli optimum, comparison of nearby diameters, and the asymptotic pattern of zero padding in optimized binomial arrays. |
| [Small Lyapunov ratios and the structure of Berry–Esseen maximizers](paper/small-lyapunov-structure/small-lyapunov-structure.pdf) | Common-lattice approximation for changing arrays approaching $C_{\mathrm E}$ as $L\to0$, and exact padded i.i.d. binary structure of every full finite-horizon maximizer with sufficiently small $L$. The structural threshold is ineffective. |

The [overview](docs/overview.md) explains the arguments and dependencies.
The [bounds note](docs/bounds.md) states the constants and the classes covered
by the sharp inequality. The [numerical observations](docs/numerical-observations.md)
describe the stored binomial searches, their coverage, and the structural
conjectures they suggest. Each paper includes a bibliography and TeX source.

The small-Lyapunov argument adapts He and Cheng's method for sufficiently
large i.i.d. samples. Its imported estimates are cited to
[their first version](https://arxiv.org/abs/2609.06358v1); the
[current version](https://arxiv.org/abs/2609.06358v3) presents their i.i.d.
eventual-sharpness theorem. The independent-array argument also uses the
finite paper's common-diameter theorem.

## Check or build

Install [uv](https://docs.astral.sh/uv/) and, to build the PDFs,
[Tectonic](https://tectonic-typesetting.github.io/en-US/). The locked
verification environment uses Python 3.12. From this directory:

```sh
./verify.sh --no-latex
./build.sh
```

The first command checks hashes, stored partitions and bounds, and the
default certificate checks. The [verification guide](verify/README.md),
[finite-array supplement](certificates/finite/README.md), and
[Bernoulli supplement](certificates/structure/README.md) explain what is
checked and how to reevaluate the numerical inputs. The second command
builds all seven papers under `dist/`. Running `./verify.sh` also builds them.

## Citation and reuse

Author: Logan Bell. Cite a manuscript by its author, title, and date
(6 October 2026), with collection version `0.7.0` to identify the files
used. [CITATION.cff](CITATION.cff) supplies collection metadata.
The manuscripts and explanatory text use CC BY 4.0; code, configuration,
and certificate data use MIT. See [LICENSE](LICENSE).
