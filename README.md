# The Berry–Esseen constant for independent summands

**11 September 2026 · collection version 0.4.0**

For independent centered real random variables with finite third absolute
moments, put

$$
B^2=\sum_i\mathbb E X_i^2>0,\qquad
L=\frac{\sum_i\mathbb E|X_i|^3}{B^3},\qquad
\Delta=\sup_x\left|\mathbb P\left\{\frac{\sum_iX_i}{B}\le x\right\}-\Phi(x)\right|.
$$

Here $\Phi$ is the standard normal distribution function. The least
constant in $\Delta\le C L$ for every such finite family is
$C_{\mathrm{ind}}$. This collection proves

$$
C_{\mathrm E}=\frac{3+\sqrt{10}}{6\sqrt{2\pi}}
=0.409732\ldots\ \le C_{\mathrm{ind}}<0.44995.
$$

[Esseen (1956)](https://doi.org/10.1080/03461238.1956.10414946) proved
the lower bound. The conjecture $C_{\mathrm{ind}}=C_{\mathrm E}$ remains open.

| Paper | Results |
| --- | --- |
| [An upper bound on the Berry–Esseen constant for independent summands](paper/independent-berry-esseen.pdf) | $C_{\mathrm{ind}}<0.44995$, with analytic estimates and an interval certificate. |
| [Berry–Esseen bounds for two summands and Bernoulli arrays](paper/finite-arrays/finite-arrays.pdf) | $C_2<49/120$ for two arbitrary summands; the $C_{\mathrm E}$ bound for common-diameter Bernoulli sums, three Bernoulli summands of arbitrary diameters, and specified classes of larger arrays. |
| [The sharp Berry–Esseen bound for independent arrays with small Lyapunov ratio](paper/small-lyapunov/small-lyapunov.pdf) | A universal $L_0>0$ such that $\Delta\le C_{\mathrm E}L$ whenever $0<L\le L_0$, for arbitrary independent summands. |
| [Extremizers for the independent Berry–Esseen inequality](paper/extremizers/extremizers.pdf) | Variational constraints on extremal distributions and exclusions for classes of infinite sums. |

The small-Lyapunov paper adapts the method of
[He and Cheng (2026)](https://arxiv.org/abs/2609.06358v1) for sufficiently
large i.i.d. samples and uses the finite paper's common-diameter theorem.
Its proof gives no numerical value of $L_0$.

The [overview](docs/overview.md) explains the arguments and dependencies.
The [bounds note](docs/bounds.md) collects the constants and the precise
classes covered by the sharp inequality. Each paper includes a bibliography
and TeX source.

## Check or build

Install [uv](https://docs.astral.sh/uv/) and, to build the PDFs,
[Tectonic](https://tectonic-typesetting.github.io/en-US/). The verification
environment uses Python 3.12. From this directory:

```sh
./verify.sh --no-latex
./build.sh
```

The first command checks the supplied records and runs the default
certificate checks. The [upper-bound guide](verify/README.md) and
[finite-results supplement](certificates/finite/README.md) distinguish
these checks from full numerical reevaluation and give the corresponding
commands. The second command builds all four papers under `dist/`.
Running `./verify.sh` performs the default checks and the PDF builds.

## Citation and reuse

Author: Logan Bell. Cite a manuscript by its author, title, and date
(11 September 2026), with collection version `0.4.0` to identify the files
used. [CITATION.cff](CITATION.cff) supplies collection metadata.
The manuscripts and explanatory text use CC BY 4.0; code, configuration,
and certificate data use MIT. See [LICENSE](LICENSE).
