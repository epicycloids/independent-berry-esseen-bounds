# The Berry–Esseen Constant for Independent Summands

**6 October 2026 · Collection version 0.7.0**

Navigation: [Overview & Reading Guide](docs/overview.md) · [Summary of Bounds](docs/bounds.md) · [Numerical Studies & Conjectures](docs/numerical-observations.md) · [Verification Guide](verify/README.md)

---

Let $X_1, \dots, X_n$ be independent centered real random variables with positive total variance and finite third absolute moments. Set

```math
V = \sum_{i=1}^n \mathbb{E} X_i^2 > 0, \qquad
R = \sum_{i=1}^n \mathbb{E}|X_i|^3 < \infty, \qquad
L = \frac{R}{V^{3/2}}, \qquad
\Delta = \sup_{x \in \mathbb{R}} \left|\mathbb{P}\!\left(\frac{\sum_{i=1}^n X_i}{\sqrt{V}} \le x\right) - \Phi(x)\right|,
```

where $\Phi$ is the standard normal distribution function and $L$ is the Lyapunov ratio. The universal Berry–Esseen constant for independent summands is

```math
C_{\mathrm{ind}} = \sup \frac{\Delta}{L},
```

where the supremum is taken over all $n \ge 1$ and all such arrays $(X_1, \dots, X_n)$.

In 1956, [Esseen](https://doi.org/10.1080/03461238.1956.10414946) showed that sums of identically distributed two-point random variables with probabilities $p_{\mathrm{E}} = (4 - \sqrt{10})/2$ and $1 - p_{\mathrm{E}}$ achieve the limiting ratio

```math
C_{\mathrm{E}} = \frac{3 + \sqrt{10}}{6\sqrt{2\pi}} = 0.40973218370239634299\ldots,
```

establishing the lower bound $C_{\mathrm{ind}} \ge C_{\mathrm{E}}$. Whether $C_{\mathrm{ind}} = C_{\mathrm{E}}$ holds for all independent arrays remains open.

This repository presents a collection of eight manuscripts, together with computer-assisted proof certificates and verification suites, establishing the unconditional upper bound

```math
C_{\mathrm{E}} \le C_{\mathrm{ind}} < 0.44988,
```

proving the sharp inequality $\Delta \le C_{\mathrm{E}} L$ in several finite and asymptotic regimes, and characterizing the structure of extremal distributions.

## Manuscripts

Each manuscript includes its compiled PDF, LaTeX source, and bibliography under `paper/`.

| Manuscript | Major results | Principal methods |
| --- | --- | --- |
| 1. [An upper bound on the Berry–Esseen constant for independent summands](paper/independent-berry-esseen.pdf) ([TeX](paper/main.tex)) | Unconditional bound $C_{\mathrm{ind}} \le U_* < 0.44988$ for all finite independent arrays. | Zero-bias characteristic-function differential inequality, signed comparison with a product of cosines, Prawitz smoothing, and rigorous interval subdivision over moment parameters. |
| 2. [Berry–Esseen bounds for two summands and Bernoulli arrays](paper/finite-arrays/finite-arrays.pdf) ([TeX](paper/finite-arrays/main.tex)) | $C_2 < 49/120 < C_{\mathrm{E}}$ for two arbitrary summands; $\Delta \le C_{\mathrm{E}} L$ for common-diameter two-point arrays, for at most three two-point summands, and for arrays with variance concentrated in two summands. | Quadratic moment duality (reducing summands to at most three atoms), variational contact exclusions, symmetric polynomial reductions for two-point probabilities, and interval certificates. |
| 3. [The sharp Berry–Esseen bound for independent arrays with small Lyapunov ratio](paper/small-lyapunov/small-lyapunov.pdf) ([TeX](paper/small-lyapunov/main.tex)) | Sharp inequality $\Delta \le C_{\mathrm{E}} L$ for every independent array with $L \le \exp(-\exp(97))$. | Penalized counterexample selection, uniform jitter smoothing, quantitative stability in Esseen's moment inequality, two-point support clustering, and conditional Fourier comparisons adapting [He and Cheng (2026)](https://arxiv.org/abs/2609.06358v1). |
| 4. [Zero-bias transport and extremal Bernoulli arrays](paper/transport/transport.pdf) ([TeX](paper/transport/main.tex)) | Equality in $2V W_1(S, S^*) \le R$ holds if and only if all nonzero summands share a single two-point law; sharp limit $\lim_{m \to \infty} C_m \le C_{\mathrm{E}}$ for repeated blocks, with equality only for Esseen two-point summands. | Fixed-variance Choquet representation over three-point extreme laws, conditional-variance identity for common-diameter Bernoulli arrays, and a sharp lattice span–transport inequality. |
| 5. [Bernoulli extremizers for the independent Berry–Esseen ratio](paper/bernoulli-structure/bernoulli-structure.pdf) ([TeX](paper/bernoulli-structure/main.tex)) | Every $n$-summand common-diameter maximizer consists of i.i.d. two-point summands and zeros; $C_{2,\mathrm{bin}} = C_1 \approx 0.370352$; local common-diameter replacement strictly improves the ratio; $O(m^{-2})$ binomial phase expansion shows at most four zero summands are eventually optimal. | Elementary-symmetric probability fibers, certified curvature bounds, rejective-sampling conditional $L^1$ estimates, and seventh-order Edgeworth–Euler–Maclaurin expansions with quadratic-field separation. |
| 6. [Small Lyapunov ratios and the structure of Berry–Esseen maximizers](paper/small-lyapunov-structure/small-lyapunov-structure.pdf) ([TeX](paper/small-lyapunov-structure/main.tex)) | Arrays with $L \to 0$ and $\Delta/L \to C_{\mathrm{E}}$ concentrate on a common lattice and saturate zero-bias transport; every global maximizer with $L < \varepsilon_0$ (for an ineffective absolute $\varepsilon_0 > 0$) consists of i.i.d. two-point summands and zeros. | Stability analysis of Shevtsova's smoothing inequality, joint local limit theorems for lattice counts and rounding residuals, relative conditional first-moment bounds, and discrete Green-function concavity. |
| 7. [Extremizers for the independent Berry–Esseen inequality](paper/extremizers/extremizers.pdf) ([TeX](paper/extremizers/main.tex)) | Every summand of a maximizer has at most four support points and admits a moment-preserving three-point replacement; summands of unrestricted maximizers are convolution-indecomposable; if $C_{\mathrm{ind}} > C_{\mathrm{E}}$, compactified maximizers have no Gaussian or Poisson components, are nonatomic when infinite, and obey weighted variance summability for $L \le 18481/25000$. | Second-order probability variations, Cox–Kemperman third-moment superadditivity, Lévy–Khintchine compactification, contact-propagation Hessians, and multi-scale grouping with randomized lattice rounding. |
| 8. [On the structural conjecture for extremizers of the independent Berry–Esseen inequality](paper/structural-conjecture/structural-conjecture.pdf) ([TeX](paper/structural-conjecture/main.tex)) | Separates sum-law selection from minimum third-moment factorization ($\mathcal{M}_n(\mu)$); proves that a maximizer of $C_n$ is a padded i.i.d. two-point array if and only if its sum law is binomial of order $m \le n$; establishes a sharp separated-support moment family and an endogenous characteristic-function criterion ($\mathcal{J} = 0$); synthesizes proved regimes, certified local obstructions, and the remaining open global selection step. | Exact third-moment splitting identity, convex-order fiber compactness, probability generating functions, size-biased moment sum-of-squares identities, and endogenous Fourier coherence. |

## Documentation and Verification Guides

| Guide | Contents |
| --- | --- |
| [Overview & Reading Guide](docs/overview.md) | Mathematical narrative connecting the eight manuscripts, proof strategies, dependency structure, and open problems. |
| [Summary of Bounds](docs/bounds.md) | Explicit constants, moment regimes, and classes of independent arrays where $\Delta \le C_{\mathrm{E}} L$ is established. |
| [Numerical Studies & Conjectures](docs/numerical-observations.md) | High-precision binomial optimizations through $m = 3000$, zero-summand phase patterns, a certified suboptimal four-summand local maximizer, and open finite-$n$ conjectures. |
| [Verification Guide](verify/README.md) | Overview of the verification suite, moment-domain partition certificates, and commands for numerical reevaluation. |
| [Scalar Certificates](certificates/scalar/README.md) | Verification of polynomial, logarithmic, directional, disk, and cosine-comparison scalar bounds for Manuscript 1. |
| [Finite-Array Certificates](certificates/finite/README.md) | Verification of two-summand contact exclusions, three-Bernoulli partitions, and Bernstein odds certificates for Manuscript 2. |
| [Bernoulli Structure Certificates](certificates/structure/README.md) | Verification of common-diameter curvature covers and two-Bernoulli stationary exclusions for Manuscript 5. |
| [Interval Evaluation Kernel](verify/kernel/README.md) | Architecture of the outward-rounded interval and ball-arithmetic evaluators under `verify/kernel/`. |

## Running Verification and Building the Manuscripts

The verification suite requires [uv](https://docs.astral.sh/uv/) (which manages the pinned Python 3.12 environment, including `numpy`, `scipy`, `python-flint` for Arb ball arithmetic, and `sympy`). Building the PDF manuscripts from LaTeX source additionally requires [Tectonic](https://tectonic-typesetting.github.io/en-US/).

From the repository root:

```sh
./verify.sh --no-latex   # Verify file hashes, partitions, and default certificates
./build.sh               # Compile the manuscripts into dist/
./verify.sh              # Run all default checks and compile the manuscripts
```

Details on optional full arithmetic reevaluation of the leaf enclosures are provided in the [Verification Guide](verify/README.md), [Scalar Certificates](certificates/scalar/README.md), [Finite-Array Certificates](certificates/finite/README.md), and [Bernoulli Structure Certificates](certificates/structure/README.md).

## Citation and License

Author: Logan Bell. Cite individual manuscripts by author, title, and manuscript date (6 October 2026 for Manuscripts 1–7; 18 October 2026 for Manuscript 8), referencing collection version `0.7.0`. Repository metadata is provided in [CITATION.cff](CITATION.cff). Manuscripts and documentation are licensed under CC BY 4.0; code, configuration, and certificate data are licensed under the MIT License. See [LICENSE](LICENSE).
