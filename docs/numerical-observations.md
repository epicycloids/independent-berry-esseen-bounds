# Numerical Studies and Structural Conjectures

Navigation: [Repository Root](../README.md) · [Overview & Reading Guide](overview.md) · [Summary of Bounds](bounds.md) · [Verification Guide](../verify/README.md)

---

For a fixed number of independent centered summands $n$, let

```math
V = \sum_{i=1}^n \mathbb{E} X_i^2 > 0, \qquad
R = \sum_{i=1}^n \mathbb{E}|X_i|^3 < \infty, \qquad
L = \frac{R}{V^{3/2}}, \qquad
\mathcal{C}(X_1, \dots, X_n) = \frac{d_{\mathrm{K}}\!\left(\sum_{i=1}^n X_i / \sqrt{V}, \, Z\right)}{L},
```

where $Z \sim N(0, 1)$ and $d_{\mathrm{K}}$ denotes Kolmogorov distance. Denote by $C_n$ the supremum of $\mathcal{C}(X_1, \dots, X_n)$ over all $n$-tuples of independent centered real random variables (allowing deterministic zero summands).

This guide summarizes the numerical investigations stored in [`data/empirical-observations.json`](../data/empirical-observations.json) and the certified four-summand local maximizer in [`certificates/observations/local-trap.json`](../certificates/observations/local-trap.json), and explains how they motivate the structural theorems and open conjectures of the collection. A comprehensive mathematical treatment of the structural conjecture—including the exact-fiber third-moment functional $\mathcal{M}_n(\mu)$, the separated-support moment and characteristic-function criteria, and the analysis of what remains open—is given in [Manuscript 8](../paper/structural-conjecture/structural-conjecture.pdf) ([TeX](../paper/structural-conjecture/main.tex)).

---

## 1. Binomial Optima and Zero Summands

Consider an array consisting of $m$ identically distributed centered two-point summands and $n - m$ zero summands:

```math
X_i = \frac{B_i - p}{\sqrt{mpq}} \quad (1 \le i \le m), \qquad
X_i = 0 \quad (m < i \le n),
```

where $B_1, \dots, B_m$ are independent $\operatorname{Bern}(p)$ random variables with $0 < p < 1$ and $q = 1 - p$. This array has total variance $V = 1$, common support diameter $h = 1/\sqrt{mpq}$, and Lyapunov ratio $L = (p^2 + q^2)/\sqrt{mpq}$. Writing $K \sim \operatorname{Bin}(m, p)$ and $z_k = (k - mp)/\sqrt{mpq}$, its Berry–Esseen ratio is

```math
R_m(p) = \frac{\sqrt{mpq}}{p^2 + q^2} \max_{0 \le k \le m} \bigl\{ \mathbb{P}(K \le k) - \Phi(z_k), \, \Phi(z_k) - \mathbb{P}(K < k) \bigr\}.
```

Let

```math
B_m = \max_{0 < p < 1} R_m(p), \qquad P_n = \max_{1 \le m \le n} B_m.
```

Because adjoining zero summands leaves both the sum and $L$ unchanged, any ratio $B_m$ achieved by $m$ active summands remains achievable for every $n \ge m$. By Theorem 1.1 of [Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf), $P_n$ is the exact maximum of $\mathcal{C}(X_1, \dots, X_n)$ over all $n$-summand two-point arrays with a common support diameter.

### Numerical evaluation of $B_m$ and $P_n$

The dataset [`data/empirical-observations.json`](../data/empirical-observations.json) records numerical evaluations of $B_m$ and $P_n$ for $1 \le m \le 3000$:

- **Exhaustive threshold evaluation ($1 \le m \le 100$):** A grid of $4096$ probabilities in $(0, 1/2]$ (covering $(0, 1)$ by reflection) evaluates every atom $0 \le k \le m$ and both one-sided CDF branches. Local maxima are bracketed and refined to $60$ decimal places, with stationarity residuals at most $2.61 \times 10^{-59}$.
- **Central-threshold scan ($101 \le m \le 3000$):** Seven central threshold indices $k$ around $mp$ are evaluated for each $m$. For $101 \le m \le 450$, a denser control check evaluates all central threshold cells, both CDF branches, and $17$-point sub-meshes per cell; for ratios exceeding $0.4045$, the bounds of [Schulz (2016)](https://d-nb.info/1197702695/34) exclude noncentral thresholds.

Selected values of the padded binomial envelope $P_n = B_{m_n}$ (with $j_n = n - m_n$ zero summands) are shown below:

| Summands $n$ | Active summands $m_n$ | Zero summands $j_n = n - m_n$ | Optimal probability $p$ | Support diameter $h$ | Ratio $P_n = B_{m_n}$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 0 | 0.390958631281333 | 2.049327 | 0.370351523474379 |
| 2 | 1 | 1 | 0.390958631281333 | 2.049327 | 0.370351523474379 |
| 3 | 3 | 0 | 0.431607387750433 | 1.165657 | 0.398156326353304 |
| 10 | 10 | 0 | 0.429865759681911 | 0.638771 | 0.406507736905768 |
| 50 | 49 | 1 | 0.414895409236786 | 0.289945 | 0.409145180658100 |
| 100 | 99 | 1 | 0.417421147624576 | 0.203806 | 0.409451697372009 |
| 412 | 409 | 3 | 0.418879394858355 | 0.100221 | 0.409665276377194 |
| 413 | 409 | 4 | 0.418879394858355 | 0.100221 | 0.409665276377194 |

Notably, $n = 412$ and $n = 413$ are the first values of $n$ for which $P_n$ uses three and four zero summands, respectively: both retain $m = 409$, whose refined ratio $B_{409} \approx 0.409665276377194310$ exceeds $B_{412} \approx 0.409665275927097300$ by approximately $4.50 \times 10^{-10}$.

---

## 2. Asymptotic Phase Explanation of Zero Padding

Why does setting $j_n = n - m_n \in \{0, 1, 2, 3, 4\}$ summands to zero improve the binomial Berry–Esseen ratio at certain values of $n$?

As $m \to \infty$, the optimizing probability $p$ converges to $p_{\mathrm{E}} = (4 - \sqrt{10})/2$, while the optimal fractional offset $c = mp - k$ between the mean $mp$ and the integer threshold $k$ converges to $c_* = (20 + \sqrt{10})/72$. Because $k$ must be an integer, $mp - k$ cannot simultaneously equal $c_*$ while $p = p_{\mathrm{E}}$ for every integer $m$.

In Section 4 of [Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf), a uniform Edgeworth–Euler–Maclaurin expansion establishes that

```math
B_m = C_{\mathrm{E}} - \frac{Q_*}{m} - \frac{H_* + \alpha \operatorname{dist}(m p_{\mathrm{E}} - \delta_0, \mathbb{Z})^2}{m^2} + O(m^{-3}),
```

where $H_*$ is a constant independent of $m$, and

```math
Q_* = \frac{\phi(0)(867 + 287\sqrt{10})}{25920}, \qquad
\alpha = \frac{\phi(0)(60 + 19\sqrt{10})}{30}, \qquad
\delta_0 = \frac{2590 + 577\sqrt{10}}{12960}.
```

Consequently, when comparing $m = n - j$ against $n$ for bounded $j \ge 0$, choosing $j$ zero summands minimizes the second-order loss

```math
\eta j + \operatorname{dist}(\theta_n + j p_{\mathrm{E}}, \mathbb{Z})^2, \qquad
\eta = \frac{Q_*}{\alpha} = \frac{2510 - 747\sqrt{10}}{8640} \approx 0.01709,
```

where $\theta_n = \operatorname{round}(n p_{\mathrm{E}} - \delta_0) - (n p_{\mathrm{E}} - \delta_0) \in [-1/2, 1/2]$. The linear term $\eta j$ measures the $O(n^{-1})$ penalty for reducing the active sample size from $n$ to $n - j$, while the quadratic term measures the phase mismatch of the integer threshold.

Because $5\eta$ exceeds the maximum of the lower envelope of the first five parabolas $j \in \{0, 1, 2, 3, 4\}$ over $\theta \in [-1/2, 1/2]$, at most four zero summands are ever optimal for sufficiently large $n$. By equidistribution of the irrational rotation $n p_{\mathrm{E}} \bmod 1$, the asymptotic frequencies of $j_n = 0, 1, 2, 3, 4$ are

```math
p_{\mathrm{E}} - \omega, \quad p_{\mathrm{E}} - \omega, \quad \sqrt{10} - 3, \quad \omega, \quad \omega, \qquad
\text{where } \omega = \frac{7773 - 2447\sqrt{10}}{8640} \approx 0.004045.
```

Thus three and four zero summands each occur at roughly $0.4\%$ of all large $n$. Carrying the expansion to $O(m^{-3})$ and applying a quadratic-field norm bound in $\mathbb{Q}(\sqrt{10})$ shows that for all sufficiently large $n$, the maximizing active count $m_n$ is unique and the gaps between successive record active counts belong to $\{2, 3, 5\}$.

---

## 3. Local Maximizers and Structural Diagnostics

### A certified strict local maximizer with unequal diameters

Equal support diameters produce atom collisions in the convolution sum $S = \sum_{i=1}^n h_i(B_i - p_i)$. Perturbing $h = (h_1, \dots, h_n)$ splits a colliding atom into distinct values, creating a directional cusp in the Berry–Esseen ratio. When all active diameters are equal, this cusp locks the diameters together.

However, atom collisions can also occur among *unequal* diameters with rational ratios, creating strict local maxima that are globally suboptimal. Consider $n = 4$ two-point summands with support diameters $h = (1, 1, 1, 2)$, upper probabilities $(p, p, p, r)$ where

```math
p = 0.57275160060602450695\ldots, \qquad
r = 0.33747296608998501783\ldots,
```

and raw threshold $k = 2$. At this threshold, four binary configurations $(1,1,0,0)$, $(1,0,1,0)$, $(0,1,1,0)$, and $(0,0,0,1)$ collide at the sum value $2$. These four vectors form the vertices of a tetrahedron in the hyperplane $\{v \in \mathbb{R}^4 : v \cdot h = 2\}$.

The verifier [`verify/observations/local_trap.py`](../verify/observations/local_trap.py) checks the interval certificate [`certificates/observations/local-trap.json`](../certificates/observations/local-trap.json) using 320-bit Arb ball arithmetic in a box of radius $10^{-28}$ around $(p, r)$, certifying that:

1. The box contains a unique stationary point of the $(p, p, p, r)$ probability family, with negative-definite $4 \times 4$ probability Hessian and Berry–Esseen ratio
   ```math
   \mathcal{C} \in 0.28248894967964153979405802 \pm 5.54 \times 10^{-27}.
   ```
2. The first-order diameter derivative is a positive multiple of $\beta \cdot d - \max_{v} v \cdot d$, where $\beta$ lies strictly in the interior of the tetrahedron of colliding configurations with barycentric weights $(0.31153, 0.31153, 0.31153, 0.06541)$. Consequently, every diameter perturbation $d$ transverse to overall scaling strictly decreases the ratio at linear rate $-c \|d_\perp\|$.
3. Every single-summand variational contact majorant holds globally on $\mathbb{R}$.

Thus this unequal-diameter configuration is a strict local maximizer among four-summand two-point arrays (and satisfies the single-coordinate variational conditions), yet its ratio $0.282489\ldots$ is far below the ratio $\sqrt{2}/4 \approx 0.353553$ of two symmetric Bernoulli summands and two zeros, and below the four-summand optimum $P_4 = B_4 \approx 0.400892$. This example illustrates why the global classification theorems in Manuscripts 5 and 6 combine local variational conditions with global replacement arguments.

### Binomial arrays with a compound-Poisson component

[`data/empirical-observations.json`](../data/empirical-observations.json) also records numerical comparisons for centered sums of the form $S_{\mathrm{raw}} = \operatorname{Bin}(M, p) + g(N_+ - N_-)$, where $N_\pm \sim \operatorname{Pois}(\ell_\pm)$ are independent Poisson variables representing a limit of rare summands of step size $g$, with variance $V = Mpq + g^2(\ell_+ + \ell_-)$ and third-moment sum $R = Mpq(p^2 + q^2) + g^3(\ell_+ + \ell_-)$. Across all tested configurations ($M \in \{1, 2, 3, 5, 10, 20\}$, $g \in \{1/2, 1, 2\}$), every positive-intensity candidate has a strictly smaller ratio than the pure binomial sum ($\ell_+ = \ell_- = 0$) at the same probability $p$, consistent with the elimination of Poisson components in Manuscript 7.

### Comparison with zero-bias transport

By Theorem 1.1 of [Manuscript 4](../paper/transport/transport.pdf), the normalized zero-bias transport ratio

```math
\mathcal{B}_{\mathrm{add}}(X_1, \dots, X_n) = \frac{2V W_1(S, S^*)}{R} \le 1
```

equals $1$ if and only if all nonzero summands are identically distributed two-point variables. While $\mathcal{B}_{\mathrm{add}}(X^{(n)}) \to 1$ along any sequence with $L_n \to 0$ and $\mathcal{C}(X^{(n)}) \to C_{\mathrm{E}}$ ([Manuscript 6](../paper/small-lyapunov-structure/small-lyapunov-structure.pdf), Theorem 1.1), the pointwise inequality $\mathcal{C} \le C_{\mathrm{E}} \mathcal{B}_{\mathrm{add}}$ does not hold for arbitrary laws. Indeed, for $n = 1$, modifying $Z \sim N(0, 1)$ by keeping its normal density outside $[-a, a]$ and replacing its restriction to $[-a, a]$ by three atoms at $\{-a, 0, a\}$ with matching mass, zero mean, and second moment yields $\mathcal{C}(X_a) \sim a/6$ and $\mathcal{B}_{\mathrm{add}}(X_a) \sim 5a^2/18$ as $a \downarrow 0$ ([Manuscript 8](../paper/structural-conjecture/structural-conjecture.pdf), Proposition 3.5), so that $\mathcal{C}(X_a) / \mathcal{B}_{\mathrm{add}}(X_a) \sim 3/(5a) \to \infty$ even though $\mathcal{C}(X_a) \to 0$.

---

## 4. Open Finite-$n$ Conjectures

The theoretical results in Manuscripts 2–8 and the numerical evidence above point to two open conjectures for every fixed number of summands $n \ge 1$ (together with the uniqueness questions analyzed in [Manuscript 8](../paper/structural-conjecture/structural-conjecture.pdf), Section 6):

1. **Finite-$n$ Value Conjecture:** For every $n \ge 1$,
   ```math
   C_n = P_n = \max_{1 \le m \le n} B_m.
   ```
2. **Finite-$n$ Extremizer Conjecture:** For every $n \ge 1$, every independent array $(X_1, \dots, X_n)$ attaining $C_n$ consists of $m \le n$ identically distributed two-point summands and $n - m$ zero summands, up to a common scale factor, reflection, and permutation of coordinates.

What is currently proved toward these conjectures:

- Both conjectures hold for $n = 1$ ([Bentkus and Kirša (1989)](https://doi.org/10.1007/BF00966194)).
- Both conjectures hold at every $n \ge 1$ when restricted to two-point summands with a common support diameter ([Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf), Theorem 1.1), and for $n = 2$ when restricted to two-point summands with arbitrary diameters ([Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf), Theorem 1.2).
- Every global maximizer of $C_n$ has at most four support points per summand and admits an equally optimal representation with at most three support points per summand ([Manuscript 7](../paper/extremizers/extremizers.pdf), Theorem 3.2).
- Both conjectures hold for any global maximizer of $C_n$ whose Lyapunov ratio satisfies $L < \varepsilon_0$ ([Manuscript 6](../paper/small-lyapunov-structure/small-lyapunov-structure.pdf), Theorem 1.2). Because adjoining zero summands increases $n$ without reducing $L$, extending this classification from small $L$ to all $n$ remains open.
- By Theorem 2.3 of [Manuscript 8](../paper/structural-conjecture/structural-conjecture.pdf), the Finite-$n$ Extremizer Conjecture for a maximizer of $C_n$ is equivalent to the statement that its sum $S = X_1 + \cdots + X_n$ has a centered binomial distribution of order $m \le n$.
