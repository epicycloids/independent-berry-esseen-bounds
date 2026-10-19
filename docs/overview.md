# Mathematical Overview and Reading Guide

Navigation: [Repository Root](../README.md) · [Summary of Bounds](bounds.md) · [Numerical Studies & Conjectures](numerical-observations.md) · [Verification Guide](../verify/README.md)

---

## 1. The Problem and Central Open Questions

Let $X_1, \dots, X_n$ be independent centered real random variables with positive total variance $V = \sum_{i=1}^n \mathbb{E} X_i^2 > 0$ and finite third absolute moments $R = \sum_{i=1}^n \mathbb{E}|X_i|^3 < \infty$. Write

```math
S = \sum_{i=1}^n X_i, \qquad
L = \frac{R}{V^{3/2}}, \qquad
\Delta = \sup_{x \in \mathbb{R}} \left|\mathbb{P}\!\left(\frac{S}{\sqrt{V}} \le x\right) - \Phi(x)\right|,
```

where $\Phi$ is the standard normal distribution function. The classical Berry–Esseen theorem asserts that the ratio $\Delta / L$ is bounded above by a universal constant $C_{\mathrm{ind}}$ across all finite independent arrays. Similarly, for each fixed integer $n \ge 1$, let $C_n$ denote the supremum of $\Delta / L$ over arrays of at most $n$ independent summands (permitting deterministic zero summands), so that $C_1 \le C_2 \le \cdots$ and $C_{\mathrm{ind}} = \sup_{n \ge 1} C_n$.

In 1956, [Esseen](https://doi.org/10.1080/03461238.1956.10414946) determined the largest asymptotic value of $\Delta / L$ among identically distributed summands. Taking $X_i$ to be identically distributed centered two-point random variables with upper atom probability $p_{\mathrm{E}} = (4 - \sqrt{10})/2$, one obtains $\lim_{n \to \infty} \Delta / L = C_{\mathrm{E}}$, where

```math
C_{\mathrm{E}} = \frac{3 + \sqrt{10}}{6\sqrt{2\pi}} = 0.40973218370239634299\ldots.
```

Consequently, $C_{\mathrm{ind}} \ge C_{\mathrm{E}}$. Three natural questions organize the study of the independent Berry–Esseen inequality:

1. **The Sharp Constant Conjecture:** Does $\Delta \le C_{\mathrm{E}} L$ hold for every finite independent array, i.e., is $C_{\mathrm{ind}} = C_{\mathrm{E}}$?
2. **The Finite-$n$ Value Conjecture:** For each $n \ge 1$, let $B_m$ denote the maximum of $\Delta / L$ over $m$ identically distributed two-point summands, and let $P_n = \max_{1 \le m \le n} B_m$ be the maximum achieved by $m \le n$ identically distributed two-point summands together with $n - m$ zero summands. Does $C_n = P_n$ hold for every $n \ge 1$? (Since $B_m < C_{\mathrm{E}}$ for every $m$ by [Schulz (2016)](https://d-nb.info/1197702695/34) and $\lim_{m \to \infty} B_m = C_{\mathrm{E}}$, the identity $C_n = P_n$ for all $n$ would immediately imply $C_{\mathrm{ind}} = C_{\mathrm{E}}$.)
3. **The Extremizer Classification Question:** For a fixed $n \ge 1$, is every array attaining $C_n$ composed of $m \le n$ identically distributed two-point summands and $n - m$ zero summands, up to a common scale factor, reflection, and coordinate permutation?

## 2. Architecture of the Collection

The seven manuscripts in this collection attack these questions from three complementary directions:

- **Quantitative upper bounds (Manuscripts 1–3):** Establishing the unconditional bound $C_{\mathrm{ind}} < 0.44988$ for all arrays, proving $C_2 < 49/120 < C_{\mathrm{E}}$ for two arbitrary summands, proving $\Delta \le C_{\mathrm{E}} L$ for common-diameter two-point arrays and for three arbitrary two-point summands, and proving $\Delta \le C_{\mathrm{E}} L$ unconditionally whenever $L \le \exp(-\exp(97))$.
- **Extremal structure for two-point arrays and small Lyapunov ratios (Manuscripts 4–6):** Characterizing equality in the zero-bias Wasserstein transport inequality, proving that common-diameter two-point maximizers are identically distributed, analyzing local diameter perturbations and the $O(m^{-2})$ phase expansion governing optimal zero summands, and proving that global maximizers with sufficiently small $L$ are exactly i.i.d. two-point arrays with zero summands.
- **Variational theory of general extremizers (Manuscript 7):** Proving that every summand of any maximizer is supported on at most four points (and can be replaced by a three-point summand with the same variance and third absolute moment), establishing convolution indecomposability, and restricting hypothetical counterexamples above $C_{\mathrm{E}}$.

| Suggested reading order | Manuscript | Role in the collection |
| --- | --- | --- |
| Standalone quantitative bound | [1. Unrestricted upper bound](../paper/independent-berry-esseen.pdf) | Proves $C_{\mathrm{ind}} < 0.44988$ for all finite independent arrays. |
| Finite-array foundations | [2. Two summands and Bernoulli arrays](../paper/finite-arrays/finite-arrays.pdf) | Proves $C_2 < 49/120$ and $\Delta \le C_{\mathrm{E}} L$ for common-diameter Bernoulli arrays (used in Manuscripts 3 and 7). |
| Small-$L$ sharp bound | [3. Small Lyapunov ratio](../paper/small-lyapunov/small-lyapunov.pdf) | Proves $\Delta \le C_{\mathrm{E}} L$ for $L \le \exp(-\exp(97))$ (used in Manuscript 7). |
| Transport rigidity | [4. Zero-bias transport](../paper/transport/transport.pdf) | Characterizes equality in $2V W_1(S, S^*) \le R$ and in fixed-block asymptotic limits. |
| Binary extremizers & phase | [5. Bernoulli extremizers](../paper/bernoulli-structure/bernoulli-structure.pdf) | Classifies common-diameter maximizers (used in Manuscript 6), proves $C_{2,\mathrm{bin}} = C_1$, and derives the binomial padding pattern. |
| Asymptotic & exact small-$L$ structure | [6. Small-Lyapunov structure](../paper/small-lyapunov-structure/small-lyapunov-structure.pdf) | Proves lattice concentration as $L \to 0$, $\Delta/L \to C_{\mathrm{E}}$, and exact i.i.d. binary structure of global maximizers for $L < \varepsilon_0$. |
| General variational theory | [7. Extremizers](../paper/extremizers/extremizers.pdf) | Proves the four-point and three-point support theorems, indecomposability, and structural restrictions on hypothetical violations of $C_{\mathrm{ind}} = C_{\mathrm{E}}$. |

---

## 3. Quantitative Upper Bounds

### Manuscript 1: Unrestricted Upper Bound ($C_{\mathrm{ind}} < 0.44988$)

[Manuscript 1](../paper/independent-berry-esseen.pdf) establishes $C_{\mathrm{ind}} \le U_* < 0.44988$, improving upon the bound $0.5583$ of [Shevtsova (2013)](https://www.mathnet.ru/eng/ia252).

For a centered random variable $X$ of variance $v > 0$, its zero-biased distribution $X^*$, introduced by [Goldstein and Reinert (1997)](https://doi.org/10.1214/aoap/1043862419), satisfies $\mathbb{E}[X g(X)] = v \mathbb{E} g'(X^*)$ for smooth test functions $g$. Normalizing $V = 1$, the characteristic function $f_S(t) = \mathbb{E} e^{itS}$ satisfies $f_S'(t) = -t f_{S^*}(t)$, which yields the integral representation

```math
f_S(t) - e^{-t^2/2} = e^{-t^2/2} \int_0^t u \bigl(f_S(u) - f_{S^*}(u)\bigr) e^{u^2/2}\,du.
```

Classical proofs bound $|f_S - f_{S^*}|$ using the $L^1$-Wasserstein estimate $W_1(X_j, X_j^*) \le \mathbb{E}|X_j|^3 / (2v_j)$ obtained independently by [Goldstein (2010)](https://doi.org/10.1214/10-AOP527) and [Tyurin (2009)](https://arxiv.org/abs/0912.0726v1). Manuscript 1 introduces three refinements that substantially tighten this Fourier comparison across unequal variances:

1. **Third-moment excess and concave aggregation:** For a standardized summand $W_j = X_j / \sqrt{v_j}$ with $\rho_j = \mathbb{E}|W_j|^3$, the identity $\mathbb{E}[(|W_j| - 1)^2(|W_j| + 1/2)] = \rho_j - 1 \ge 0$ vanishes if and only if $W_j$ is symmetric Rademacher ($\pm 1$ with equal probability), for which $f_{W_j}(x) = \cos x$ and $f_{W_j^*}(x) = \operatorname{sinc} x$. Bounding the real and imaginary parts of $f_{W_j} - f_{W_j^*}$ in terms of $\rho_j - 1$ and taking the least concave majorant aggregates these summand bounds in terms of $L$ and $\tau = \sum_j v_j^{3/2}$.
2. **Signed cosine-product reference:** Following the symmetric Bernoulli comparison studied by [Mattner and Shevtsova (2019)](https://alea.impa.br/articles/v16/16-19.pdf), the difference $f_S(t) - e^{-t^2/2}$ is decomposed through the reference product $C_v(t) = \prod_j \cos(t\sqrt{v_j})$. Because $C_v(t) - e^{-t^2/2}$ is real and nonpositive prior to the first zero of the cosines, its sign can be retained inside Prawitz's smoothing integral, while $|f_S(t) - C_v(t)|$ is controlled by a first-order linear differential inequality whose forcing term vanishes at $\rho_j = 1$.
3. **Isolation of the largest-variance summand:** Tracking the largest summand variance $d = \max_j v_j$ and its third absolute moment $b$ sharpens both the characteristic-function product decay and the Cauchy–Schwarz forcing estimates.

Analytic bounds cover $L \le L_- \approx 0.0014$ (via [Shevtsova (2012)](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf)) and $L \ge L_+ \approx 1.21$ (via Cantelli's inequality). Over the intermediate range $[L_-, L_+]$, a rigorous interval computation subdivides the feasible domain of $(L, d, b, \tau)$ into closed boxes and certifies $\Delta / L \le U_* < 0.44988$.

### Manuscript 2: Two Summands and Bernoulli Arrays

[Manuscript 2](../paper/finite-arrays/finite-arrays.pdf) proves exact and sharp bounds for small arrays and two-point families:

- **Two arbitrary summands:** $C_2 < 49/120 = 0.408333\ldots < C_{\mathrm{E}}$.
- **Common-diameter two-point arrays:** $\Delta \le C_{\mathrm{E}} L$ for any number of independent two-point summands whose nondegenerate members share the same distance $h$ between their two support points, extending [Schulz's binomial theorem (2016)](https://d-nb.info/1197702695/34) to unequal probabilities.
- **At most three two-point summands:** $\Delta \le C_{\mathrm{E}} L$ for up to three independent two-point summands with arbitrary support diameters, and $C_{\mathrm{E}} L - \Delta > 3/80$ for at most two two-point summands.
- **Variance concentration in two summands:** Combining these two-summand bounds with a Stein–Ornstein–Uhlenbeck Gaussian comparison shows that $\Delta \le C_{\mathrm{E}} L$ holds for any independent array in which two two-point summands carry at least $22/23$ of the total variance, or two arbitrary summands carry at least $439/440$ of the total variance.

The proof of $C_2 < 49/120$ fixes the threshold and individual variances and applies strong quadratic moment duality ([Winkler (1988)](https://doi.org/10.1287/moor.13.4.581); [Bertsimas and Popescu (2005)](https://doi.org/10.1137/S1052623401399903)) to replace each summand by an extreme maximizer supported on at most three points. The variational contact equations and second-order variations then reduce a hypothetical maximizer with $C_2 \ge 49/120$ to finitely many combinatorial threshold configurations, each of which is ruled out analytically or by rigorous interval arithmetic.

### Manuscript 3: Sharp Bound for Small Lyapunov Ratio

[Manuscript 3](../paper/small-lyapunov/small-lyapunov.pdf) proves that for every finite independent array,

```math
L \le \exp(-\exp(97)) \quad \Longrightarrow \quad \Delta \le C_{\mathrm{E}} L.
```

This adapts the method introduced by [He and Cheng (2026)](https://arxiv.org/abs/2609.06358v1) for large i.i.d. samples to independent summands with arbitrary, nonidentical distributions. Because adjoining deterministic zero summands increases $n$ without changing $S$ or $L$, the smallness condition for general independent arrays is governed by the Lyapunov ratio $L$ rather than $n$.

The proof proceeds in four stages with fully explicit constants:

1. **Penalized counterexample selection:** If a violation occurred below the cutoff, maximizing $\Delta / L$ minus a smooth penalty on $L$ (localized by Shevtsova's remainder bound $\Delta \le C_{\mathrm{E}} L + 0.2538 L^{4/3}$) would yield an attained maximizer whose individual summands satisfy support-contact equations and uniform support bounds $|X_i| \le 8L\sqrt{V}$.
2. **Uniform jitter and approximate periodicity:** A finite pigeonhole construction selects a smoothing scale $h$ at which the characteristic function of $S$ convolved with $\operatorname{Unif}[-h/2, h/2]$ is close to its one-term Edgeworth expansion.
3. **Support confinement via Esseen's moment inequality:** Quantitative stability in Esseen's lattice moment inequality forces $h$ near the Esseen spacing and confines the support of every summand either to two narrow clusters separated by approximately $h$ or to a small neighborhood of zero carrying negligible total variance.
4. **Local comparison with a common-diameter Bernoulli array:** Writing each two-cluster summand as $h_i(B_i - p_i) + U_i$, the array is compared directly with the common-diameter Bernoulli array $\bar{h}(B_i - p_i)$, which satisfies $\Delta \le C_{\mathrm{E}} L$ by Manuscript 2. Depending on the total variance $\Lambda$ of the diameter deviations and intra-cluster residuals, either a conditional one-sided clipping inequality (for small $\Lambda$) or a conditional Fourier minorant (for larger $\Lambda$) shows that the perturbations strictly reduce $\Delta / L$ below $C_{\mathrm{E}}$.

---

## 4. Extremal Structure for Two-Point Arrays and Small Lyapunov Ratios

### Manuscript 4: Zero-Bias Transport and Repeated Blocks

[Manuscript 4](../paper/transport/transport.pdf) investigates the $L^1$-Wasserstein inequality

```math
2V W_1(S, S^*) \le R,
```

which underlies zero-bias proofs of the Berry–Esseen theorem. Theorem 1.1 proves that equality holds if and only if all nonconstant summands have the *same* centered two-point distribution. The single-summand equality characterization uses Choquet's theorem on the compact convex set of centered probability measures with variance at most one: every extreme point has at most three atoms, the barycenter decomposition automatically forces almost every extreme component to have variance *equal* to one, and strict inequality for three-point distributions reduces the support to two-point distributions. Across multiple summands, unequal support diameters cause strict cancellation under convolution, while for a common diameter an exact conditional-variance identity (Proposition 4.1) expresses the transport deficit $R/(2V) - W_1(S, S^*)$ as a positive weighted sum of variances of the success probabilities $p_i$.

Section 5 combines this transport rigidity with a sharp lattice inequality: for any centered variable $X$ supported on a coset of $h\mathbb{Z}$ with variance $v$ and third moment $\kappa = \mathbb{E} X^3$,

```math
|\kappa| + 3hv \le (3 + \sqrt{10}) \, 2v W_1(X, X^*),
```

with equality if and only if $X$ is a centered binomial variable $h(\operatorname{Bin}(n, p) - np)$ with $p \in \{p_{\mathrm{E}}, 1 - p_{\mathrm{E}}\}$. When a fixed finite block $(X_1, \dots, X_n)$ is repeated $m$ times independently, Esseen's asymptotic formula and these two inequalities yield $\lim_{m \to \infty} C_m \le C_{\mathrm{E}}$, with equality if and only if all nonzero summands in the block are identically distributed Esseen two-point variables.

### Manuscript 5: Bernoulli Extremizers and Binomial Phase

[Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf) studies the Berry–Esseen ratio over independent two-point summands:

- **Common-diameter classification (Theorem 1.1):** For every $n \ge 1$, the maximum of $\Delta / L$ over $n$ independent two-point summands with a common support diameter equals $P_n = \max_{1 \le m \le n} B_m$, and every maximizing array consists of $m$ identically distributed two-point summands and $n - m$ zero summands.
- **Two arbitrary two-point summands (Theorem 1.2):** Over arrays of at most two two-point summands with arbitrary diameters, the maximum ratio is $C_{2,\mathrm{bin}} = C_1 \approx 0.3703515235$, attained only when one summand is zero.
- **Local replacement of unequal diameters (Theorem 1.3):** Fix $c_0 > 0$. For two-point arrays $X_i = g_i(B_i - p_i)$ with $\Delta / L \ge c_0$, total count variance $A = \sum_i p_i(1 - p_i) \ge A_0$, weighted diameter spread $\sigma = \bigl(\sum_i p_i(1 - p_i) e_i^2\bigr)^{1/2} \in (0, \sigma_0]$, and maximum relative diameter deviation $\eta = \max_i |e_i| \le \eta_0$ (where $e_i = g_i/h_0 - 1$ and $h_0$ is the variance-weighted mean diameter), replacing all diameters by $h_0$ while keeping every probability $p_i$ fixed strictly increases $\Delta / L$ by at least $\kappa \sigma^2 / (\sigma + \eta) > 0$, with $A_0, \sigma_0, \eta_0, \kappa > 0$ depending only on $c_0$.
- **Binomial phase and optimal zero summands (Theorem 4.1):** Why can $P_n = \max_{1 \le m \le n} B_m$ be strictly larger for $m < n$ than at $m = n$? A uniform seventh-order Edgeworth–Euler–Maclaurin expansion shows that

```math
B_m = C_{\mathrm{E}} - \frac{Q_*}{m} - \frac{H_* + \alpha \operatorname{dist}(m p_{\mathrm{E}} - \delta, \mathbb{Z})^2}{m^2} + O(m^{-3}).
```

Because the integer threshold $k$ must be an integer while the optimal probability $p$ stays within $O(m^{-1})$ of $p_{\mathrm{E}}$, reducing the active count from $n$ to $n - j$ trades the linear penalty $\eta j$ against the quadratic phase distance $\operatorname{dist}(\theta_n + j p_{\mathrm{E}}, \mathbb{Z})^2$. The lower envelope of these five parabolas proves that for all sufficiently large $n$, the optimal number of zero summands $j_n = n - m_n$ belongs to $\{0, 1, 2, 3, 4\}$ with explicit irrational limiting frequencies, and successive record active sizes $m_n$ have gaps in $\{2, 3, 5\}$.

### Manuscript 6: Structure of Maximizers at Small Lyapunov Ratio

[Manuscript 6](../paper/small-lyapunov-structure/small-lyapunov-structure.pdf) proves two structural theorems for general independent arrays as $L \to 0$:

- **Asymptotic structure near $C_{\mathrm{E}}$ (Theorem 1.1):** For any sequence of independent arrays with $L_n \to 0$ and $\Delta_n / L_n \to C_{\mathrm{E}}$, the summands concentrate in $W_2$ and $W_3$ around centered distributions supported on cosets of a common lattice $h_n \mathbb{Z}$ with $h_n \sim R_n / (r_* V_n)$, the normalized zero-bias transport ratio $2V_n W_1(S_n, S_n^*)/R_n$ converges to $1$, and a subset of $m_n \sim \beta_*^2 / L_n^2$ summands carrying $1 - o(1)$ of the variance and third absolute moments is $W_3$-close to identically distributed Esseen two-point laws.
- **Exact structure of small-$L$ global maximizers (Theorem 1.2):** There exists an absolute (ineffective) threshold $\varepsilon_0 > 0$ such that every global maximizer of $\Delta / L$ over all $n$-summand independent arrays with $L < \varepsilon_0$ consists *exactly* of identically distributed two-point summands and zero summands.

The central obstacle in proving Theorem 1.2 is that asymptotic approximation errors cannot rule out arbitrarily small rounding residuals ("phase noise") modulo the lattice span $h_n$. The proof overcomes this by establishing a *relative* conditional first-moment bound $\mathbb{E}(|T| \mid K = k) \ge c \mathbb{E}|T|$ valid down to arbitrarily small $\mathbb{E}|T| > 0$. This shows that replacing the summands by their exact lattice roundings strictly increases $\Delta / L$ whenever the rounding residual $T$ is nonzero, forcing exact lattice support; discrete Green-function concavity of the support-contact function then reduces each summand to two adjacent lattice points, where Manuscript 5's common-diameter theorem completes the classification.

---

## 5. Variational Theory of General Extremizers

[Manuscript 7](../paper/extremizers/extremizers.pdf) studies the variational properties of general maximizing arrays and of hypothetical counterexamples to $C_{\mathrm{ind}} = C_{\mathrm{E}}$:

1. **Four-point and three-point support reductions (Theorem 3.2):** At any local maximum of $\Delta / L$ with respect to two-sided signed probability perturbations of a summand (preserving mass and mean while allowing variance and third moment to vary), the second logarithmic derivative of $\Delta / L$ restricts the rank of the moment vectors $(1, u, u^2, |u|^3)$ on the support, proving that *every* summand of a maximizer has at most four support points. At a global maximum, any four-point summand lies on a flat probability segment whose endpoint replaces it by a three-point summand with the exact same variance and third absolute moment.
2. **Convolution indecomposability (Lemma 3.3 and Corollary 3.4):** By the strict inequality $\mathbb{E}|U + W|^3 > \mathbb{E}|U|^3 + \mathbb{E}|W|^3$ for nonconstant independent centered variables ([Cox and Kemperman (1983)](https://doi.org/10.1214/aop/1176993520); [Ruzankin (2014)](https://doi.org/10.1016/j.spl.2014.08.006)), splitting a decomposable summand into its independent convolution factors preserves the sum $S$ and total variance $V$ while strictly decreasing $R$. Thus every nonzero summand of an unrestricted maximizer (or of a fixed-$n$ maximizer with at least one zero summand) is convolution-indecomposable.
3. **Compactification and structure above $C_{\mathrm{E}}$ (Theorem 1.2):** If $C_{\mathrm{ind}} > C_{\mathrm{E}}$, Manuscript 3's small-$L$ theorem prevents maximizing sequences from degenerating to $L \to 0$. Passing to the limit via the classical accompanying-Poisson / Lévy–Khintchine compactification yields an attained maximizer consisting of countably many independent three-point summands $X_i$ plus an independent infinitely divisible remainder $Y$ (with Gaussian variance $a$ and Lévy measure $\nu$). Variational contact propagation eliminates the Lévy jump component ($\nu = 0$) and forces $\operatorname{diam}(\operatorname{supp} X_i) \to 0$; a lattice Edgeworth comparison and zero-bias Fourier bound eliminate the Gaussian component ($a = 0$); and a Sperner antichain bound combined with a collective probability variation proves that $S$ is nonatomic whenever infinitely many $X_i$ are nonconstant.
4. **Weighted variance summability and excluded regularity (Theorem 1.3 and Section 10):** Under the additional condition $L \le 18481/25000$, multi-scale grouping combined with Cramér moderate deviations and randomized lattice rounding proves $\sum_i v_i H_i^{-\gamma} < \infty$ for all $\gamma < \sqrt{5} - 1$ (and all $\gamma < \sqrt{21} - 3$ when all but finitely many summands are two-point), ruling out two-sided $C^{1,1}$ regularity of $F_S$ at the maximizing threshold as well as dense or nonoverlapping rare Bernoulli tails.
