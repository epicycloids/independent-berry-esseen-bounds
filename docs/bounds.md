# Summary of Bounds and Regimes

Navigation: [Repository Root](../README.md) · [Overview & Reading Guide](overview.md) · [Numerical Studies & Conjectures](numerical-observations.md) · [Verification Guide](../verify/README.md)

---

## Normalization and Definitions

For a finite family of independent centered real random variables $`X_1, \dots, X_n`$ with positive total variance and finite third absolute moments, we write

```math
\begin{aligned}
V &= \sum_{i=1}^n \mathbb{E} X_i^2 > 0, \\
R &= \sum_{i=1}^n \mathbb{E}|X_i|^3 < \infty, \\
L &= \frac{R}{V^{3/2}}, \\
\Delta &= \sup_{x \in \mathbb{R}} \left|\mathbb{P}\!\left(\frac{\sum_{i=1}^n X_i}{\sqrt{V}} \le x\right) - \Phi(x)\right|,
\end{aligned}
```

where $`\Phi`$ is the standard normal distribution function. Here $`L`$ is the Lyapunov ratio and $`\Delta / L`$ is the Berry–Esseen ratio. We denote by $`C_{\mathrm{ind}} = \sup \Delta / L`$ the universal constant over all finite independent arrays, and by $`C_n`$ the supremum over arrays of at most $`n`$ independent summands.

A *two-point* (or *centered Bernoulli*) summand is a centered random variable supported on at most two points, which takes the form $`X_i = h_i(B_i - p_i)`$ with $`B_i \sim \operatorname{Bern}(p_i)`$, $`p_i \in (0, 1)`$, and support diameter $`h_i > 0`$ (or $`X_i = 0`$ in the degenerate case).

## 1. Global Bounds on $`C_{\mathrm{ind}}`$

| Bound | Value | Reference |
| --- | --- | --- |
| Lower bound | $`C_{\mathrm{ind}} \ge C_{\mathrm{E}} = \dfrac{3 + \sqrt{10}}{6\sqrt{2\pi}} = 0.40973218370239634299\ldots`$ | [Esseen (1956)](https://doi.org/10.1080/03461238.1956.10414946) |
| Previous published upper bound | $`C_{\mathrm{ind}} \le 0.5583`$ | [Shevtsova (2013)](https://www.mathnet.ru/eng/ia252) |
| Upper bound in this collection | $`C_{\mathrm{ind}} \le U_* < 0.44988`$ | [Manuscript 1](../paper/independent-berry-esseen.pdf) |

The unconditional upper bound $`C_{\mathrm{ind}} \le U_*`$ is established by partitioning the positive half-line of Lyapunov ratios $`L > 0`$ into three regimes at the exact dyadic rational endpoints $`L_- = 6456360425798343 / 2^{52} \approx 0.0014`$ and $`L_+ = 1362338887279575 / 2^{50} \approx 1.21`$:

| Lyapunov ratio regime | Method | Upper bound on $`\Delta / L`$ |
| --- | --- | --- |
| $`0 < L \le L_-`$ | Asymptotic bound $`\Delta / L \le C_{\mathrm{E}} + 0.3413 L^{1/3}`$ ([Shevtsova (2012)](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf), Corollary 4.18) | $`< 0.447913037296`$ |
| $`L_- \le L \le L_+`$ | Certified interval cover of the moment domain $`(L, d, b, \tau)`$ using signed cosine-reference Prawitz smoothing | $`\le U_* = \dfrac{1013039690742047}{2251799813685248} = 0.449879995808\ldots`$ |
| $`L \ge L_+`$ | Variance-only Cantelli bound $`\Delta \le \sup_{x > 0}\bigl(\Phi(x) - \frac{x^2}{1+x^2}\bigr) < 0.540936541549`$ | $`< 0.447055000001`$ |

See the [Verification Guide](../verify/README.md) for the structure of the 1,670-band moment partition in [`certificates/cover.json`](../certificates/cover.json).

---

## 2. Regimes Where the Sharp Inequality $`\Delta \le C_{\mathrm{E}} L`$ Holds

Although the unrestricted conjecture $`C_{\mathrm{ind}} = C_{\mathrm{E}}`$ remains open on the intermediate interval $`\exp(-\exp(97)) < L < 1.321`$, the sharp inequality $`\Delta \le C_{\mathrm{E}} L`$ (or a strictly smaller bound) is established in the following regimes:

| Class of independent arrays | Bound established | Source |
| --- | --- | --- |
| **Small Lyapunov ratio:** arbitrary independent summands with $`0 < L \le \ell_0 := \exp(-\exp(97))`$ | $`\Delta \le C_{\mathrm{E}} L`$ (optimal) | [Manuscript 3](../paper/small-lyapunov/small-lyapunov.pdf), Theorem 1.1 |
| **Large Lyapunov ratio:** arbitrary independent summands with $`L \ge 1.321`$ | $`\Delta < 0.54093655 < C_{\mathrm{E}} L`$ | Cantelli's inequality ($`0.54093655 / 1.321 < C_{\mathrm{E}}`$) |
| **At most two arbitrary summands ($`n \le 2`$):** arbitrary centered laws $`X_1, X_2`$ | $`C_2 < \dfrac{49}{120} = 0.408333\ldots < C_{\mathrm{E}}`$ | [Manuscript 2](../paper/finite-arrays/finite-arrays.pdf), Theorem 1.1 |
| **At most two two-point summands:** $`X_1, X_2`$ supported on at most two points each | $`C_{\mathrm{E}} L - \Delta > \dfrac{3}{80}`$ and $`\dfrac{\Delta}{L} \le C_1 \approx 0.370352`$ | [Manuscript 2](../paper/finite-arrays/finite-arrays.pdf), Theorem 1.3; [Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf), Theorem 1.2 |
| **At most three two-point summands ($`n \le 3`$):** arbitrary support diameters $`h_1, h_2, h_3 > 0`$ | $`\Delta \le C_{\mathrm{E}} L`$ | [Manuscript 2](../paper/finite-arrays/finite-arrays.pdf), Theorem 1.2(i) |
| **Common-diameter two-point arrays:** any $`n \ge 1`$ two-point summands with equal support diameters $`h_i = h`$ | $`\Delta \le C_{\mathrm{E}} L`$ (optimal) | [Manuscript 2](../paper/finite-arrays/finite-arrays.pdf), Theorem 1.2(ii) |
| **Two dominant two-point summands:** arbitrary array containing two two-point summands with $`v_1 + v_2 \ge \frac{22}{23} V`$ | $`\Delta \le C_{\mathrm{E}} L`$ | [Manuscript 2](../paper/finite-arrays/finite-arrays.pdf), Corollary 1.4 |
| **Two dominant arbitrary summands:** arbitrary array containing two summands with $`v_1 + v_2 \ge \frac{439}{440} V`$ | $`\Delta \le C_{\mathrm{E}} L`$ | [Manuscript 2](../paper/finite-arrays/finite-arrays.pdf), Corollary 1.4 |
| **Two-block moment condition:** two summands (or grouped blocks with third-moment increase $`e \le R/500`$) having $`v_1 + v_2 \ge \frac{11}{16} V`$ and $`b_1 + b_2 \le \frac{1}{2} R`$ | $`\Delta < C_{\mathrm{E}} L`$ | [Manuscript 2](../paper/finite-arrays/finite-arrays.pdf), Section 7 |

---

## 3. Exact Extremizers and Asymptotic Limits

Beyond upper bounds on $`\Delta / L`$, Manuscripts 4–7 establish exact formulas and structural classifications for extremal arrays:

1. **One summand and two binary summands:** By [Bentkus and Kirša (1989)](https://doi.org/10.1007/BF00966194), the one-summand constant is achieved by a two-point distribution:

   ```math
   C_1 = \max_{0 < p < 1} \frac{\sqrt{p(1 - p)}}{p^2 + (1 - p)^2} \left\{\Phi\!\left(\sqrt{\frac{p}{1 - p}}\right) - p\right\} = 0.370351523474378920\ldots.
   ```

   [Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf) (Theorem 1.2) proves that the two-summand two-point supremum $`C_{2,\mathrm{bin}}`$ equals $`C_1`$, and that every attaining two-summand two-point array has one zero summand.
2. **Common-diameter two-point maximizers:** [Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf) (Theorem 1.1) proves that for every $`n \ge 1`$, the maximum $`C_n^{\mathrm{com}}`$ of $`\Delta / L`$ over at most $`n`$ two-point summands with a common support diameter satisfies

   ```math
   C_n^{\mathrm{com}} = P_n := \max_{1 \le m \le n} B_m,
   ```

   where $`B_m`$ is the $`m`$-summand i.i.d. binomial optimum, and every maximizing array consists of $`m`$ identically distributed two-point summands and $`n - m`$ zero summands.
3. **Global maximizers at small Lyapunov ratio:** [Manuscript 6](../paper/small-lyapunov-structure/small-lyapunov-structure.pdf) (Theorem 1.2) proves that there exists an absolute constant $`\varepsilon_0 > 0`$ (not explicitly computed) such that any array attaining the unrestricted $`n`$-summand supremum $`C_n`$ with $`L < \varepsilon_0`$ consists of identically distributed two-point summands and zero summands.
4. **Repeated finite blocks:** Given a fixed finite independent block $`(X_1, \dots, X_n)`$ with variance $`V`$, third absolute moment sum $`R`$, third signed moment sum $`\kappa = \sum_{i=1}^n \mathbb{E} X_i^3`$, and maximal sum-lattice span $`h \ge 0`$ ($`h = 0`$ if nonlattice), [Manuscript 4](../paper/transport/transport.pdf) (Theorem 1.2) proves that for $`m`$ independent copies of the block,

   ```math
   \lim_{m \to \infty} \frac{\Delta_m}{L_m} = \frac{|\kappa| + 3hV}{6\sqrt{2\pi}\,R} \le C_{\mathrm{E}} \frac{2V W_1(S, S^*)}{R} \le C_{\mathrm{E}},
   ```

   with equality at $`C_{\mathrm{E}}`$ if and only if all nonzero $`X_i`$ are identically distributed Esseen two-point variables.
