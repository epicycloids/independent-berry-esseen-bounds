# Bounds and normalization

For a finite independent family of centered real random variables, write

```math
V=\sum_j\mathbb E X_j^2>0,\qquad
L=\frac{\sum_j\mathbb E|X_j|^3}{V^{3/2}},\qquad
\Delta=\sup_x\left|\mathbb P\!\left(
\frac{\sum_jX_j}{\sqrt V}\le x\right)-\Phi(x)\right|.
```

Here $\Phi$ is the standard normal distribution function, and every third
absolute moment is finite. The constant $C_{\mathrm{ind}}$
is the supremum of $\Delta/L$ over these families; $C_2$ is the same
supremum restricted to families of at most two summands. The companion
papers also use $\ell$ for $L$.

| Bound on $C_{\mathrm{ind}}$ | Source |
| --- | --- |
| $C_{\mathrm{ind}}\ge C_{\mathrm E}=(3+\sqrt{10})/(6\sqrt{2\pi})=0.409732\ldots$ | [Esseen (1956)](https://doi.org/10.1080/03461238.1956.10414946) |
| $C_{\mathrm{ind}}\lt 0.44988$ | [Upper-bound paper](../paper/independent-berry-esseen.pdf) |

The conjecture $C_{\mathrm{ind}}=C_{\mathrm E}$ remains open.

## The general upper bound

The computation covers an intermediate interval $[L_-,L_+]$, where
$L_-\approx0.0014$ and $L_+\approx1.21$. The exact rational endpoints
are specified in the [paper](../paper/independent-berry-esseen.pdf) and
[certificate](../certificates/cover.json).

| Range | Estimate | Upper bound for $\Delta/L$ |
| --- | --- | --- |
| $0\lt L\le L_-$ | $C_{\mathrm E}+0.3413L^{1/3}$ | $0.447913037296$ |
| $L_-\le L\le L_+$ | Interval evaluation of the smoothing bounds | $U_*$ below |
| $L\ge L_+$ | $\Delta\lt 0.540936541549$ | $0.447055000001$ |

The first estimate is Corollary 4.18 of
[Shevtsova (2012)](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf),
valid for independent summands when $L\le0.01$. The variance-only bound
in the last row follows from Cantelli's inequality by maximizing
$\Phi(x)-x^2/(1+x^2)$ for $x\gt 0$.

The largest bound in the intermediate interval is

```math
U_*=\frac{1013039690742047}{2251799813685248}
=0.449879995808387445\ldots<0.44988.
```

Thus $C_{\mathrm{ind}}\le U_*$. The
[verification guide](../verify/README.md) describes the finite partition,
the arithmetic, and the checks supplied with the calculation.

## Finite classes satisfying the sharp inequality

The [finite-array paper](../paper/finite-arrays/finite-arrays.pdf) proves:

| Class | Bound |
| --- | --- |
| At most two arbitrary summands | $C_2\lt 49/120=0.408333\ldots\lt C_{\mathrm E}$ |
| At most three summands, each supported on at most two points, with arbitrary support diameters | $\Delta\le C_{\mathrm E}L$ |
| Any number of summands, each supported on at most two points, with a common positive support diameter | $\Delta\le C_{\mathrm E}L$; the constant is optimal over this class |
| Arrays with two Bernoulli coordinates carrying at least $22/23$ of total variance | $\Delta\le C_{\mathrm E}L$ |
| Arrays with two arbitrary coordinates carrying at least $439/440$ of total variance | $\Delta\le C_{\mathrm E}L$ |

Here a Bernoulli coordinate means a centered variable supported on at most
two points. A common support diameter is the same distance between the
two support points of every nondegenerate summand.

The paper also gives conditions involving the third absolute moments of
a selected pair. For example, after normalizing total variance to one,
two coordinates satisfying

```math
\operatorname{Var}(X_i)+\operatorname{Var}(X_j)\ge\frac{11}{16},
\qquad \mathbb E|X_i|^3+\mathbb E|X_j|^3\le\frac L2
```

imply the sharp inequality for the full array. If coordinates are grouped
into new summands, the increase in their third absolute moments must also
be included.

The Bernoulli variance threshold follows from the uniform two-coordinate
margin $C_{\mathrm E}L-\Delta\gt 3/80$ at total variance one. The
[finite-results supplement](../certificates/finite/README.md) describes
the interval calculations supporting the finite theorems.

## Small and large Lyapunov ratios

The [small-Lyapunov paper](../paper/small-lyapunov/small-lyapunov.pdf) proves

```math
0<L\le\ell_0:=\exp(-\exp(97))
\quad\Longrightarrow\quad\Delta\le C_{\mathrm E}L
```

for every independent array. The cutoff follows from explicit estimates in the proof. Esseen's Bernoulli arrays
have $L\to0$ and $\Delta/L\to C_{\mathrm E}$, so the constant is optimal
even in this restricted range.

The argument adapts [He and Cheng's method (2026)](https://arxiv.org/abs/2609.06358v1)
and uses the common-diameter theorem. Its quantitative selection also
uses Shevtsova's Corollary 4.18, specifically
$\Delta\le C_{\mathrm E}L+0.2538L^{4/3}$ for $L\le10^{-3}$.

The variance-only bound $\Delta\lt 0.540936541549$ gives
$\Delta\lt C_{\mathrm E}L$ for $L\ge1.321$: one has $C_{\mathrm E}\gt 0.4097$
and $1.321\cdot0.4097=0.5412137$. The unresolved part of the sharp
conjecture therefore lies in $\ell_0\lt L\lt 1.321$, outside the classes listed
above. The [extremizer paper](../paper/extremizers/extremizers.pdf) gives
further structural restrictions on possible maximizing violations.

## Bernoulli maximizers and transport

The [Bernoulli paper](../paper/bernoulli-structure/bernoulli-structure.pdf)
proves that every fixed-size common-diameter maximizer has identically
distributed active coordinates, with any remaining coordinates zero. The
optimal number of active coordinates can depend on the size of the array.
For at most two Bernoulli coordinates with arbitrary diameters, the exact
optimum is

```math
C_{2,\mathrm{binary}}=C_1=0.370351523474378920\ldots.
```

Every attaining array has one active coordinate. If both variances are
required to be positive, the same supremum is approached but not attained.
The equality $C_2=C_1$ for two arbitrary summands is not established.

The paper's diameter-replacement theorem assumes a sufficiently large
Bernoulli variance sum and small relative diameter variation. Its binomial
phase analysis explains the asymptotic choice of active count when zero
padding is allowed. The [numerical observations](numerical-observations.md)
give the stored finite-horizon candidates and the coverage of their searches.
The conjecture that these candidate values equal the unrestricted $C_n$
remains open.

The [general structure paper](../paper/small-lyapunov-structure/small-lyapunov-structure.pdf)
proves an absolute, ineffective $\varepsilon_0\gt0$ such that every full
finite-horizon global maximizer with $L\lt\varepsilon_0$ has identical
nonzero two-point coordinates and zeros. This structural threshold is
separate from the explicit value cutoff $\exp(-\exp(97))$. For arbitrary
changing arrays with $L\to0$ and $\Delta/L\to C_{\mathrm E}$, the same
paper proves common-lattice approximation and asymptotic equality in
additive zero-bias transport; exact binary structure requires globality.

For a fixed finite block of summands, let $V$ be its total variance,
$R$ its third absolute moment sum, $\kappa$ its third signed moment sum,
and $h$ the maximal lattice span of its sum, with $h=0$ in the nonlattice
case. The [transport paper](../paper/transport/transport.pdf) proves that
repeating the block $M$ times gives

```math
\lim_{M\to\infty}\frac{\Delta_M}{L_M}
=\frac{|\kappa|+3hV}{6\sqrt{2\pi}\,R}\le C_{\mathrm E}.
```

Equality holds precisely for identical nonconstant Esseen Bernoulli
coordinates, together with possible zero coordinates. This limit requires
the block to remain fixed as $M$ increases.
