# Bounds and normalization

For a finite independent family of centered real random variables, write

$$
V=\sum_j\mathbb E X_j^2>0,\qquad
L=\frac{\sum_j\mathbb E|X_j|^3}{V^{3/2}},\qquad
\Delta=\sup_x\left|\mathbb P\!\left(
\frac{\sum_jX_j}{\sqrt V}\le x\right)-\Phi(x)\right|,
$$

where $\Phi$ is the standard normal distribution function and every third
absolute moment is finite. The constant $C_{\mathrm{ind}}$ is the supremum
of $\Delta/L$ over these families. The companion papers write $\ell$ for $L$.

| Bound on $C_{\mathrm{ind}}$ | Source |
| --- | --- |
| $C_{\mathrm{ind}}\ge C_{\mathrm E}=(\sqrt{10}+3)/(6\sqrt{2\pi})=0.4097321837023963\ldots$ | [Esseen (1956)](https://doi.org/10.1080/03461238.1956.10414946); see also [Shevtsova (2012), Section 1](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf) |
| $C_{\mathrm{ind}}\le0.5583$ | [Shevtsova (2013)](https://www.mathnet.ru/eng/ia252) |
| $C_{\mathrm{ind}}<0.454$ | [Upper-bound paper](../paper/independent-berry-esseen.pdf) |

The conjecture $C_{\mathrm{ind}}=C_{\mathrm E}$ remains open.

## The explicit upper-bound calculation

The interval calculation uses endpoints $L_-\approx0.0014$ and
$L_+\approx1.21$. Their exact rational values are stored in
[the certificate](../certificates/cover.json) and given in the
[paper](../paper/independent-berry-esseen.pdf).

| Range | Estimate used | Bound for $\Delta/L$ |
| --- | --- | --- |
| $0<L\le L_-$ | $C_{\mathrm E}+0.3413L^{1/3}$ | Less than 0.447913037296 |
| $L_-\le L\le L_+$ | Interval calculation of the paper's smoothing bounds | Recorded interval bound $U_*$, defined below |
| $L\ge L_+$ | $\Delta<0.540936541549$ | Less than 0.447055000001 |

The estimate for $0<L\le L_-$ is Corollary 4.18 of
[Shevtsova (2012)](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf),
which applies to independent summands when $L\le0.01$.
The bound on $\Delta$ follows from Cantelli's inequality by maximizing
$\Phi(x)-x^2/(1+x^2)$ over $x>0$.

The calculation covers $[L_-,L_+]$ with 1,683 closed intervals. Its largest
recorded upper bound for $\Delta/L$ is the exact rational number

$$
U_*=\frac{4089268438506339}{9007199254740992}<0.454.
$$

Together, the three ranges give $C_{\mathrm{ind}}\le U_*<0.454$. The
[verification guide](../verify/README.md) describes the arithmetic and
the commands for checking the calculation.

## Classes satisfying the conjectured inequality

The [finite-results paper](../paper/finite-arrays/finite-arrays.pdf) proves
the following statements with $L$ as defined above.

| Class of independent centered summands | Bound |
| --- | --- |
| Variables supported on at most two points, with one common positive support diameter for all nondegenerate summands | $\Delta\le C_{\mathrm E}L$; the constant is optimal over this class |
| Two variables, each supported on at most two points, with arbitrary diameters and positive total variance | $\Delta<C_{\mathrm E}L$ |
| At each coordinate, an arbitrary mixture of a fixed centered Esseen law and its reflection, with the same scale for every coordinate | $\Delta\le C_{\mathrm E}L$ |

In the last row, the variance-one Esseen law gives probability
$p=(4-\sqrt{10})/2$ to $\sqrt{(1-p)/p}$ and probability $1-p$ to
$-\sqrt{p/(1-p)}$. The scale is a common factor $c>0$ multiplying this law
or its reflection; each coordinate may have its own mixture weight.

A common support diameter means that every nondegenerate summand has the
same distance between its two support points. Lattice support alone allows
different integer multiples of a common span and is a broader condition.
The two-summand result does not determine the optimal constant for that
class. The [finite-results supplement](../certificates/finite/README.md)
describes the algebraic and interval computations used in these proofs.

## The sharp bound at small and large $L$

The [small-Lyapunov paper](../paper/small-lyapunov/small-lyapunov.pdf) proves
that a universal $L_0>0$ exists such that every independent array with
$0<L\le L_0$ satisfies $\Delta\le C_{\mathrm E}L$. Its proof adapts
[He and Cheng's i.i.d. argument (2026)](https://arxiv.org/abs/2609.06358v1)
and uses the common-diameter theorem above. The proof gives no numerical
value of $L_0$, so the explicit computation for $C_{\mathrm{ind}}<0.454$
continues to use Shevtsova's estimate at its lower endpoint.

At the other end, the variance-only bound
$\Delta<0.540936541549$ implies $\Delta<C_{\mathrm E}L$ whenever
$L\ge1.321$, since $C_{\mathrm E}>0.4097$ and
$1.321\cdot0.4097=0.5412137$. Thus the unresolved part of the sharp
conjecture lies in $L_0<L<1.321$. The small-$L$ theorem changes neither
Esseen's lower bound nor the numerical upper bound on
$C_{\mathrm{ind}}$.
