# How the bound is obtained

Let $X_1,\ldots,X_n$ be independent centered real random variables with
total variance one and finite third absolute moments. Write

$$
S=\sum_jX_j,\qquad v_j=\mathbb E X_j^2,\qquad
L=\sum_j\mathbb E|X_j|^3,\qquad
\tau=\sum_jv_j^{3/2},\qquad
\Delta=\sup_x|\mathbb P(S\le x)-\Phi(x)|.
$$

Here $\Phi$ is the standard normal distribution function.

The conjecture asks whether $\Delta/L$ is always at most
$C_{\mathrm E}=(\sqrt{10}+3)/(6\sqrt{2\pi})$. The
[paper](../paper/independent-berry-esseen.pdf) proposes the upper bound
$0.454$, using the following estimates and a finite interval calculation.
Independent mathematical review and a full numerical reevaluation remain
outstanding. The zero-bias and Fourier comparison follows
[Tyurin (2009)](https://arxiv.org/abs/0912.0726v1).

## Moments and characteristic functions

Write $f_X(t)=\mathbb E e^{itX}$. Independence gives
$f_S(t)=\prod_j f_{X_j}(t)$. A moment inequality bounds the modulus of each factor by
$\sqrt{1-2q_j(t)}$, where $0\le q_j(t)<1/2$. The paper bounds
$-\tfrac12\sum_j\log(1-2q_j(t))$ from below using both $L$ and
$\tau$, including the nonnegative third-moment excess
$\sum_j(\mathbb E|X_j|^3-v_j^{3/2})=L-\tau$.
The same argument bounds products with any one factor omitted.

Choose a summand of maximal variance, and denote its variance by $d$
and third absolute moment by $b$. The remaining summands have total
variance $1-d$, third-moment sum $L-b$, and variance-power sum
$\tau-d^{3/2}$. The product estimate and zero-bias comparison use these
parameters to bound the remaining summands separately.

## Zero-bias comparison

The zero-bias transformation was introduced by
[Goldstein and Reinert (1997)](https://doi.org/10.1214/aoap/1043862419).
For a centered variable $X$ of variance $v>0$, its zero-bias law $X^*$ satisfies

$$
v\mathbb E g'(X^*)=\mathbb E[Xg(X)]
$$

for absolutely continuous test functions with finite expectations.
For the normalized sum, $f_S'(t)=-t f_{S^*}(t)$. Hence

$$
f_S(t)-e^{-t^2/2}
=e^{-t^2/2}\int_0^t u\bigl(f_S(u)-f_{S^*}(u)\bigr)e^{u^2/2}\,du.
$$

The difference $f_S-f_{S^*}$ is a sum of single-summand errors multiplied
by products with one factor omitted. The paper bounds its real and
imaginary parts using third-moment excess. A concave majorant combines
these estimates across unequal summands by Jensen's inequality.

Some frequency intervals also use the real reference function
$\prod_j\cos(t\sqrt{v_j})$. Bounds for its distance from $f_S$ retain
part of the sign information that is lost when every Fourier term is
replaced by its modulus. A differential inequality controls that distance;
moment inequalities bound the reference function uniformly over each
parameter box.

## Smoothing and the finite cover

Prawitz's smoothing inequality converts the Fourier estimates into a
bound for the difference between the distribution functions. The calculation
partitions the threshold variable $x$ into closed intervals, bounds the
integrals on whole frequency intervals, and uses Cantelli's inequality
outside a bounded range of $x$. On each threshold interval it takes the
minimum of the bounds obtained from several choices of smoothing parameters.

For each closed interval of $L$, necessary moment inequalities restrict
$(d,b,\tau)$ to a containing box. Subdivision continues until every box
is either excluded by those inequalities or assigned a uniform upper bound
strictly below $0.454$. Where $b$ is subdivided separately, the calculation
retains a cover of its entire feasible interval.

The completed records cover $0.0014\le L\le1.21$, interpreted at the exact
stored binary64 endpoints. Analytic inequalities cover the two complementary
ranges. Their bounds and the exact largest recorded upper bound appear in the
[bounds note](bounds.md). The [verification guide](../verify/README.md)
describes how the supplied evidence can be checked.
