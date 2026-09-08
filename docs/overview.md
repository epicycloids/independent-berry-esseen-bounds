# A guide to the argument

Let a finite array of independent centered real variables have total
variance one and finite third absolute moments. Write

$$
S=\sum_jX_j,\qquad L=\sum_j\mathbb E|X_j|^3,\qquad
\Delta=\sup_x|\mathbb P(S\le x)-\Phi(x)|.
$$

Here $\Phi$ is the standard normal distribution function.
The constant under study is the supremum of the ratio $\Delta/L$ over all
finite such arrays. Individual variables may have different laws and
different variances. The conjecture is that the supremum equals Esseen's
constant, approximately 0.40973218. This edition proposes an upper bound
strictly below 0.475. Independent mathematical review remains outstanding.

This guide follows the argument in
[the paper](../paper/independent-berry-esseen.pdf).

## Characteristic functions and zero-bias comparison

For a random variable $X$, write $f_X(t)=\mathbb E e^{itX}$ for its
characteristic function. Independence gives
$f_S(t)=\prod_j f_{X_j}(t)$. For $t\ge0$, Section 2 constructs bounds
of the form
$|f_{X_j}(t)|\le\sqrt{1-2q_j(t)}$, with $0\le q_j(t)<1/2$.
The inequality $-\tfrac12\log(1-2q)\ge q+q^2$ then gives

$$
|f_S(t)|\le\exp\left(-\sum_jq_j(t)-\sum_jq_j(t)^2\right).
$$

Moment inequalities bound the sums in the exponent from below. The same
estimates bound products with one factor omitted.

For a centered variable $X$ of variance $v>0$, its zero-bias law $X^*$ is
defined by

$$
v\mathbb E g'(X^*)=\mathbb E[Xg(X)]
$$

for absolutely continuous test functions for which the expectations exist.
The proof bounds $f_S-f_{S^*}$ using the individual differences
$f_{X_j}-f_{X_j^*}$ and the products with one factor omitted.
Centered summands of variance zero can be omitted.

## Third-moment excess

For a centered variance-one variable $W$ with finite third absolute moment,
set $Y=|W|$ and $\rho=\mathbb E|W|^3$. Then

$$
\mathbb E[(Y-1)^2(Y+1/2)]=\rho-1\ge0.
$$

Equality holds precisely when $W$ takes the values $-1$ and $1$ with
equal probability. Its characteristic function is then $\cos x$, and
that of its zero-bias law is $\operatorname{sinc}x=\sin(x)/x$, with
value one at zero.

Sections 3–4 bound the differences of $\operatorname{Re}f_W$ from cosine
and $\operatorname{Re}f_{W^*}$ from sinc in terms of $\rho-1$. They obtain
a one-sided bound for $\operatorname{Re}(f_W-f_{W^*})$ and a separate
bound for $\operatorname{Im}(f_W-f_{W^*})$. Concavity combines the summand
estimates with weights proportional to $(\mathbb E X_j^2)^{3/2}$.

## A maximal-variance summand

Let $d$ be the largest summand variance, $b$ its third absolute moment,
and $\tau=\sum_j(\mathbb E X_j^2)^{3/2}$. Section 5 bounds that summand
and the sum of the remaining summands separately, both in the
characteristic-function product and in the zero-bias error. Necessary
constraints on $(d,b,\tau)$ restrict the parameter domain while including
every array under consideration.

Since the sum has variance one, the zero-bias identity gives
$f_S'(t)=-t f_{S^*}(t)$. Comparing this with the derivative of the Gaussian
characteristic function gives, for $t\ge0$,

$$
f_S(t)-e^{-t^2/2}
=e^{-t^2/2}\int_0^t u\bigl(f_S(u)-f_{S^*}(u)\bigr)e^{u^2/2}\,du.
$$

## Smoothing

Section 6 applies Prawitz's smoothing inequality to turn the
characteristic-function bounds into a bound on $\Delta/L$. The integrals
are split into frequency intervals and bounded over each interval.

## Interval calculation

For each closed interval of $L$, the calculation in Section 7 subdivides
a box of possible $(d,b,\tau)$ values. A box is accepted only when its
computed upper bound for $\Delta/L$ is finite and below 0.475 throughout
the box. It is discarded when the moment constraints prove it infeasible;
otherwise, it is split into two closed boxes. The complete tree is
retained. After subdivision, the recorded maxima are checked against the
finer candidate bound 0.474999998. Published small-moment estimates and a
variance-only inequality cover the two complementary ranges.

The [bounds note](bounds.md) gives the resulting numbers. The
[verification guide](../verify/README.md) distinguishes checking saved
records, reconstructing the partition, and recomputing numerical leaf bounds.
