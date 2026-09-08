# Real and imaginary parts of the zero-bias error

This note identifies the earlier formula implemented by
`directional_bounds.directional_cells`. The real- and imaginary-part
estimates and their later refinements are derived in
[the paper](../../paper/independent-berry-esseen.pdf).

For a centered variable $Z$ of variance $v>0$, write
$f_Z(x)=\mathbb E e^{ixZ}$. Its zero-bias law $Z^*$ is defined by
$v\mathbb E g'(Z^*)=\mathbb E[Zg(Z)]$ for absolutely continuous test functions
with the required expectations. For a centered variance-one variable $W$, put
$\rho=\mathbb E|W|^3$ and $e=\rho-1\ge0$. Write
$\operatorname{sinc}x=\sin(x)/x$, with value one at zero. The moment excess
$e$ gives the estimates

$$
|\operatorname{Re}f_W(x)-\cos x|\le e|x|^3/6,
\qquad
|\operatorname{Re}f_{W^*}(x)-\operatorname{sinc}x|\le e|x|/2,
$$

and

$$
|\operatorname{Im}(f_W(x)-f_{W^*}(x))|
\le (|x|/2+|x|^3/6)\sqrt{e(\rho+5/3)}.
$$

For a finite collection of centered summands with positive variances, let
$B=\sum_j\mathbb E|X_j|^3$, $\tau=\sum_jv_j^{3/2}$, and
$v_j=\mathbb E X_j^2\le a$. Assume finite third absolute moments;
centered summands of variance zero can be omitted. Standardize each summand
as $W_j=X_j/\sqrt{v_j}$ and put
$\rho_j=\mathbb E|X_j|^3/v_j^{3/2}$.
At $t\ge0$, put $x=t\sqrt a$ and $q(x)=1+x^2/3$.
Scaling the one-variable estimates gives a factor $tv_j^{3/2}/2$ in each
summand. The weighted mean of $\rho_j$ with weights $v_j^{3/2}/\tau$
equals $B/\tau$.
After bounding each $t\sqrt{v_j}$ by $x$, concavity in $\rho_j$ gives

$$
\sum_jv_j|f_{X_j}(t)-f_{X_j^*}(t)|\le tP_{\rm dir}/2,
$$

where

$$
P_{\rm dir}=
\min\left\{B,
\sqrt{\left(\frac{2x\tau}{3}+q(x)(B-\tau)\right)^2
      +q(x)^2(B-\tau)(B+5\tau/3)}\right\}.
$$

For an empty collection, the error sum is zero; set $P_{\rm dir}=0$.

`directional_cells` evaluates an outward upper bound for this expression
using interval endpoints. Its argument `dhi` is the upper variance bound
$a$. This formula uses the two-sided real-part estimate and retains a cross
term when that estimate is squared.

The one-sided estimate described in [EVEN_SUPPORT.md](EVEN_SUPPORT.md)
removes that cross term. The final aggregate formula also uses a polynomial
bound for the cosine/sinc difference. The recorded Gaussian evaluator calls
`paired_bounds.cosine_cells`; the [source guide](README.md) shows the class
inheritance and method overrides. `directional_cells` remains an earlier
implementation, not the final formula used by that evaluator.
