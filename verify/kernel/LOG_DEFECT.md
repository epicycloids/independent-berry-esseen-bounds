# Logarithmic bounds for characteristic-function products

This note defines the quantities used by `log_product_cells` and
`log_range_cells`. The derivation is in the product-bound section of
[the paper](../../paper/independent-berry-esseen.pdf).

For centered independent summands, put $v_j=\mathbb E X_j^2$,
$b_j=\mathbb E|X_j|^3$, $V_0=\sum_jv_j$, $B=\sum_jb_j$, and
$\tau=\sum_jv_j^{3/2}$. Let $a$ bound every individual variance, and write
$f_j(t)=\mathbb E e^{itX_j}$. For $t\ge0$ and $\kappa=0.099162$, set

$$
h_j=\frac{v_jt^2}{2}-\kappa(b_j+v_j^{3/2})t^3,
\qquad q_j=\max(h_j,0),
\qquad H=\frac{V_0t^2}{2}-\kappa(B+\tau)t^3.
$$

The single-summand bound gives $|f_j(t)|\le\sqrt{1-2q_j}$, with

$$
0\le q_j\le m(t,a):=
\max_{0\le x\le t\sqrt a}(x^2/2-2\kappa x^3)
\le\frac1{216\kappa^2}<\frac12.
$$

The inequality $-\tfrac12\log(1-2q)\ge q+q^2$ therefore yields

$$
\prod_j|f_j(t)|\le\exp\left(-H-\sum_jq_j^2\right).
$$

For $V_0>0$, put $A_0=t^2/2-2\kappa\sqrt a\,t^3$. If $A_0\ge0$, the
moment-excess argument in the paper gives

$$
\sum_jq_j^2\ge Q:=
\left(\max\left\{0,
 A_0\frac{\tau}{\sqrt{V_0}}-\kappa t^3(B-\tau)
\right\}\right)^2.
$$

Set $Q=0$ when $A_0<0$. After omitting index $i$,

$$
\sum_{j\ne i}q_j\ge\sum_{j\ne i}h_j=H-h_i\ge H-m,
\qquad
\sum_{j\ne i}q_j^2\ge\max(0,Q-m^2).
$$

Thus the full product is bounded by $\min\{1,e^{-H-Q}\}$, and every
product with one factor omitted is bounded by
$\min\{1,e^{-H+m-\max(0,Q-m^2)}\}$. The empty product is one.

`interval_bounds.log_product_cells` encloses these expressions over each
frequency and parameter interval with $0\le\tau\le\texttt{tauhi}$.
`stable_bounds.log_range_cells` uses the supplied interval
$\texttt{taulo}\le\tau\le\texttt{tauhi}$. Both account for the dependence of
$H$ and $Q$ on the same $\tau$ when bounding the exponent. Their argument
`dhi` supplies the variance bound $a$; `deleted=True` selects the bound with
one factor omitted. If the lower endpoint of $V_0$ is nonpositive, the
functions return the upper bound one.

The higher-order logarithmic estimates discussed in the paper were not
used in the recorded calculation. These source functions use the quadratic
correction above.
