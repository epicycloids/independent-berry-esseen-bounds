# A one-sided bound for the real zero-bias error

This note connects the quadratic inequality in
[the paper](../../paper/independent-berry-esseen.pdf) to `even_cells` and its
later refinement `cosine_cells`. The characteristic-function and zero-bias
notation is defined in [DIRECTIONAL_STABILITY.md](DIRECTIONAL_STABILITY.md).

Let $W$ be centered with variance one, let $W^*$ have its zero-bias law, and
put $e=\mathbb E|W|^3-1$. For $0<x\le2$, define

$$
F_x(y)=\cos(xy)-y\sin(xy)/x,
\qquad A_x=\cos x+(1+x^2)\operatorname{sinc}x.
$$

The quadratic inequality proved in the paper is

$$
F_x(y)\ge F_x(1)-\frac{A_x}{2}(y^2-1),\qquad y\ge0.
$$

The right side is tangent at $y=1$ when written as a function of $y^2$.
Put $Y=|W|$. Applying the paper's zero-bias identity with
$h(w)=\sin(xw)/x$ gives
$\mathbb E\cos(xW^*)=\mathbb E[W\sin(xW)]/x$, and hence

$$
\mathbb E F_x(Y)
=\operatorname{Re}(f_W(x)-f_{W^*}(x)).
$$

Since $\mathbb E Y^2=1$, the quadratic term in the minorant has expectation
zero. Define $g=\operatorname{sinc}x-\cos x\ge0$ and $c=x/2+x^3/6$.
The minorant and the two-sided moment-excess estimate give

$$
-g\le\operatorname{Re}(f_W(x)-f_{W^*}(x))\le-g+ce,
\qquad
|\operatorname{Re}(f_W(x)-f_{W^*}(x))|^2\le g^2+c^2e^2.
$$

The maximum squared value on the displayed real-part interval is
$\max\{g^2,(ce-g)^2\}$. Because $g,c,e\ge0$, this is at most
$g^2+c^2e^2$. Combining this last bound with the imaginary-part estimate
and applying concavity gives the aggregate formulas below.

## The expression called `E3`

For summand variances $v_j\le a$, let $B=\sum_j\mathbb E|X_j|^3$,
$\tau=\sum_jv_j^{3/2}$, $x=t\sqrt a$, and $q(x)=1+x^2/3$, with $t\ge0$.
The diagnostic label `EVEN_SUPPORT.md E3` in `even_bounds.probe` refers to

$$
P_{\rm even}(t;B,a,\tau)=
\min\left\{B,
\sqrt{(2x\tau/3)^2+q(x)^2(B-\tau)(2B+2\tau/3)}\right\}.
$$

This is the expression enclosed by `even_bounds.even_cells`. It uses
$\alpha(x)=2x/3$. The resulting bound on the sum of individual zero-bias
errors is

$$
\sum_jv_j|f_{X_j}(t)-f_{X_j^*}(t)|
\le\frac t2 P_{\rm even}(t;B,a,\tau).
$$

For $0<x\le2$, this uses
$2(\operatorname{sinc}x-\cos x)/x\le2x/3$.
At zero the bound follows by continuity. For $x>2$, the square-root term
is at least $B$, so the minimum uses the classical bound $tB/2$.

## The final aggregate formula

On $0\le x\le2$, the paper sharpens the upper bound for
$2(\operatorname{sinc}x-\cos x)/x$ to the following polynomial, with the
quotient defined by continuity at zero. The definition for $x>2$ retains
the classical branch of the minimum:

$$
\alpha(x)=
\begin{cases}
2x/3-x^3/15+x^5/420,&0\le x\le2,\\
2x/3,&x>2,
\end{cases}
$$

The resulting aggregate formula is

$$
P(t;B,a,\tau)=
\min\left\{B,
\sqrt{\alpha(x)^2\tau^2+q(x)^2(B-\tau)(2B+2\tau/3)}\right\}.
$$

`paired_bounds.alpha_upper` encloses the piecewise polynomial, and
`paired_bounds.cosine_cells` uses it to enclose $P$. For $x>2$, the minimum
selects the classical bound $B$. Thus `E3` is the earlier expression with
$\alpha(x)=2x/3$, rather than the final paper formula.

The recorded evaluator obtains its aggregate estimate from
`CosineWeights.prefactor`, which calls `cosine_cells`. For the separate task
of bounding a single characteristic-function modulus,
`PairedWeights.single_cells` takes the minimum of the original bound and
`paired_bounds.paired_cells`. The [source guide](README.md) identifies the
relevant method overrides.
