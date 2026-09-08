# The remaining gap

The conjectured constant is
$C_{\mathrm E}=(\sqrt{10}+3)/(6\sqrt{2\pi})$.
The proof candidate below 0.475 awaits independent mathematical review.
The refinements below have not produced a smaller completed unrestricted
bound in this edition.

## Retain more terms in the product logarithm

For $t\ge0$, the [logarithmic-bounds note](../verify/kernel/LOG_DEFECT.md)
gives factor bounds $|f_{X_j}(t)|\le\sqrt{1-2q_j}$. The paper retains
the quadratic term in the expansion

$$
-\tfrac12\log(1-2q)
=q+\sum_{p=2}^{\infty}\frac{2^{p-1}}p q^p,
\qquad 0\le q<\tfrac12.
$$

Every omitted term is nonnegative. Section 8 of
[the paper](../paper/independent-berry-esseen.pdf) uses a norm inequality
and Hölder interpolation to derive lower bounds $Q_p$ for
$\sum_j q_j^p$. Keeping finitely many of these terms also sharpens bounds
for products with one factor omitted, after subtracting an upper bound
for that factor's contribution to each power sum.

A stronger numerical bound would require interval implementations with
outward rounding and a complete cover at the new target.

## Use both endpoints of the real-part interval

For a centered variance-one variable $W$ and its zero-bias law $W^*$,
write $f_W$ and $f_{W^*}$ for their characteristic functions. The paper proves

$$
-g\le\operatorname{Re}(f_W(x)-f_{W^*}(x))\le-g+ce,
\qquad g=\operatorname{sinc}(x)-\cos x,\quad
c=x/2+x^3/6,\quad e=\mathbb E|W|^3-1
$$

when $0\le x\le2$, with $\operatorname{sinc}(x)=\sin(x)/x$ and its
continuous value one at zero. The current aggregation uses
$g^2+c^2e^2$ to bound the square of this real part. The sharper
single-summand expression is $\max\{g^2,(ce-g)^2\}$.

To aggregate this expression across unequal summands, one needs a
concave majorant or additional moment information that separates its
two branches.

## Preserve relative phase

The current zero-bias comparison takes the modulus of each complex summand
error before summing, discarding cancellation between those errors.
Retaining two summands separately, or controlling an aggregate
phase, could reduce that loss. Such a result would need to work uniformly
over both lattice and nonlattice arrays and over unequal variances.

## Sharpen the moment domain or the numerical enclosure

After normalizing the total variance to one, the certificate restricts
the largest summand variance $d$, its third absolute moment $b$, and
$\tau=\sum_j(\mathbb E X_j^2)^{3/2}$ by necessary moment inequalities.
Some allowed tuples may not be realizable by any array. Further necessary
inequalities could exclude these tuples.

Finer subdivision of the parameter boxes or frequency intervals can reduce
overestimation introduced by interval evaluation. Trying more smoothing
parameters may lower the bound obtained from the same smoothing inequality.
