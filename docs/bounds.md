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
of $\Delta/L$ over these families.

| Bound on $C_{\mathrm{ind}}$ | Status and source |
| --- | --- |
| $C_{\mathrm{ind}}\ge C_{\mathrm E}=(\sqrt{10}+3)/(6\sqrt{2\pi})=0.4097321837023963\ldots$ | [Esseen (1956)](https://doi.org/10.1080/03461238.1956.10414946); see also [Shevtsova (2012), Section 1](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf) |
| $C_{\mathrm{ind}}\le0.5583$ | Published upper bound, [Shevtsova (2013)](https://www.mathnet.ru/eng/ia252) |
| $C_{\mathrm{ind}}<0.454$ | Upper bound proposed in this draft |

Independent mathematical review remains outstanding, and the interval
calculation has not been recomputed in full. The conjecture
$C_{\mathrm{ind}}=C_{\mathrm E}$ remains open.

## Three ranges of $L$

The interval calculation uses endpoints $L_-\approx0.0014$ and
$L_+\approx1.21$. Their exact rational values are stored in
[the certificate](../certificates/cover.json) and given in the
[paper](../paper/independent-berry-esseen.pdf).

| Range | Estimate used | Upper bound for $\Delta/L$ |
| --- | --- | --- |
| $0<L\le L_-$ | $C_{\mathrm E}+0.3413L^{1/3}$ | Less than 0.447913037296 |
| $L_-\le L\le L_+$ | Interval calculation of the paper's smoothing bounds | At most $U_*$, defined below |
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

The bounds in the table hold uniformly on their stated ranges. Their
maximum gives the proposed upper bound for $C_{\mathrm{ind}}$. See the
[verification guide](../verify/README.md) for the completed checks and
arithmetic assumptions.
