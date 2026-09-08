# Bounds and their scope

All constants below use the classical normalization

$$
L=\frac{\sum_j\mathbb E|X_j|^3}{(\sum_j\mathbb E X_j^2)^{3/2}},
\qquad
\Delta=\sup_x\left|\mathbb P\!\left(
\frac{\sum_jX_j}{\sqrt{\sum_j\mathbb E X_j^2}}\le x\right)-\Phi(x)\right|.
$$

Here $\Phi$ is the standard normal distribution function. The arrays are
finite, independent, and centered, with positive total variance and finite
third absolute moments.

## Global comparison

| Setting | Bound | Status and source |
| --- | --- | --- |
| Unrestricted independent arrays | $C_{\mathrm{ind}}\ge C_{\mathrm E}=0.40973218370239634299\ldots$ | Esseen's lower bound; see [Esseen (1956)](https://doi.org/10.1080/03461238.1956.10414946) and [Shevtsova (2012)](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf) |
| Unrestricted independent arrays | $C_{\mathrm{ind}}\le0.5583$ | Published comparison from [Shevtsova (2013)](https://www.mathnet.ru/eng/ia252) |
| Unrestricted independent arrays | $C_{\mathrm{ind}}<0.474999998$ | This edition's proof candidate; independent mathematical review outstanding |
| Identically distributed independent variables | Upper bound 0.4690 | Published restricted-class comparison from [Shevtsova (2013)](https://www.mathnet.ru/eng/ia252) |
| Binomial sums | Sharp constant $C_{\mathrm E}$ | [Schulz (2016)](https://d-nb.info/1197702695/34) |

Here
$C_{\mathrm E}=(\sqrt{10}+3)/(6\sqrt{2\pi})$.
The table includes selected published comparisons.

The candidate leaves a gap of approximately 0.0652678143 above the lower
bound. The maximum of the finite upper enclosures below is strictly less
than 0.474999998, giving the proposed bound below 0.475 for the supremum.

## The three moment ranges

Let $L_-$ and $L_+$ be the exact binary64 endpoints stored in
[the certificate](../certificates/cover.json). Their printed decimal values
are 0.006 and 1.15. The exact fractions are given in Section 7 of the
[paper](../paper/independent-berry-esseen.pdf).

| Range | Technique | Upper bound for $\Delta/L$ |
| --- | --- | --- |
| $0<L\le L_-$ | Published bound $C_{\mathrm E}+0.3413L^{1/3}$ | Less than 0.471750509536 |
| $L_-\le L\le L_+$ | Completed interval cover using the paper's combined estimates | At most the largest recorded endpoint, 0.47499999772796003 |
| $L\ge L_+$ | Variance-only inequality $\Delta<0.540936541549$ | Less than 0.470379601347 |

The first formula is Corollary 4.18 of
[Shevtsova (2012)](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf),
applicable for $L\le0.01$. The variance-only constant follows by maximizing
$\Phi(x)-x^2/(1+x^2)$ for $x>0$, using Cantelli's inequality. The short
[scalar checker](../verify/scalars.py) encloses both complements with Arb.

The middle range has 545 bands checked over boxes of moment parameters,
215 bands accepted by a uniform bound, and seven by the variance-only
bound.

The [verification guide](../verify/README.md) gives the scope of the original
calculation and subsequent replay. [Remaining questions](questions.md)
discusses additional logarithmic terms and a tighter real-part bound.
