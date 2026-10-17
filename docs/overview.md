# Reading the papers

The problem is to determine the least constant $C_{\mathrm{ind}}$ in the
Berry–Esseen inequality for independent real summands. The
[bounds note](bounds.md) defines the normalized third-moment sum $L$ and
normal-approximation error $\Delta$. The conjectured constant is Esseen's
$C_{\mathrm E}=(3+\sqrt{10})/(6\sqrt{2\pi})$.

The seven papers have different purposes. For the sharp conjecture, begin with
the finite-array and small-Lyapunov papers. The transport and Bernoulli
papers explain when two-point laws and common spacings give equality, and
what additional hypotheses permit replacing unequal spacings. The
extremizer paper studies full maximizers and a hypothetical violation of
the conjecture. The general structure paper treats changing arrays with
small Lyapunov ratio. The unrestricted upper-bound computation can be
read independently. The [numerical observations](numerical-observations.md)
explain the candidate family, the observed zero padding, and the open
reduction from arbitrary independent laws to that family.

## The unrestricted upper bound

The [upper-bound paper](../paper/independent-berry-esseen.pdf) proves
$C_{\mathrm{ind}}\lt 0.44988$. It combines characteristic-function estimates
and Prawitz smoothing with the zero-bias comparison developed by
[Tyurin (2009)](https://arxiv.org/abs/0912.0726v1). The zero-bias
transformation was introduced by
[Goldstein and Reinert (1997)](https://doi.org/10.1214/aoap/1043862419);
the sharp mean-distance estimate was obtained independently by
[Goldstein (2010)](https://doi.org/10.1214/10-AOP527) and Tyurin.
Comparison with symmetric two-point variables of matching variances uses
[Mattner and Shevtsova (2019)](https://alea.impa.br/articles/v16/16-19.pdf).
The reference characteristic function is real, so the smoothing integral
can retain the sign of its difference from the Gaussian characteristic
function.

Analytic estimates cover small and large $L$. An interval computation covers
the intermediate range by subdividing the moment parameters and bounding
the smoothing integral on every resulting box. The
[verification guide](../verify/README.md) distinguishes inspection of the
stored certificate from reevaluation of its numerical inputs.

## Finite arrays

The [finite-array paper](../paper/finite-arrays/finite-arrays.pdf) proves
$C_2\lt 49/120$ for two arbitrary summands. It proves the $C_{\mathrm E}$ bound
for three Bernoulli summands of arbitrary diameters and for any number
with a common diameter. The common-diameter theorem extends
[Schulz's binomial theorem (2016)](https://d-nb.info/1197702695/34)
to unequal probabilities.

For general laws, the argument fixes the means and variances and optimizes
at a selected threshold. A maximizing law can be chosen with at most three
support points. The moment-set reduction and its quadratic-majorant formula
use the methods of [Winkler (1988)](https://doi.org/10.1287/moor.13.4.581)
and [Bertsimas and Popescu (2005)](https://doi.org/10.1137/S1052623401399903).
The remaining finite-dimensional inequalities are proved analytically or
by interval arithmetic.

Extending a two-summand bound to a larger array requires controlling both
the third moments after grouping and the approximation error after
convolution with the remaining sum. The paper gives sufficient conditions
in terms of variance and third moments. Its signed convolution estimate
retains cancellation between positive and negative distribution-function
errors. The [finite supplement](../certificates/finite/README.md) describes
the computational parts of these arguments.

## Small Lyapunov ratio

The [small-Lyapunov paper](../paper/small-lyapunov/small-lyapunov.pdf) proves

```math
L\le\exp(-\exp(97))\quad\Longrightarrow\quad\Delta\le C_{\mathrm E}L
```

for arbitrary independent summands. The proof follows He and Cheng's
combination of Fourier estimates, optimization over summand laws, and local
comparison near two-point distributions. Imported lemmas are cited to
[version 1](https://arxiv.org/abs/2609.06358v1); their
[current manuscript](https://arxiv.org/abs/2609.06358v3) proves eventual
sharpness for i.i.d. samples. For independent arrays, $L$ measures smallness:
adjoining zero summands changes the number of coordinates without changing
the error or the ratio.

A hypothetical small-$L$ violation yields an attained maximizer after a
penalty is subtracted. An explicit remainder estimate of
[Shevtsova (2012)](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf)
controls its third-moment sum. Coordinate variations bound the support
radii. Fourier estimates then select an approximate common spacing, and
Esseen's lattice moment inequality confines the supports to two narrow
intervals, apart from summands with controlled total variance.

Bernoulli labels record which interval contains each summand. Conditional
estimates control the residual displacements. At zero residual variance,
the array has a common diameter and the finite theorem applies. With small
positive residual variance, a one-sided loss dominates the perturbation
errors. With larger residual variance, a conditional Fourier estimate
improves the smoothed bound. Explicit choices in these estimates give the
cutoff displayed above.

## Transport and Bernoulli maximizers

The [transport paper](../paper/transport/transport.pdf) characterizes equality
in the additive bound for the Wasserstein distance between a sum and its
zero-bias transform. Equality holds precisely when all nonconstant summands
are identical centered two-point laws. The proof reduces general laws to
laws on at most three points, proves strictness for genuine three-point
laws, and traces equality through convolution and mixtures.

A lattice-span inequality connects this transport bound to the leading
term in Esseen's asymptotic formula. When a fixed finite block is repeated,
the limiting Berry–Esseen ratio is at most $C_{\mathrm E}$; equality forces
the block's nonconstant coordinates to be identical Esseen Bernoulli laws.
This is an asymptotic statement for a fixed block. It supplies no uniform
estimate for blocks changing with the number of repetitions.

The [Bernoulli paper](../paper/bernoulli-structure/bernoulli-structure.pdf)
classifies every maximizer with a fixed number of coordinates and a common
diameter: its active coordinates are identically distributed. It also
proves that two Bernoulli coordinates of arbitrary diameters have exactly
the one-coordinate optimum, about $0.3703515235$.

For unequal diameters in larger arrays, separating support points can lose
mass at the threshold. A local replacement theorem quantifies that loss
when the total Bernoulli variance is large and the relative diameter
variation is small. Within the optimized binomial family, a phase expansion
explains why zeros sometimes improve a fixed-horizon optimum and proves
that at most four zeros are eventually needed. This classification concerns
binomial candidates; identifying the unrestricted optimum with this family
is a further conjecture. The [Bernoulli supplement](../certificates/structure/README.md)
gives the finite certificates used in the classification arguments.

## Changing arrays and exact structure

The [general structure paper](../paper/small-lyapunov-structure/small-lyapunov-structure.pdf)
assumes $L\to0$ and $\Delta/L\to C_{\mathrm E}$ for changing finite arrays.
It proves that the variance and third-moment sum are concentrated in
coordinates close to one common orientation of Esseen's two-point law.
It also constructs centered coordinatewise laws on a common lattice with
summed squared and cubed Wasserstein errors $o(h^2)$ and $o(h^3)$, where
$h\sim R/(r_{\ast}V)$ and $r_{\ast}=p_{\ast}^2+(1-p_{\ast})^2$ for
$p_{\ast}=(4-\sqrt{10})/2$. These lattice laws may have more than two atoms.
The additive zero-bias transport bound approaches equality.

For an actual global maximizer among all independent laws at its own
finite horizon, the conclusion becomes exact: below an absolute,
ineffective Lyapunov threshold, every original representation consists
of identical nonzero two-point laws and zeros. The finite probability
need not equal $p_{\ast}$. This is a classification of full maximizers;
the separate quantitative small-Lyapunov theorem bounds the value for
every array and gives an explicit cutoff.

The main obstruction to exact structure is phase noise carried by masses
too small for an absolute approximation error to detect. Full marginal
contact equations first confine every support phase. A conditional
first-moment estimate then measures the lost endpoint mass relative to
the actual phase norm. An integer-lattice competitor at the same horizon
has a strictly larger ratio whenever that norm is positive. Contact
curvature reduces the resulting exact lattice laws to adjacent two-point
laws and synchronizes their probabilities.

## Restrictions on hypothetical extremizers

Assume $C_{\mathrm{ind}}\gt C_{\mathrm E}$. The
[extremizer paper](../paper/extremizers/extremizers.pdf) obtains an attained
maximum in a class of countably many independent factors and an infinitely
divisible remainder. The third-moment quantity belongs to the
representation; it is not the third absolute moment of the sum.

Every original marginal of a full finite-horizon maximizer has at most
four atoms. A maximizing representation can be chosen with at most three
atoms per marginal while preserving each marginal's variance and third
absolute moment. The same conclusions hold for the attained relaxed
maximizers considered in the paper. Splitting a nontrivial convolution
factor strictly lowers the original third-moment sum; hence the originals
of an unrestricted relaxed maximizer are convolution indecomposable.
At a fixed horizon, this splitting requires an available zero coordinate.

The variational equations eliminate the Gaussian and
jump parts and force the factor support radii to tend to zero. An infinitely
active sum is nonatomic. The same equations exclude separated finite blocks
and rare Bernoulli tails satisfying the paper's nonoverlap conditions.

For a three-point maximizing representation under the additional bound
$\beta\le18481/25000$ at total variance one, the paper proves
weighted-variance summability. It also excludes specified
one-sided regularity at a maximizing threshold. Dense sequences of rare
Bernoulli gaps force precisely such regularity, yielding further tail
exclusions. These results restrict possible extremizers; they do not prove
that an unrestricted extremizer is a Bernoulli array or contains a large
common-diameter core.
