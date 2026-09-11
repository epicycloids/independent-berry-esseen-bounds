# Reading the papers

The problem is to determine the least constant $C_{\mathrm{ind}}$ in the
Berry–Esseen inequality for independent real summands. The
[bounds note](bounds.md) defines the normalized third-moment sum $L$ and
normal-approximation error $\Delta$, and collects the numerical results.
The conjectured constant is Esseen's
$C_{\mathrm E}=(3+\sqrt{10})/(6\sqrt{2\pi})$.

The four papers address the numerical upper bound, finite arrays, small
Lyapunov ratio, and extremal distributions. For the sharp conjecture,
start with the finite-array paper and then the small-Lyapunov paper. The
upper-bound computation can be read independently of both.

## The unrestricted upper bound

The [upper-bound paper](../paper/independent-berry-esseen.pdf) proves
$C_{\mathrm{ind}}<0.44995$. It combines characteristic-function estimates
and Prawitz smoothing with the zero-bias comparison developed by
[Tyurin (2009)](https://arxiv.org/abs/0912.0726v1). The zero-bias
transformation was introduced by
[Goldstein and Reinert (1997)](https://doi.org/10.1214/aoap/1043862419);
the sharp mean-distance estimate was obtained independently by
[Goldstein (2010)](https://doi.org/10.1214/10-AOP527) and Tyurin.
The proof also compares the
sum with a sum of symmetric two-point variables having the same variances,
as in [Mattner and Shevtsova (2019)](https://alea.impa.br/articles/v16/16-19.pdf).
The reference characteristic function is real, so the smoothing integral
can retain the sign of its difference from the Gaussian characteristic
function.

Analytic estimates cover small and large $L$. In the intermediate range,
an interval computation bounds the smoothing integral over a finite
partition of the moment parameters. The
[verification guide](../verify/README.md) describes the certificate and
distinguishes checks of stored records from numerical reevaluation.

## Finite arrays and extension to larger sums

The [finite-array paper](../paper/finite-arrays/finite-arrays.pdf) proves
$C_2<49/120$ for two arbitrary summands. It proves the conjectured bound
for three Bernoulli summands of arbitrary diameters and for any number
of Bernoulli summands with a common diameter. The latter result extends
[Schulz's binomial theorem (2016)](https://d-nb.info/1197702695/34) to
unequal probabilities.

The general-law reduction fixes each summand's mean and variance and
optimizes at a selected threshold. An extreme maximizing law has at most
three support points. The proof develops the moment-set argument and its
quadratic-majorant formula, using the methods of
[Winkler (1988)](https://doi.org/10.1287/moor.13.4.581) and
[Bertsimas and Popescu (2005)](https://doi.org/10.1137/S1052623401399903).
The remaining finite-dimensional inequalities are treated analytically
and by interval arithmetic.

Applying the two-summand theorem to a larger array requires bounds for
both the change in third absolute moments under grouping and the error
after convolution with the remaining sum. The finite-array paper derives
sufficient conditions from these bounds, including the variance thresholds
in the [bounds note](bounds.md). Its signed convolution estimate retains
cancellation between positive and negative distribution-function errors.

The [finite-results supplement](../certificates/finite/README.md) gives
the computational inputs and verification commands.

## The sharp bound for small Lyapunov ratio

The [small-Lyapunov paper](../paper/small-lyapunov/small-lyapunov.pdf) proves
that there is a universal $L_0>0$ for which $\Delta\le C_{\mathrm E}L$
whenever $0<L\le L_0$. The summands need not share a distribution, variance,
or support. The proof is by compactness and supplies no numerical value
of $L_0$.

[He and Cheng (2026)](https://arxiv.org/abs/2609.06358v1) prove the sharp
bound for all sufficiently large i.i.d. samples. The independent-array
argument follows their combination of Fourier estimates, optimization
over summand laws, and comparison near two-point distributions. Here $L$
takes the place of sample size: adjoining zero summands changes the
number of coordinates without changing $L$ or $\Delta$.

Assuming violations with $L\to0$, the proof selects maximizing arrays
and constrains their summand supports to small neighborhoods of a
common-diameter Bernoulli configuration. A Bernoulli label records which
neighborhood contains each summand. Conditional estimates for the
remaining displacement reduce the inequality to the common-diameter
theorem from the finite-array paper. No new numerical calculation is
needed for this reduction.

## Restrictions on extremizers

Assume that $C_{\mathrm{ind}}>C_{\mathrm E}$. The
[extremizer paper](../paper/extremizers/extremizers.pdf) studies
representations that attain the ratio $C_{\mathrm{ind}}$ in its enlarged
class of independent convolutions. A representation consists of countably
many independent factors and an infinitely divisible remainder, with
separate Gaussian and jump terms. The additive third-moment quantity is
part of the representation; it is not the third absolute moment of the sum.

A maximizing representation can be chosen with at most three support
points in each independent factor. For representations with this property,
the variational equations force the Gaussian and jump terms to vanish
and the factor support radii to tend to zero. If infinitely many factors
are nonconstant, their sum is nonatomic. A condition on the finite support
configurations that can complete to a maximizing threshold also excludes
infinitely many exactly separated finite blocks. The paper defines the
enlarged class and the separation condition and proves these assertions.
