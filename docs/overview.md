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
$C_{\mathrm{ind}}<0.44988794$. It combines characteristic-function estimates
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

The finite paper also derives necessary conditions for a maximum over
the full independent class with one three-state coordinate and two
two-state coordinates, when the threshold configurations are exactly
$(0,0,1)$, $(0,1,0)$, and $(1,0,0)$, with only $(0,0,0)$ below.
The two Bernoulli probabilities must agree; the remaining conditions
express the masses as affine functions of total variance and give a
quadratic inequality in that variance.

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

```math
L\le\exp(-\exp(20004))\quad\Longrightarrow\quad
\Delta\le C_{\mathrm E}L.
```

The summands may have different distributions, variances, and supports.
The very small cutoff results from a conservative choice of constants in
the proof.

[He and Cheng (2026)](https://arxiv.org/abs/2609.06358v1) prove the sharp
bound for all sufficiently large i.i.d. samples. The independent-array
argument follows their combination of Fourier estimates, optimization
over summand laws, and comparison near two-point distributions. Here $L$
measures smallness: adjoining zero summands changes the number of
coordinates without changing $L$ or $\Delta$.

A hypothetical violation at small $L$ yields an attained maximizer after
subtracting a penalty that increases with $L$. An explicit remainder
estimate of [Shevtsova (2012)](https://publikacio.uni-eszterhazy.hu/3231/1/AMI_39_from241to307.pdf)
keeps its Lyapunov ratio below twice the original cutoff. Varying one
summand at a time then bounds every support radius by $8L$.

Fourier estimates select an approximate common spacing and give an explicit error
for the sum after adding a uniform variable. The equality case of
Esseen's lattice moment inequality and the variational equations confine
every support to two narrow intervals, apart from small summands whose
total variance is controlled.

The local comparison theorem records which interval contains each
summand with a Bernoulli label. Conditional estimates control the
remaining displacement while permitting unequal probabilities and
unequal spacings. At zero residual variance the array is a common-diameter Bernoulli
sum. For small positive residual variance, a one-sided loss dominates
the perturbation errors. For larger residual variance, a conditional
Fourier estimate bounds interval probabilities on both sides of a
threshold and improves the bound obtained by adding the uniform variable. These comparisons reduce the sharp
inequality to the common-diameter theorem in the finite-array paper.

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
are nonconstant, their sum is nonatomic. The paper also excludes infinitely
many finite blocks whose support diameters are summable and whose distinct
block sums are separated by at least twice the total diameter of all later
blocks. The paper states the precise separation condition and proves these
assertions.

Further restrictions hold when the sum of third absolute moments,
denoted by $\beta$, satisfies $\beta\le18481/25000$ at total variance
one. Write $H_i$ for the diameter of the support of a nonconstant
factor and $v_i$ for its variance. Then

```math
\sum_i\frac{v_i}{H_i^\gamma}<\infty
\qquad(0<\gamma<\sqrt5-1).
```

In particular, $\sum_i\mathbb E|X_i|<\infty$. If all but finitely many
factors have two states, the weighted-variance conclusion holds for
$0<\gamma<\sqrt{21}-3$. The proof combines restrictions on the support
shapes with grouping estimates and a lattice comparison that accounts
for the changes in variance and third moments.

These estimates also restrict the distribution function $F$ at a
threshold $z$ where $F(z)-\Phi(z)=C_{\mathrm{ind}}\beta>0$. Under the
same bound on $\beta$, a maximizing sum cannot satisfy
$F(z+t)=F(z)+\phi(z)t+O(t^2)$ as $t\to0$ through both signs.
For a sum with only finitely many three-state factors, even the
one-sided expansion

```math
F(z)-F(z-h)=\phi(z)h+O(h^{1+\alpha}),\qquad h\downarrow0,
```

is impossible when $\alpha>(\sqrt{21}-3)/6$.
