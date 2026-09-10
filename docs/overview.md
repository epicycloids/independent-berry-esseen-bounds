# Reading the papers

For independent centered real variables with finite third absolute moments,
write

$$
B^2=\sum_i\mathbb E X_i^2>0,\qquad
\ell=\frac{\sum_i\mathbb E|X_i|^3}{B^3},\qquad
\Delta=\sup_x\left|\mathbb P\left\{\frac{\sum_iX_i}{B}\le x\right\}-\Phi(x)\right|.
$$

The conjecture is $\Delta\le C_{\mathrm E}\ell$ for every finite independent
array, where $C_{\mathrm E}=(3+\sqrt{10})/(6\sqrt{2\pi})$. Esseen's sequence
of identically distributed two-point summands attains this constant in the
limit.

For the sharp conjecture, start with the small-Lyapunov paper. Its only
result taken from the finite-results companion is the theorem for Bernoulli
sums with a common support diameter. The companion also treats two summands
with arbitrary diameters and reductions for general arrays. Both papers
can be read independently of the unrestricted upper-bound computation.
The [bounds note](bounds.md) collects the numerical values and restricted
classes.

## The sharp bound for small Lyapunov ratio

The [small-Lyapunov paper](../paper/small-lyapunov/small-lyapunov.pdf) proves
that there is a universal $\ell_0>0$ such that

$$
0<\ell\le\ell_0\quad\Longrightarrow\quad
\Delta\le C_{\mathrm E}\ell.
$$

The summands need not have the same law, variance, or support. The proof is
by compactness and gives no numerical value of $\ell_0$.

[He and Cheng (2026)](https://arxiv.org/abs/2609.06358v1) prove the sharp
bound for i.i.d. samples above a universal sample-size threshold. Their
proof combines Fourier smoothing, optimization over summand laws, and a
comparison near a two-point law. The independent-array proof follows this
strategy, with the Lyapunov ratio replacing sample size as the asymptotic
parameter. Sample size alone cannot serve this purpose for arbitrary
independent arrays: adjoining zero summands changes neither $\ell$ nor
$\Delta$.

The proof assumes that violations exist with $\ell\to0$ and selects
maximizing arrays by adding a smooth penalty in $\ell$. For these arrays,
Fourier estimates and variations of the individual summand laws constrain
all support points. After division by $B\ell$, each support lies either
in two shrinking intervals with a nearly common separation, or in a
shrinking interval about zero. Equality in Esseen's lattice moment
inequality identifies the limiting separation. The positive maximizing
discrepancy fixes the sign of the third moment. These conclusions include
support points whose probabilities tend to zero.

For each summand in the first group, a Bernoulli label records which
interval it occupies. A common translation and scaling writes the sum as
$K+W$, where $K$ is the sum of those labels. The residual $W$ includes the
unequal separations, variation within the intervals, and all the summands
in the second group. Estimates for $W$ conditional on $K$ give an exact
comparison with the common-diameter Bernoulli sum when
$0<\operatorname{Var}(W)\to0$. When $\operatorname{Var}(W)$ is bounded below
by a positive constant, conditional Fourier estimates give positive mass
between the integer levels and a bound strictly below $C_{\mathrm E}$ for
$\Delta/\ell$.
The common-diameter theorem covers $W=0$.

The argument depends on the finite companion's certificate-assisted
common-diameter theorem. It does not require a new numerical calculation.
The general sharp conjecture remains open for Lyapunov ratios bounded away
from zero; the [bounds note](bounds.md) states the remaining range.

## Finite classes and reductions

The [finite-results paper](../paper/finite-arrays/finite-arrays.pdf) proves
$\Delta\le C_{\mathrm E}\ell$ when every summand has at most two support
points and all nondegenerate summands have the same support diameter.
The probabilities may differ, and the constant is optimal over this class.
This extends [Schulz's binomial theorem (2016), Theorem 1](https://d-nb.info/1197702695/34).
For two summands, each supported on at most two points, the paper proves
the strict inequality $\Delta<C_{\mathrm E}\ell$ with arbitrary diameters.
It does not determine the optimal constant for two summands.

For the common-diameter proof, write each summand as $h(B_i-p_i)$, where
$B_i$ is Bernoulli with success probability $p_i$. At a fixed mean and
integer threshold, a negative minimum of the Berry–Esseen deficit can be
reduced to a binomial sum and one Bernoulli variable. A curvature inequality
excludes this minimum when the binomial group has at least two variables;
the two-summand theorem supplies the remaining case. The latter theorem
checks the two interior atoms directly and uses a separate estimate at the
support endpoints.

For independent Bernoulli variables $B_i$ and positive weights $h_i$,
[Yehuda and Yehudayoff's product-measure antichain bound, Theorems 1 and 2](https://arxiv.org/html/2509.01160v1)
implies that every atom of $\sum_i h_iB_i$ has probability at most
$\max_k\mathbb P(N=k)$, where $N=\sum_iB_i$. The proof uses an increasing
coupling of Bernoulli vectors conditioned on $N=k$, constructed by
[Broman, van de Brug, Kager, and Meester (2012), Lemma 2.2](https://alea.impa.br/articles/v9/09-17.pdf).
The finite paper adds an equality condition: if $N$ has a unique mode
strictly between its support endpoints, equality forces all weights to
agree. It also extends the comparison to short intervals.

Every strict counterexample to the general conjecture can be replaced by
one with the same number of summands and at most three support points per
summand. The proof fixes each summand's mean and variance and maximizes its
contribution at a selected threshold. A cubic penalty makes the maximum
attainable; an extreme maximizing law has at most three support points.
The paper proves the moment-set argument and an attained quadratic-majorant
formula for this objective. General results on these methods are given by
[Winkler (1988)](https://doi.org/10.1287/moor.13.4.581) and
[Bertsimas and Popescu (2005)](https://doi.org/10.1137/S1052623401399903).
Two finite problems remain: the sharp bound for two-point arrays with
arbitrary diameters, and a three-point mixture inequality. Proving both
would establish the conjecture for general independent arrays.
The [computational supplement](../certificates/finite/README.md) describes
the certificate calculations used in the finite-class proofs.

## The unrestricted upper bound

The [upper-bound paper](../paper/independent-berry-esseen.pdf) proves
$C_{\mathrm{ind}}<0.454$. It combines moment bounds for characteristic
functions, a zero-bias comparison following
[Tyurin (2009)](https://arxiv.org/abs/0912.0726v1), and Prawitz smoothing.
It also compares the sum with a sum of symmetric two-point variables having
the same variances, an approximation studied by
[Mattner and Shevtsova (2019)](https://alea.impa.br/articles/v16/16-19.pdf).
Because the reference characteristic function is real, its difference from
the Gaussian characteristic function can enter the smoothing integral with
its sign. Analytic estimates cover small and large third-moment ratios;
an interval calculation covers the intermediate range. This paper writes $L$ for
the quantity denoted by $\ell$ above.

The [verification guide](../verify/README.md) describes the partition,
arithmetic assumptions, and commands for checking the calculation.
