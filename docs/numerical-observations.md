# Numerical observations and structural conjectures

6 October 2026

The finite-horizon independent Berry–Esseen problem asks which centered
summand laws maximize the normal-approximation error divided by their
normalized third-moment sum. The available computations give a detailed
candidate: identical two-point laws, with some coordinates set to zero.
Their optimal active count changes irregularly with the horizon. Phase
analysis explains that pattern within the candidate family; the remaining
structural question is whether arbitrary independent laws can do better.

For independent centered real laws, write

```math
V=\sum_i\mathbb E X_i^2>0,\qquad
R=\sum_i\mathbb E|X_i|^3,\qquad
L=R/V^{3/2},\qquad
C=\frac{d_K(\sum_iX_i/\sqrt V,Z)}{L}.
```

Here $Z$ is standard normal. The price $R$ belongs to the separate
summands: replacing them by their grouped sum can change it. Let $C_n$ be
the supremum over $n$ summands, allowing deterministic zeros. The
[bounds note](bounds.md) states the known bounds on
$C_{\mathrm{ind}}=\sup_n C_n$.

## Binary cores and zero padding

Take $m$ independent Bernoulli variables $B_i$ with upper probability $p$,
put $q=1-p$, and use

```math
X_i=\frac{B_i-p}{\sqrt{mpq}},\quad i\le m,
\qquad X_i=0,\quad m<i\le n.
```

The total variance is one. Each active support has gap
$h=1/\sqrt{mpq}$, and $L=(p^2+q^2)h$. For
$K\sim\operatorname{Bin}(m,p)$ and $z_k=(k-mp)/\sqrt{mpq}$, the
candidate ratio is

```math
R_m(p)=\frac{\sqrt{mpq}}{p^2+q^2}
\max_{0\le k\le m}
\{\mathbb P(K\le k)-\Phi(z_k),
\Phi(z_k)-\mathbb P(K<k)\}.
```

The two terms check the closed and left sides of every atom. Define
$B_m=\max_pR_m(p)$ and $P_n=\max_{m\le n}B_m$. Padding adds coordinates
without changing the ratio, so a record core remains available at later
horizons. The computations estimate these binomial quantities; they do
not search all independent marginal laws.

For every $1\le m\le100$, a grid of 4096 probabilities in $(0,1/2]$
uses reflection to cover the two orientations. Every atom and both sides
are evaluated. Visible local peaks are refined within their neighboring
grid intervals, and the selected stationary roots are then refined at
60 decimal places. The largest stored stationarity residual among these
100 refinements is about $2.61\times10^{-59}$. Small residuals support
the reported root values; the grid search does not certify the number of
roots in every interval.

The scan continues through $m=3000$ by checking seven central
cells for a closed CDF threshold per core. A separate control at
$101\le m\le450$ checks all central threshold cells, both CDF sides, cell endpoints and
stationary-root brackets on a mesh of 17 points per cell. The control
agrees with the selected cells at the first appearances of three and
four zeros. [Schulz's binomial bounds](https://d-nb.info/1197702695/34)
exclude the noncentral competitors for these returned ratios, which
exceed .4045; the finite mesh still supplies no interval root-count
certificate. Thus the numerical coverage above 100 remains a
central-cell search.

The following rows show the best observed padded cores. Displayed
probabilities and ratios are rounded from the 60-decimal refinements;
the support gaps are formula conversions at the displayed core parameters.

| Horizon $n$ | Active core $m$ | Zeros | Upper probability $p$ | Gap $h$ | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 | 1 | 1 | .390958631281333 | 2.049327 | .370351523474379 |
| 3 | 3 | 0 | .431607387750433 | 1.165657 | .398156326353304 |
| 10 | 10 | 0 | .429865759681911 | .638771 | .406507736905768 |
| 50 | 49 | 1 | .414895409236786 | .289945 | .409145180658100 |
| 100 | 99 | 1 | .417421147624576 | .203806 | .409451697372009 |
| 412 | 409 | 3 | .418879394858355 | .100221 | .409665276377194 |
| 413 | 409 | 4 | .418879394858355 | .100221 | .409665276377194 |

The first observed three-zero and four-zero horizons are 412 and 413.
Both retain core 409. Its refined ratio is
.40966527637719430986; the refined ratio of the 412-coordinate core is
.40966527592709729990. The difference is small enough that refinement
matters when deciding whether 412 sets a new record.

The [data file](../data/empirical-observations.json) contains the scan rows,
the 100 refined roots, the event refinements at 407–413, selected padded
horizons, and the comparisons discussed below. Its column descriptions
distinguish a core's own ratio from the active count selected by the
padded record envelope.

## The phase explanation

The optimizing probability approaches Esseen's
$p_{\ast}=(4-\sqrt{10})/2$. The selected binomial threshold has
$mp-k\to c_{\ast}=(20+\sqrt{10})/72$. A central threshold therefore cannot
remain at its preferred phase for every integer $m$ while keeping $p$
at its preferred value.

For globally optimized binomial cores, the uniform phase expansion in the
[Bernoulli proof, Section “Binomial phase and the optimal number of zero coordinates”](../paper/bernoulli-structure/bernoulli-structure.pdf) gives

```math
B_m=C_{\mathrm E}-\frac{Q_{\ast}}{m}
-\frac{H_{\ast}+\alpha\operatorname{dist}(mp_{\ast}-\delta_0,\mathbb Z)^2}{m^2}
+O(m^{-3}),
```

where $H_{\ast}$ is independent of $m$ and

```math
C_{\mathrm E}=\frac{3+\sqrt{10}}{6\sqrt{2\pi}},\quad
Q_{\ast}=\frac{\phi(0)(867+287\sqrt{10})}{25920},\quad
\alpha=\frac{\phi(0)(60+19\sqrt{10})}{30},\quad
\delta_0=\frac{2590+577\sqrt{10}}{12960}.
```

Schulz's dissertation supplies the binomial bound, central-threshold
reduction and cubic phase polynomial. The uniform expansion permits comparison
of cores a bounded number of horizons apart. Put
$\eta=Q_{\ast}/\alpha=(2510-747\sqrt{10})/8640$ and
$\theta_n=\operatorname{round}(np_{\ast}-\delta_0)-(np_{\ast}-\delta_0)$.
To leading order, padding by $j$ minimizes

```math
\eta j+\operatorname{dist}(\theta_n+jp_{\ast},\mathbb Z)^2.
```

The increasing term prices a smaller core, while the second term
measures its phase loss. The lower envelope has only five active choices,
$j=0,1,2,3,4$. Consequently the optimized IID family eventually uses at
most four zeros. Their respective frequencies are
$p_{\ast}-\omega,p_{\ast}-\omega,\sqrt{10}-3,\omega,\omega$, with
$\omega=(7773-2447\sqrt{10})/8640\simeq.00404$.
Three and four zeros each occur only about four times per thousand
horizons asymptotically.

Finite-horizon phase boundaries move at order $1/n$. The next expansion
and a quadratic-field separation argument give a unique winning IID
core at every sufficiently large horizon and exact eventual record gaps
$\lbrace2,3,5\rbrace$. These are qualitative mathematical results within the
binomial family, with no effective first horizon. The limiting unshifted
policy alone does not decide every finite horizon. At 147 and 148, for
example, it predicts three and four zeros; the drift correction and the
stored high-precision central candidates give one and two. At 409–414,
all three comparisons give the run $0,1,2,3,4,0$.

## Local stability and objective comparisons

Equal gaps create collisions: many Bernoulli outcomes contribute to one
sum atom. A small gap perturbation separates those outcomes, and a
threshold retaining the whole atom must follow an extreme displaced
outcome. This produces a cusp in the ratio as the gaps vary. Probability
variations have a smooth second-order Hessian instead. Both effects
help explain why a search can settle on a stable binary core.

Unequal gaps can also produce a stable collision. The explicit example
uses four Bernoulli variables with gaps $(1,1,1,2)$, probabilities
$(p,p,p,r)$, and closed raw threshold 2, where
$p=.57275160060602450695\ldots$ and
$r=.33747296608998501783\ldots$. The outcomes at the attaining threshold
are the pair outcomes of the three gap-1 variables and the gap-2 singleton.
The [fixed-box verifier](../verify/observations/local_trap.py) rechecks the
[published interval certificate](../certificates/observations/local-trap.json)
at 320 bits in a probability box of radius $10^{-28}$ about the stored
centers. It proves that the box contains one exact stationary root,
that the selected CDF threshold is unique, that the complete
four-probability Hessian is negative definite, and that both
one-coordinate cubic majorants hold on the whole real line. The ratio
is enclosed by
$0.28248894967964153979405802\pm5.54\times10^{-27}$.

The local gap argument is geometric. The four colliding outcome vectors
are $(1,1,0,0)$, $(1,0,1,0)$, $(0,1,1,0)$, and $(0,0,0,1)$.
They form a tetrahedron in the hyperplane perpendicular to
$h=(1,1,1,2)$. The certificate gives positive weights on all four
vertices, approximately $(.31153,.31153,.31153,.06541)$, whose weighted
average is the first-order gap coefficient $\beta$. For a gap direction
$d$, the directional ratio derivative is a positive factor times
$\beta\cdot d-\max_v v\cdot d$. This is strictly negative unless
$d$ is proportional to $h$: an interior weighted average is below the
largest vertex value whenever those values differ. Positive outcome
masses and the strict CDF margins exclude nearby thresholds retaining
only part of the separated collision. On a chart transverse to scale,
the gap loss is therefore at least $c\lVert\delta h\rVert$, uniformly for nearby
probabilities; at the original gap pattern, stationarity and the negative
Hessian give a probability loss of at least $c'\lVert\delta p\rVert^2$.
Together these estimates prove a strict local maximum among four
independent two-point laws, modulo common scale.

Two symmetric Bernoulli variables plus two zeros have ratio
$\sqrt2/4\gt .35355$, and the same verifier certifies the strict comparison.
Local necessary systems can therefore retain an unequal-gap array that
loses an elementary whole-array comparison. The
[Bernoulli paper](../paper/bernoulli-structure/bernoulli-structure.pdf)
uses such objective comparisons in its classification theorems.

A second control tests a finite binary core with a signed
compound-Poisson block. For independent Poisson counts of intensities
$\ell_+,\ell_-$, its raw sum, variance and limiting original price are

```math
S_{\rm raw}=\operatorname{Bin}(M,p)+g(N_+-N_-),\qquad
V=Mpq+g^2(\ell_++\ell_-),\qquad
R=Mpq(p^2+q^2)+g^3(\ell_++\ell_-).
```

Center $S_{\rm raw}$ by its mean and use the ratio $V^{3/2}d_K/R$.
The Poisson term is charged by the cubic intensity price inherited from
many rare original summands. The pilot examines
$M\in\lbrace1,2,3,5,10,20\rbrace$ and $g\in\lbrace1/2,1,2\rbrace$, with positive,
negative and both-sign intensity stages, with $p\in[.05,.5]$ and positive
total intensity in $[10^{-6},12]$. All 54 stochastic stages reached their
evaluation caps of 80 or 160 evaluations. Six finalists receive 60-decimal rechecks with analytic
Poisson-tail bounds at their reported parameters.

Every reported mixed best loses to its zero-intensity core at the same
reported $p$. The highest mixed value, at $M=20,g=1$, is about
.4083120738147323, versus .4083121893578447 for that same-probability
core. Apparent gains over coarse probability controls came from refining
$p$. The comparison isolates the effect of adding intensity; comparison
only with a coarse grid does not. The saved runs give candidate values
within these profiles; reaching an evaluation cap does not establish a
profile's maximum.

The [transport paper](../paper/transport/transport.pdf) supplies another
structural diagnostic. For the sum $S$ and its zero-bias law $S^*$,

```math
B_{\rm add}=\frac{2VW_1(S,S^*)}{R}\le1,
```

with equality precisely for padded IID binary originals. A finite
inequality $C\le C_{\mathrm E}B_{\rm add}$ would connect the diagnostic
directly to the objective, but such a universal bridge is false.

An exact one-law control retains a standard Gaussian outside $[-a,a]$
and replaces its interior by three atoms at $-a,0,a$, preserving mean
and variance. Write $P=\mathbb P(|G|\le a)$ and
$r=\mathbb E[G^2;|G|\le a]/(a^2P)$; the three masses are
$Pr/2,P(1-r),Pr/2$. Its discrepancy and transport satisfy

```math
D=\frac{P(1-r)}2,\qquad
W_1(X_a,X_a^*)=\frac{aP}{2}\{(1-r)^2+r^2\},\qquad
\frac{C}{B_{\rm add}}
=\frac{1-r}{2a\{(1-r)^2+r^2\}}\sim\frac3{5a}.
```

As $a\to0$, $C\sim a/6$ while $B_{\rm add}\sim5a^2/18$.
Although $C$ tends to zero, $C/B_{\rm add}$ tends to infinity, ruling
out a bound $C\le K B_{\rm add}$ with an absolute finite $K$.
Exact transport saturation identifies a representation; deriving it
from Berry–Esseen global optimality requires an additional argument.

## The remaining conjectures

The numerical record envelope and the proved reductions suggest two
finite-horizon statements of different strength:

1. **Value:** $C_n=P_n=\max_{m\le n}B_m$ for every horizon with zeros
   allowed.
2. **Representation:** every attaining independent array is padded IID
   binary, up to scale, reflection and permutation.

The computations examine the right side of the value conjecture.
Equality of values would still leave the possibility of other attaining
representations. Requiring every individual variance to be positive
preserves the padded supremum by small-noise approximation, but can
remove its attainment.

The common-gap classification proves the representation statement in
that class. A further theorem shows that if each row of a changing
sequence is an actual full fixed-horizon global maximizer and $L\to0$,
then every original representation is eventually exactly padded IID
binary. It imposes no initial binary-support or common-gap assumption.
What remains is to connect an arbitrary full global maximizer to those
structural regimes, retaining whole-law competitors and the original
additive price. In particular, increasing the horizon by padding does
not itself make $L$ small.

The value conjecture would imply the sharp independent constant
conjecture $C_{\mathrm{ind}}=C_{\mathrm E}$, because every finite
binomial optimum is strictly below $C_{\mathrm E}$ and the binomial
optima approach it.
