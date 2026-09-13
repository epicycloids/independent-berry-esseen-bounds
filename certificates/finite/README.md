# Certificates for two summands and Bernoulli arrays

The [paper](../../paper/finite-arrays/main.tex) proves a strict bound of
49/120 for two arbitrary independent summands. It proves Esseen's bound
for arrays consisting of at most three Bernoulli coordinates, or entirely
of Bernoulli coordinates with a common support diameter, and $C_{\mathrm E}L-\Delta>3/80$ for arrays of at most two Bernoulli
coordinates with positive total variance. Here $L$ is the Lyapunov ratio
and $\Delta$ is the normal-approximation error. Its
[contact appendix](../../paper/finite-arrays/contact-proofs.tex) gives the
case reductions and inequalities used by the computational proofs.
A triple means a three-point summand and a pair a two-point summand.
The appendix defines the support-index event labels used below.

[current-results.json](current-results.json) identifies the current exact
partitions, arithmetic definitions, receipt hashes, and recorded counts.
The partitions and decision recipes are in
[current-finite.zip](current-finite.zip). The common-diameter,
maximal-atom, and product-odds calculations are recorded in
[results.json](results.json).

| Calculation | Completed arithmetic evidence |
| --- | --- |
| Two Bernoulli summands: $C_{\mathrm E}L-\Delta>3/80$ | 63,792 leaves covering the two unit parameter cubes and the support-endpoint unit square, evaluated with separately arranged formulas at 320 bits. |
| Three Bernoulli coordinates, arbitrary diameters | Eight recipes, 2,624,314 nodes and 1,312,161 accepted leaves. The bounds were reevaluated at 320 bits, using the producer's mathematical formulas and a separate tree traversal. |
| Triple/pair `110`, unique representation (`lower110`) | 15,081 independently evaluated leaves. |
| Triple/pair `210`, unique low/high representation (`low210`) | 113,322 independently evaluated leaves. |
| Triple/pair `221`, upper-gap collision (`collision221`) | 113,154 independently evaluated leaves. |
| Triple/pair `210`, lower-gap collision (`collision210`) | 154,848 independently evaluated leaves. |
| Two triples `221` or `320`, two representations (`tt221_negative`, `tt221_positive`) | 68,963 and 90,233 independently evaluated leaves in its two charts. The positive chart uses the signed divided-difference inequality proved in the appendix. |
| Equality at 49/120 | 427 lower-event leaves and 34,084 upper-gap-collision leaves satisfy strict negative rejection inequalities at 384 bits. |
| Two triples `321`, three or two representations (`three_owner321`, `adjacent_owner321`, `endpoint_owner321`) | Complete fixed scalar grids, with direct independent evaluations at 320 bits. |
| Common-diameter curvature | 224,419 accepted boxes at 192 bits; the primary program was also rerun. |
| Scalar maximal-atom inequality | 129,024 slabs; a separate coefficient-sum implementation checked them at 256 and 384 bits. |
| Product-odds inequality | Exact Bernstein coefficients in the quadratic field Q(sqrt(5)), with a separate algebraic calculation. |

“Independently evaluated” in the contact rows means a different formulation
of the necessary equations, using direct moments of the marginal laws
where indicated in the paper. These implementations use the same Arb
arithmetic library. A repeated evaluation using the same formulas is
identified as such in the table.

From the public repository root, run:

```sh
uv run --frozen python verify/finite/check.py
```

The default checks hashes, exact tree geometry, leaf boxes, all eight
Bernoulli recipes, and the complete scalar-grid index sets. It does not
recompute numerical inequalities. It rejects missing leaves, repeated partition paths, and overlapping
interiors; closed boxes may share boundary faces. These checks run in `./verify.sh --no-latex`.

To recompute the frozen current finite inequalities, without creating or
subdividing cells:

```sh
uv run --frozen python verify/finite/current/check.py --arithmetic
```

A single recorded family can be selected with `--case`, for example
`--case margin`, `--case tt221_positive`, or `--case scalars`. The full
arithmetic run is substantial. The verifier accepts a leaf only if the
outward interval calculation proves the required strict inequality or a
necessary contact condition is impossible. Enclosures printed in receipts
are not trusted as arithmetic inputs.

The broader command below also reconstructs the common-diameter curvature,
scalar atom, and exact odds calculations. The curvature calculation uses
its supplied subdivision algorithm; its complete run succeeds only when
all cells are certified. Exhausting a limit is a failure.

```sh
uv run --frozen python verify/finite/check.py --full --independent
```

The locked environment uses Python 3.12, python-flint 0.8.0, and SymPy
1.14.0. Arb and SymPy are described by
[Johansson (2017)](https://doi.org/10.1109/TC.2017.2690633) and
[Meurer et al. (2017)](https://doi.org/10.7717/peerj-cs.103).
