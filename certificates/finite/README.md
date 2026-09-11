# Certificates for two summands and Bernoulli arrays

The [paper](../../paper/finite-arrays/main.tex) proves a strict bound of
49/120 for two arbitrary independent summands. It proves Esseen's bound
for arrays consisting of at most three Bernoulli coordinates, or entirely
of Bernoulli coordinates with a common support diameter, and a uniform
3/80 deficit for arrays consisting of at most two Bernoulli coordinates. Its
[contact appendix](../../paper/finite-arrays/contact-proofs.tex) gives the
case reductions and inequalities used by the computational proofs.

[current-results.json](current-results.json) identifies the current exact
partitions, arithmetic definitions, receipt hashes, and recorded counts.
The partitions and decision recipes are in
[current-finite.zip](current-finite.zip). The common-diameter,
maximal-atom, and product-odds calculations are recorded in
[results.json](results.json).

| Calculation | Completed arithmetic evidence |
| --- | --- |
| Two-Bernoulli deficit greater than 3/80 | 63,792 leaves covering the three stated roots, evaluated with separately arranged formulas at 320 bits. |
| Three Bernoulli coordinates, arbitrary diameters | Eight recipes, 2,624,314 nodes and 1,312,161 accepted leaves. The bounds were reevaluated at 320 bits, using the producer's mathematical formulas and a separate tree traversal. |
| Triple/pair unique lower event | 15,081 independently evaluated leaves. |
| Triple/pair unique low triangular event | 113,322 independently evaluated leaves. |
| Triple/pair upper-gap collision | 113,154 independently evaluated leaves. |
| Triple/pair lower-gap collision | 154,848 independently evaluated leaves. |
| Two-triple asymmetric collision | 68,963 and 90,233 independently evaluated leaves in its two charts. The positive chart uses the signed divided-difference inequality proved in the appendix. |
| Equality at 49/120 | 427 lower-event leaves and 34,084 upper-gap-collision leaves satisfy strict negative rejection inequalities at 384 bits. |
| Three- and two-representation triangular events | Complete fixed scalar grids, with direct independent evaluations at 320 bits. |
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
recompute numerical inequalities. It also rejects missing, repeated, and
overlapping test partitions. These checks run in `./verify.sh --no-latex`.

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
