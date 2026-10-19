# Verification Guide

Navigation: [Repository Root](../README.md) · [Overview & Reading Guide](../docs/overview.md) · [Summary of Bounds](../docs/bounds.md) · [Scalar Certificates](../certificates/scalar/README.md) · [Finite-Array Certificates](../certificates/finite/README.md) · [Bernoulli Structure Certificates](../certificates/structure/README.md) · [Interval Kernel](kernel/README.md)

---

## 1. Running the Default Verification Suite

From the repository root, run:

```sh
./verify.sh --no-latex
```

This command runs in the pinned Python 3.12 environment managed by `uv` and performs the following checks:

1. **Manifest integrity (`SHA256SUMS`):** Verifies the SHA-256 digest of every file in the repository.
2. **Unit and regression tests (`verify/test_*.py`):** Verifies the partition validator, evaluator parameter bindings, nested moment-proof consumers, and failure-rejection controls.
3. **Moment-domain cover (`verify/check_saved.py` and `verify/check_replays.py`):** Reconstructs all 1,670 closed intervals of the Lyapunov ratio $`L`$ and their binary partition trees in [`certificates/cover.json`](../certificates/cover.json) and `certificates/cover.zip`, verifies every infeasibility exclusion and nested moment proof, checks the Gaussian defect interpolation tables in `certificates/gaussian-*.zip`, and confirms that every accepted region satisfies

   $`\displaystyle \frac{\Delta}{L} \le U_* = \frac{1013039690742047}{2251799813685248} < 0.44988,`$

   together with the small-$`L`$ and large-$`L`$ analytic complements (`verify/scalars.py`) and manuscript constants (`verify/paper_values.py`).
4. **Scalar majorant certificates (`verify/check_scalar_certificates.py`):** Verifies the polynomial, logarithmic, directional, disk, and cosine-comparison scalar certificates used in [Manuscript 1](../paper/independent-berry-esseen.pdf); see [Scalar Certificates](../certificates/scalar/README.md).
5. **Two-summand and three-Bernoulli certificates (`verify/finite/check.py`):** Verifies the partition trees, recipes, and exact Bernstein odds certificates used in [Manuscript 2](../paper/finite-arrays/finite-arrays.pdf); see [Finite-Array Certificates](../certificates/finite/README.md).
6. **Bernoulli structure certificates (`verify/structure/check.py`):** Verifies the common-diameter and two-Bernoulli partitions and rational benchmarks used in [Manuscript 5](../paper/bernoulli-structure/bernoulli-structure.pdf); see [Bernoulli Structure Certificates](../certificates/structure/README.md).
7. **Gaussian window certificate (`paper/extremizers/gaussian_window_certificate.py`):** Executes the exact rational arithmetic verification for Lemma 5.3 of [Manuscript 7](../paper/extremizers/extremizers.pdf).
8. **Four-summand local maximizer certificate (`verify/observations/local_trap.py`):** Re-evaluates the 320-bit Arb interval certificate [`certificates/observations/local-trap.json`](../certificates/observations/local-trap.json) described in [Numerical Studies & Conjectures](../docs/numerical-observations.md).

Running `./verify.sh` without `--no-latex` additionally compiles all eight LaTeX manuscripts into `dist/` via `./build.sh`.

---

## 2. Structure of the Moment-Domain Certificate (Manuscript 1)

After normalizing total variance to $`V = 1`$, the upper-bound calculation in [Manuscript 1](../paper/independent-berry-esseen.pdf) parameterizes each independent array by four moment quantities:

- $`L = \sum_i \mathbb{E}|X_i|^3`$, the Lyapunov ratio;
- $`d = \max_i \mathbb{E} X_i^2`$, the largest individual summand variance;
- $`b = \mathbb{E}|X_1|^3`$, the third absolute moment of a summand achieving variance $`d`$; and
- $`\tau = \sum_i (\mathbb{E} X_i^2)^{3/2}`$.

The index [`certificates/cover.json`](../certificates/cover.json) divides $`[L_-, L_+]`$ into 1,670 closed intervals of $`L`$ (called *moment bands*), whose partition trees are stored in `certificates/cover.zip`. Summary totals are recorded in [`result.json`](../result.json):

- **1,670** closed intervals of $`L`$ covering $`[L_-, L_+]`$ without gaps;
- **472,854** total tree nodes across the outer moment partitions;
- **235,659** accepted leaf boxes (comprising **143,618** individually evaluated leaf bounds and **92,041** leaves bounded by a uniform enclosure valid on the entire $`L`$ band); and
- **1,468** leaves excluded by necessary moment inequalities ($`b \ge d^{3/2}`$, $`d^{3/2} \le \tau \le \min\{L - b + d^{3/2}, \sqrt{d}\}`$, and $`\tau \le L`$).

Within **10,238** of the accepted leaf boxes (counted under `complete_owner_proofs` in `result.json`), a nested partition tree subdivides the feasible interval of $`b`$ or $`d`$ and records fixed evaluator inputs—including signed threshold partitions on $`[-2, 2]`$ or $`[-69/32, 69/32]`$, two-disk support interpolations, or split Cauchy–Schwarz forcing bounds—that certify $`\Delta / L \le U_*`$ across the leaf.

---

## 3. Recomputing Numerical Enclosures

The default `./verify.sh --no-latex` run checks the completeness of every partition tree and verifies the stored upper bounds against $`U_*`$, while reading the precomputed leaf integral enclosures from the certificate archives. To recompute the leaf integrals directly from the interval and ball-arithmetic evaluators in [`verify/kernel/`](kernel/README.md), use `verify/replay.py`:

```sh
# Recompute all leaves (including nested moment proofs) in band 0
uv run --frozen python verify/replay.py --band-index 0

# Recompute a slice of leaves in band 0
uv run --frozen python verify/replay.py --band-index 0 --first-leaf 0 --max-leaves 10
```

For each reevaluated leaf, `verify/replay.py` verifies that individually stored leaf bounds are reproduced identically and that containing-band bounds are not exceeded.

Commands for recomputing the scalar, finite-array, and Bernoulli-structure interval enclosures are documented in their respective guides:

- [Scalar Certificates (`verify/check_scalar_certificates.py --replay ...`)](../certificates/scalar/README.md)
- [Finite-Array Certificates (`verify/finite/current/check.py --arithmetic`)](../certificates/finite/README.md)
- [Bernoulli Structure Certificates (`verify/structure/check.py --arithmetic`)](../certificates/structure/README.md)
- [Four-Summand Local Maximizer (`verify/observations/local_trap.py --json`)](../docs/numerical-observations.md)
