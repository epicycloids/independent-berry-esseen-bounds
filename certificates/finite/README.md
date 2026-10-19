# Finite-Array Certificates (Manuscript 2)

Navigation: [Repository Root](../../README.md) · [Verification Guide](../../verify/README.md) · [Bernoulli Structure Certificates](../structure/README.md)

---

This directory contains the partition certificates, decision recipes, and exact polynomial records supporting [Manuscript 2: *Berry–Esseen bounds for two summands and Bernoulli arrays*](../../paper/finite-arrays/finite-arrays.pdf) ([`paper/finite-arrays/main.tex`](../../paper/finite-arrays/main.tex) and [`paper/finite-arrays/contact-proofs.tex`](../../paper/finite-arrays/contact-proofs.tex)).

## Summary of Certified Calculations

| Theorem / Proposition in Manuscript 2 | Mathematical domain and inequality | Certificate summary |
| --- | --- | --- |
| **Theorem 1.3 & Lemma 2.3** (Two-Bernoulli margin) | $`C_{\mathrm{E}} L - \Delta > 3/80`$ on the two interior-atom parameter cubes and the support-endpoint square $`[0, 99/100]^2`$. | 63,792 accepted leaves evaluated at 320-bit Arb precision (`margin`). |
| **Theorem 1.2(i)** (Three Bernoulli summands, arbitrary diameters) | $`C_{\mathrm{E}} L - \Delta \ge 0`$ across all eight ordered-diameter and interior-atom configurations. | 8 subdivision/projection recipes with 2,624,314 total nodes and 1,312,161 accepted leaves, reevaluated at 320-bit Arb precision (`bernoulli_three`). |
| **Proposition 6.4** (Triple/pair contact case `110`, unique representation) | Exclusion of necessary contact system in Section B.11 (`lower110`). | 15,081 accepted leaves (plus 427 leaves for strict exclusion at $`C = 49/120`$ at 384 bits). |
| **Proposition 6.4** (Triple/pair contact case `210`, low/high representation) | Exclusion of necessary contact system in Section B.10 (`low210`). | 113,322 accepted leaves at 320 bits. |
| **Proposition 6.4** (Triple/pair contact case `221`, upper-gap collision) | Exclusion of necessary contact system in Section B.13 (`collision221`). | 113,154 accepted leaves at 320 bits (plus 34,084 leaves for strict exclusion at $`C = 49/120`$ at 384 bits). |
| **Proposition 6.4** (Triple/pair contact case `210`, lower-gap collision) | Exclusion of necessary contact system in Sections B.14–B.16 (`collision210`). | 154,848 accepted leaves at 320 bits. |
| **Proposition 6.4** (Two-triple contact cases `221` and `320`, two representations) | Exclusion of necessary contact systems in Sections B.27–B.33 (`tt221_negative`, `tt221_positive`). | 68,963 leaves (negative middle-contact chart) and 90,233 leaves (positive middle-contact chart) at 320 bits. |
| **Proposition 6.4** (Two-triple contact case `321`, two or three representations) | Scalar exclusions in Sections B.20–B.26 (`three_owner321`, `adjacent_owner321`, `endpoint_owner321`). | Complete Cartesian grids (4,096 + 4,124 + 3,840 cells) evaluated with 320-bit Arb Taylor series. |
| **Lemma 4.3** (Common-diameter one-exception curvature) | Positivity of the curvature lower bound on $`V \ge 1/5`$, $`pq \le W \le 3`$ at coefficient $`c_* = 199/500 < C_{\mathrm{E}}`$. | 275,283 accepted boxes from 550,563 processed boxes at 192 bits (shared with [Bernoulli Structure Certificates](../structure/README.md)). |
| **Proposition A.2** (Bernoulli support-endpoint odds inequality) | Positivity of tensor Bernstein coefficients for the pair-merge and pair-chord quotients over $`\mathbb{Q}(\sqrt{5})`$. | Exact symbolic Bernstein tables in `results.json`, verified by `verify/finite/support_odds/`. |

---

## Running the Verifier

1. **Fast structural and exact-algebra check (default in `./verify.sh --no-latex`):**
   ```sh
   uv run --frozen python verify/finite/check.py
   ```
   Verifies file hashes, reconstructs every binary partition tree from the root domain (confirming that all leaf boxes partition the domain without gaps or interior overlaps), checks all eight three-Bernoulli recipe trees, and verifies the scalar grid index sets.

2. **Full interval arithmetic reevaluation of frozen partitions:**
   ```sh
   uv run --frozen python verify/finite/current/check.py --arithmetic
   ```
   Recomputes the outward-rounded Arb interval bounds on every leaf of the frozen partitions in `current-finite.zip`. A single family can be selected with `--case` (for example, `--case margin`, `--case bernoulli_three`, `--case lower110`, `--case low210`, `--case collision221`, `--case collision210`, `--case tt221_negative`, `--case tt221_positive`, `--case floor`, or `--case scalars`).

3. **Full verification including common-diameter curvature and symbolic Bernstein tables:**
   ```sh
   uv run --frozen python verify/finite/check.py --full --independent
   ```
   Runs all partition and arithmetic checks above, reconstructs the 550,563-box common-diameter curvature cover, and independently recomputes the exact $`\mathbb{Q}(\sqrt{5})`$ Bernstein coefficient tables in SymPy.
