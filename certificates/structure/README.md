# Bernoulli Structure Certificates (Manuscript 5)

Navigation: [Repository Root](../../README.md) · [Verification Guide](../../verify/README.md) · [Finite-Array Certificates](../finite/README.md)

---

This directory contains the interval certificates and rational benchmark records supporting [Manuscript 5: *Bernoulli extremizers for the independent Berry–Esseen ratio*](../../paper/bernoulli-structure/bernoulli-structure.pdf) ([`paper/bernoulli-structure/main.tex`](../../paper/bernoulli-structure/main.tex)).

## Summary of Certified Calculations

| Theorem / Lemma in Manuscript 5 | Mathematical domain and inequality | Certificate files |
| --- | --- | --- |
| **Theorem 1.1 & Lemma 2.2** (Common-diameter maximizer classification) | Positivity of the one-exception fiber curvature on $V \ge 1/5$, $pq \le W \le 3$ at coefficient $c_* = 199/500$ (550,563 processed boxes, 275,283 accepted boxes at 192 bits), and bound $\mathcal{C} < 199/500$ for two unit-diameter Bernoulli summands (13,641 processed boxes, 6,822 accepted boxes). | `common_gap.json`, `common_gap-input.json`, `benchmarks.json` |
| **Theorem 1.1** (Section 2.5: Two common-diameter summands at threshold $37/100$) | Exclusion of interior stationary points with $\mathcal{C} \ge 37/100$ on the off-diagonal probability square (4,683 processed boxes) and on the diagonal $p_1 = p_2$ (195 processed boxes). | `two_equal_gaps.json`, `two_equal_gaps-input.json` |
| **Theorem 1.2** (Section 3: Two Bernoulli summands with arbitrary diameters) | Exclusion of open stationary points with $0 < t < 1$, $0 < p, q < 1$, $V > 1369/3319684$, and $\mathcal{C} \ge 37/100$ across all three interior threshold branches (37,287 source leaves and 37,421 reevaluated descendant leaves). | `two_arbitrary_gaps-base.json`, `two_arbitrary_gaps-continuation.json`, `two_arbitrary_gaps-replay.json`, `two_arbitrary_gaps-input.json` |

In addition, `benchmarks.json` records the rational lower bounds for the explicit three-summand ($p = 43/100$, $k = 1$, ratio $> 199/500$) and five-summand ($p = 9/20$, $k = 2$, ratio $> 2/5$) binomial examples used in Section 2.4, together with the one-summand example ($p = 2/5$, ratio $> 37/100$) used in Section 2.5.

---

## Running the Verifier

1. **Fast partition and benchmark check (default in `./verify.sh --no-latex`):**
   ```sh
   uv run --frozen python verify/structure/check.py
   ```
   Verifies file hashes, exact rational benchmark inequalities, formula controls, and the complete partition trees of the arbitrary-diameter forests (including all 37,421 reevaluation subpartition leaves).

2. **Full interval arithmetic reevaluation:**
   ```sh
   uv run --frozen python verify/structure/check.py --arithmetic
   ```
   Reconstructs the common-diameter and equal-diameter interval subdivisions from scratch and reevaluates all 37,421 arbitrary-diameter descendant leaves using 192-bit Arb ball arithmetic. Individual families can be selected with `--case common-gap`, `--case two-equal-gaps`, or `--case two-arbitrary-gaps`.
