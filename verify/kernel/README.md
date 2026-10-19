# Interval Evaluation Kernel

Navigation: [Repository Root](../../README.md) · [Verification Guide](../README.md) · [Scalar Certificates](../../certificates/scalar/README.md)

---

The modules in `verify/kernel/` implement the outward-rounded IEEE-754 binary64 interval arithmetic and Arb ball-arithmetic routines used to prove the upper bound $C_{\mathrm{ind}} < 0.44988$ in [Manuscript 1](../../paper/independent-berry-esseen.pdf). The SHA-256 digest of every kernel source file and coefficient table is pinned in [`certificates/cover.json`](../../certificates/cover.json).

## Principal Components

1. **Characteristic-function product bounds (`interval_bounds.py`, `stable_bounds.py`, `moment_log.py`, `higher_log_bounds.py`, `dual_log_bounds.py`, `weighted_log_bounds.py`):**
   Bound the modulus of $\prod_j f_{X_j}(t)$ and of products with one factor omitted over boxes of the moment parameters $(L, d, b, \tau)$, combining Prawitz–Tyurin cubic/cosine bounds with quadratic and tangent-line lower bounds for $-\frac{1}{2}\log(1 - 2q_j)$.
2. **Zero-bias error bounds (`directional_bounds.py`, `even_bounds.py`, `paired_bounds.py`, `quadratic_prefactor.py`, `concave_bounds.py`, `scalar_table_bounds.py`, `disk_table_bounds.py`):**
   Bound $\sum_i v_i |f_{X_i}(t) - f_{X_i^*}(t)|$ using the third-moment excess $\rho_i - 1$, least concave majorants, and the certified affine modulus and directional tables in `scalar_table.json` and `disk_table.json`.
3. **Cosine-product reference comparison (`reference_bounds.py`, `energy_reference.py`, `energy_signed.py`, `packed_energy_reference.py`, `packed_energy_signed.py`, `cauchy_reference_adapter.py`, `cauchy_endpoint_adapter.py`, `packed_weight_adapter.py`, `packed_weight_general_adapter.py`):**
   Bound the error $|f_S(t) - \prod_i \cos(t\sqrt{v_i})|$ via the first-order differential recurrence in Section 6 of Manuscript 1 and the scalar majorant table `energy_table_complete.json`, while enclosing the signed displacement $\prod_i \cos(t\sqrt{v_i}) - e^{-t^2/2} \le 0$.
4. **Prawitz smoothing and threshold covers (`refined_integral.py`, `interval_signed.py`, `signed_gaussian.py`, `gaussian_batch.py`, `cantelli_gaussian_tail.py`, `full_b_consumer.py`):**
   Integrate the unsigned and signed Prawitz smoothing bounds over frequency grids, combine them with the Arb-enclosed Gaussian defect tables, and evaluate threshold covers on $[-2, 2]$ and $[-69/32, 69/32]$ together with Cantelli exterior-tail bounds.
5. **Moment-box and nested-proof evaluators (`stable_cover.py`, `full_b_cover.py`, `adaptive_b_consumer.py`, `../fixed_moment.py`, `../moment_proofs.py`):**
   Traverse the outer $(d, b, \tau)$ partition trees and nested $b$- or $d$-subdivision proofs, verifying moment-feasibility exclusions, two-disk complex support inequalities, and outward-rounded cell sums.
