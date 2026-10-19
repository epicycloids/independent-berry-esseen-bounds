# Scalar Certificate Verification (Manuscript 1)

Navigation: [Repository Root](../../README.md) · [Verification Guide](../../verify/README.md) · [Interval Kernel](../../verify/kernel/README.md)

---

This directory contains the scalar majorant certificates and generator scripts used in [Manuscript 1](../../paper/independent-berry-esseen.pdf) ([`paper/main.tex`](../../paper/main.tex)) to bound single-summand characteristic functions, logarithmic product corrections, zero-bias differences, and the cosine-comparison forcing term.

## Running the Checks

From the repository root, run:

```sh
uv run --frozen python verify/check_scalar_certificates.py
```

This default check (which runs automatically as part of `./verify.sh --no-latex`) performs the following verifications:

1. **Quadratic imaginary-error polynomial (`odd-quadratic.json`, `source/odd_polynomial.py`):** Regenerates the exact rational polynomials in Appendix A of Manuscript 1, converts them to the tensor Bernstein basis on $`[0, 4] \times [0, 49]`$, verifies positivity of the minimum Bernstein coefficients minus the series tail, and checks the spatial tail for $`z \ge 7`$.
2. **Logarithmic tangent lower bounds (`direct-log.json`, `source/direct_log.py`):** Uses 192-bit Arb ball arithmetic to verify the convexity/concavity transition and endpoint inequality of Lemma 2.3 for all 2,006 tangent lower bounds across $`\lambda \in \{3/4, 1, 3/2, 2, 3\}`$.
3. **Affine zero-bias modulus and directional bounds (`affine-index.json`, `disk.zip`, `directional.zip`):** Verifies the binary partition trees covering all frequency cells and spatial intervals, checks the exact analytic tail polynomials outside $`[-R, R]`$, verifies the 48-direction polygon inradius $`\eta = 997/1000`$, and confirms that the reconstructed cumulative maxima match `verify/kernel/scalar_table.json` and `verify/kernel/disk_table.json`.
4. **Cosine-comparison forcing bound (`energy.json`, `source/energy_majorant.py`, `source/energy_small.py`):** Verifies the rectangle partitions covering $`(x, y) \in [0, 1/2] \times [0, 48]`$ and $`[1/2, 1] \times [0, 16]`$, checks the exact positive tail polynomials for $`y \ge 48`$ and $`y \ge 16`$, and confirms exact agreement with `verify/kernel/energy_table_complete.json`.

### Recomputing Leaf Enclosures with Arb

The default check verifies partition completeness, exact rational tail conditions, and table reconstruction from the stored leaf enclosures. To recompute every leaf enclosure from scratch using Arb ball arithmetic:

```sh
uv run --frozen python verify/check_scalar_certificates.py --replay directional --replay disk --replay energy
```

Adding `--limit N` restricts reevaluation to the first `N` certificates in each family for a faster spot check.

---

## Files and Mathematical Roles

| File | Mathematical content |
| --- | --- |
| `odd-quadratic.json`, `source/odd_polynomial.py` | Exact rational polynomial generator and Bernstein coefficient bounds for Lemma 5.1 and Appendix A (quadratic bound on $`\mathrm{Im}(f_W - f_{W^*})`$). |
| `direct-log.json`, `source/direct_log.py` | 2,006 certified tangent lower bounds for $`-\frac{1}{2}\log(1 - 2q) - q`$ across $`\lambda \in \{3/4, 1, 3/2, 2, 3\}`$ (Lemmas 2.2–2.3). |
| `affine-index.json`, `disk.zip`, `source/disk.py` | 640 modulus certificates (75,046 leaves) across 128 frequency cells on $`[0, 2]`$ and common slopes $`\lambda \in \{5/4, 3/2, 2, 3, 4\}`$ (Section 5.2). |
| `affine-index.json`, `directional.zip`, `source/directional.py` | 7,625 directional projection certificates (1,558,869 leaves) across 61 frequency cells, 5 slopes, and 25 upper-half-plane directions (reflected by conjugation to 48 directions). |
| `energy.json`, `source/energy_majorant.py`, `source/energy_small.py` | 128 frequency-cell certificates (25,175 rectangles on $`[0, 1/2]`$ and 150,453 rectangles on $`[1/2, 1]`$) for the cosine-comparison forcing inequality in Section 6.1 and Appendix B. |
| `source/entire_target.py` | Analytic formulas and Taylor enclosures for the normalized zero-bias difference $`2(f_W(x) - f_{W^*}(x))/x`$. |

## Record Format

All rational numbers in the JSON archives are stored as exact fraction strings, and runtime lookup tables store exact IEEE-754 binary64 floats.

- **Affine certificates (`directional.zip`, `disk.zip`):** Each entry records its polynomial coefficients $`(\beta, \gamma, c_0, c_2, \Lambda)`$, base correction $`\delta`$, closed frequency interval $`[x_0, x_1]`$, working precision, spatial cutoff $`R`$, tail polynomial coefficients, and binary partition tree over $`[x_0, x_1] \times [-R, R]`$ (or $`[0, R]`$ for disk certificates). Split nodes are encoded as `[axis, cut, left_index, right_index]`, where `axis` is `0` for frequency $`x`$ and `1` for spatial variable $`w`$.
- **Cosine-comparison certificates (`energy.json`):** Each frequency cell lists closed rectangles `[xlo, xhi, ylo, yhi, lower]` whose union covers $`[x_0, x_1] \times [0, Y]`$ (with $`Y = 48`$ on $`[0, 1/2]`$ and $`Y = 16`$ on $`[1/2, 1]`$) with positive lower bounds on the desingularized quotient $`\mathcal{D}_x(y) / [x^4(y - 1)^4]`$ or $`\mathcal{D}_x(y) / (y - 1)^4`$.
