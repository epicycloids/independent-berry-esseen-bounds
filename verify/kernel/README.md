# Arithmetic source guide

These nine Python files are the arithmetic sources identified in
[the certificate](../../certificates/cover.json). Their SHA-256 digests are
checked before replay to identify the code used for the saved calculation.

| Source | Mathematical role |
| --- | --- |
| [interval_bounds.py](interval_bounds.py) | Outward-rounded arithmetic, bounds for characteristic functions, and smoothing coefficients |
| [stable_bounds.py](stable_bounds.py) | Moment constraints, logarithmic product bounds, and separate treatment of a largest-variance summand |
| [directional_bounds.py](directional_bounds.py) | An earlier aggregate estimate for the real and imaginary parts of the zero-bias error |
| [refined_integral.py](refined_integral.py) | Integration over frequency intervals, retaining the Gaussian factor on the accumulated integral |
| [even_bounds.py](even_bounds.py) | A one-sided real-part estimate and its earlier aggregate formula |
| [paired_bounds.py](paired_bounds.py) | The polynomial refinement of the aggregate estimate, and a bound for a single characteristic-function modulus |
| [batch_weights.py](batch_weights.py) | Three valid smoothing choices, evaluated together and minimized |
| [gaussian_batch.py](gaussian_batch.py) | Arb integration of a nonsingular upper bound for the Gaussian correction |
| [stable_cover.py](stable_cover.py) | Parameter partition construction and reconstruction |

The source names are implementation names. In functions that evaluate
subarray bounds, `B` denotes the sum of third absolute moments, `tau` denotes
$\sum_jv_j^{3/2}$, and `d` is an upper bound on an individual variance. It
need not be the variance of the separately retained summand. Argument
suffixes `lo` and `hi` denote interval endpoints. A function named `*_cells`
returns upper bounds on frequency intervals, called cells in the code.

## Evaluator used by replay

The public [replay wrapper](../replay.py) calls
`gaussian_batch.install_factory()`. For the certificate's `gaussian` kind,
this selects `GaussianBatchWeights`. Its inheritance chain is

```text
GaussianBatchWeights -> BatchWeights -> CosineWeights -> PairedWeights
                    -> EvenWeights -> RefinedWeights -> DirectionalWeights
                    -> StableWeights -> Weights
```

`CosineWeights.prefactor` calls `paired_bounds.cosine_cells`, which evaluates
the final aggregate formula in the paper. `PairedWeights.single_cells`
combines the original single-summand bound with `paired_cells` by taking
their minimum. `BatchWeights` evaluates three smoothing choices and takes
the smallest result; `GaussianBatchWeights` adds the Arb Gaussian-correction
bound. Earlier classes remain in the inheritance chain, but their
`prefactor` methods are overridden. In particular, neither
`directional_cells` nor `even_cells` supplies the final aggregate formula.

The companion notes define the formulas referred to by source comments:

- [LOG_DEFECT.md](LOG_DEFECT.md): logarithmic bounds for a full product and a product with one factor omitted.
- [DIRECTIONAL_STABILITY.md](DIRECTIONAL_STABILITY.md): the earlier real-and-imaginary aggregate estimate.
- [EVEN_SUPPORT.md](EVEN_SUPPORT.md): the one-sided real-part estimate, the expression called `E3`, and its later polynomial refinement.

The derivations are in [the paper](../../paper/independent-berry-esseen.pdf).
The notes connect its mathematics to the preserved function names.

## Supported commands and diagnostics

Use the commands in [the verification guide](../README.md), which also
states arithmetic assumptions and replay limits. The arithmetic files
contain diagnostic and pilot helpers in addition to the replay evaluator.
The documented commands do not call those helpers.

The `probe` helpers in `even_bounds.py` and `paired_bounds.py` import
`log_probe`, and `refined_integral.integration_checks` imports `global_probe`;
those exploratory modules are not included in this edition. All imports
required by the documented verification commands are included. Sampled
floating-point comparisons in the diagnostic helpers test particular inputs;
the completed calculation instead evaluates outward-rounded bounds over
parameter boxes.
