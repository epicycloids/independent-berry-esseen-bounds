# Verification supplement

This supplement contains the computation and checking programs supporting
the proposed bound $C_{\mathrm{ind}}<0.474999998$. Independent mathematical
review remains outstanding. The conjecture that $C_{\mathrm{ind}}$ equals
Esseen’s constant remains open; the [overview](../docs/overview.md) defines
the constant and the normalization.

## Check the saved records

With [uv](https://docs.astral.sh/uv/) installed, run from the edition root:

```sh
./verify.sh --no-latex
```

This command checks the file manifest, runs certificate-rejection tests,
validates the saved bands and tree records, recomputes the analytic endpoint
bounds with Arb, checks the recorded replay summaries against their sources
and inputs, and checks the paper's generated numerical values. It does not
repeat the parameter-box traversal or accepted-leaf quadrature, and launches
no remote work. Dependencies are locked in [uv.lock](../uv.lock);
[.python-version](../.python-version) selects Python 3.12.10.

To build the PDF, install
[Tectonic](https://tectonic-typesetting.github.io/en-US/) (the prepared
edition uses 0.15.0) and run `./build.sh`. The output is
`dist/independent-berry-esseen.pdf`. Running `./verify.sh` without an option
also builds the paper. The included PDF is a reading copy; edit the TeX and
bibliography to revise it.

## Records and their contents

| File | Contents |
| --- | --- |
| [result.json](../result.json) | Candidate status, constants, counts, runtime versions, and the digest of `cover.json` |
| [cover.json](../certificates/cover.json) | All 767 completed moment bands, exact endpoints, partition trees, accepted-leaf checksums, and arithmetic-source digests |
| [topology-replay.json](../certificates/topology-replay.json) | Recorded reconstruction of the full parameter partition |
| [selected-replay.json](../certificates/selected-replay.json) | Recorded numerical recomputation of three bands |
| [kernel/](kernel/README.md) | Nine arithmetic source files used by the evaluator |
| [SHA256SUMS](../SHA256SUMS) | Digests of the edition's other files |

A **moment band** is a closed interval of the normalized third absolute
moment sum $L$. Its stored binary64 endpoints have exact rational values;
the hexadecimal encodings and exact fractions of the full range's endpoints
are included. Printed decimals such as 0.006 and 1.15 are approximations to
those endpoints.

Each band stores its root parameter box, initial smoothing parameters,
acceptance method, and upper bound. The smoothing parameters are a frequency
cutoff $T>0$ and a number $0<s<1$ that splits the integration range at $Ts$.
The stored `T,s` are the first of three smoothing choices for the initial
moment interval. Choices for the subdivided boxes are reconstructed from
the rules in [gaussian_batch.py](kernel/gaussian_batch.py),
[batch_weights.py](kernel/batch_weights.py), and
[refined_integral.py](kernel/refined_integral.py).

A **partitioned band** divides the root box in $(d,b,\tau)$. After
normalizing the total variance to one, write $v_j=\mathbb E X_j^2$; then $d$
is the largest summand variance, $b$ is the third absolute moment of a
summand with that variance, and $\tau=\sum_jv_j^{3/2}$. The compressed tree
uses preorder traversal: `D`, `B`, and `T` split $d$, $b$, and $\tau$,
respectively; `A` accepts a box using an upper bound; `X` excludes a box
using the moment constraints. The tree's SHA-256 digest is computed from
the uncompressed instruction bytes.

Each partitioned band also stores a checksum of its accepted-leaf sequence.
Starting from 32 zero bytes, the calculation updates the digest at each
accepted leaf by hashing the preceding digest followed by
`json.dumps([box, value], separators=(',', ':')).encode()`. Here `box` is the
clipped parameter box and `value` is its computed upper bound. The final
digest, stored as `leaf_chain_sha256`, is a checksum, not a recoverable list
of those values. The individual bounds cannot be read from it. Numerical
replay reconstructs the boxes, recomputes their bounds, rebuilds the checksum,
and compares both the checksum and the largest value with the stored record.

## Scope of each check

| Check | Work performed | Work not repeated |
| --- | --- | --- |
| Saved records | File/source digests, tree syntax and counts, contiguous moment bands, declared upper endpoints, and Arb analytic complements | Parameter-box traversal and accepted-leaf quadrature |
| Geometric replay (`topology`) | Every recorded split and clipped box, every infeasibility decision, counts, and analytic complements | Accepted-leaf quadrature |
| Numerical replay (`selected` or `band`) | Geometric replay plus numerical acceptance checks; partitioned bands must reproduce their leaf checksum and maximum | Bands outside the requested selection |

The included geometric replay covers all 767 bands and 5,577,029 nodes.
The included numerical replay covers these zero-based band indices:

| Index | Printed $L$ interval | Nodes | Accepted | Infeasible |
| --- | --- | ---: | ---: | ---: |
| 0 | 0.006 to 0.00615 | 13 | 7 | 0 |
| 305 | 0.249 to 0.251 | 35 | 16 | 2 |
| 424 | 0.4800000000000001 to 0.4820000000000001 | 10,151 | 4,921 | 155 |
| Total | Three selected bands | 10,199 | 4,944 | 157 |

These bands have lower endpoints nearest 0.006, 0.25, and 0.48. All 4,944
accepted-leaf bounds were recomputed, and each band's checksum and maximum
matched its original record. The original completed calculation evaluated
every accepted-leaf bound in the full cover. The other accepted-leaf
quadratures have not been recomputed in full since that calculation.

A **uniform band** is accepted as a whole, without a parameter partition.
Its numerical replay checks acceptance against 0.475; it has no leaf checksum
and the recomputed value is not required to equal the stored upper endpoint.
A **variance-only band** uses the analytic bound $\Delta<0.54093655$ for
centered variance-one laws. Saved-record checking separately verifies that
every stored upper endpoint lies below the finer candidate value 0.474999998.
Thus uniform-band replay alone does not reproduce that finer endpoint.

## Arithmetic assumptions

The arithmetic source files use 100-bit Arb enclosures for transcendental
coefficients. The analytic-complement checker uses 192 bits. Vector interval
arithmetic assumes IEEE-754 binary64 elementary operations and square root,
rounding to nearest, gradual underflow, and no unsafe fast-math
transformations. Interval endpoints are rounded outward with `nextafter`.
Nonnegative cumulative sums and dot products are multiplied by
$1+4N\,2^{-53}$ and rounded upward, with $N=512$ frequency intervals per
range and smoothing choice. Every cumulative prefix uses the full factor.
For a nonnegative lower exponent $p<25$, the exponential routine evaluates
$P_{16}(p/64)$ downward, where $P_{16}(x)=\sum_{j=0}^{16}x^j/j!$,
takes its reciprocal upward, and performs six upward-rounded squarings.
For $p\ge25$, it returns an upper enclosure of $e^{-25}$.
Section 7 of [the paper](../paper/independent-berry-esseen.pdf) explains
these formulas and identifies their implementations.

Run Python without `-O` and without `PYTHONOPTIMIZE`: the arithmetic code uses
assertions, and the public loader refuses an optimized run. The recorded
library versions are NumPy 2.4.6, python-flint 0.8.0, and SciPy 1.15.3.
Source and input digests are checked before replay.

Independent mathematical review must assess the analytic inequalities,
coverage of the admissible moment domain, and implementation of the
arithmetic assumptions. The checks above establish their stated
computational results; they do not perform that mathematical review.

## Optional replay and compute resources

With an existing authorized Modal account, these commands request cloud
compute and may incur charges:

```sh
uv run --frozen --extra compute modal run verify/modal_replay.py --mode topology
uv run --frozen --extra compute modal run verify/modal_replay.py --mode selected
```

Run them sequentially. Each uses one CPU, at most 1 GiB memory, a 600-second
execution timeout, a 60-second startup timeout, and zero input retries.
Outputs go to `dist/replay/` and do not replace the supplied records. The
app ends after the requested run; no schedule or persistent deployment is
created.

To recompute another chosen band:

```sh
uv run --frozen --extra compute modal run verify/modal_replay.py --mode band --band-index 424
```

This mode has the same CPU, memory, startup, and retry limits, with a
3,600-second execution timeout. A full numerical replay requires selecting
every band and retaining every result; this edition has no automatic
all-band dispatch. Runtimes vary substantially by band. Estimate compute
from measured runtimes and tree sizes rather than extrapolating from the
three selected bands alone.

The evaluator can also run locally:

```sh
uv run --frozen python verify/replay.py --mode topology --output dist/topology.json
uv run --frozen python verify/replay.py --mode band --band-index 424 --output dist/band-424.json
```

Local replay can be intensive. These scripts impose no operating-system
resource cap and are separate from `./verify.sh`.

## Source map

The [arithmetic guide](kernel/README.md) identifies the formulas and the
classes used by the replay evaluator. The public wrappers validate records
and select bands; the nine arithmetic files evaluate the bounds.

The paper's numerical values are generated from `result.json` by
`paper_values.py`. After an intentional result change, regenerate them with
`--write`, review the manuscript, and refresh the edition manifest when
preparing the new version.
