# Interval evaluation code

These modules contain the interval arithmetic and evaluators used to
subdivide the feasible moment parameters and bound characteristic functions
and smoothing integrals on each parameter box.
The JSON files supply scalar coefficient tables. The
[cover](../../certificates/cover.json) records the SHA-256 of each source
and table file used by the calculation.

The checking wrappers reconstruct stored partitions or recompute numerical
bounds using these routines. The [verification guide](../README.md)
describes the commands and the scope of each check; the
[scalar guide](../../certificates/scalar/README.md) describes the coefficient
certificates.

The [fixed-input evaluator](../fixed_moment.py) checks the nested moment
proofs stored inside final boxes. It uses the unsigned and signed
smoothing bounds with their corresponding exterior-tail estimates.
The additional routines implement the split Cauchy–Schwarz estimate,
its endpoint reduction, and the convex bound on moment weights stated
in the paper. They preserve the reference center and change only the
forcing bound or error radius. They check their hypotheses and rational
contraction bounds on each parameter box.

For the intersection of two disks, the certificate specifies the
rational complex coefficient $h$ in the paper's support inequality on
each frequency and threshold rectangle. Reevaluation checks the two
disk bounds, the signed division by $L$, and every subdivision before
replacing the corresponding integrated cell contribution.
