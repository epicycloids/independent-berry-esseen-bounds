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
