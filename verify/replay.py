"""Recompute accepted-leaf integrals in one or three stored L intervals.

Uses each record's evaluator, parameters and acceptance target. Large
intervals can take substantial time. This command does not independently
verify every Gaussian quadrature node or review the analytic inequalities."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

from check_saved import ROOT, canonical, load_inputs, mixed_record, require


def install(kind):
    if kind == "quadratic_maxwell_dual5_cutoff":
        import quadratic_cover_engine as engine
    elif kind == "quadratic_maxwell_dual5_cutoff_lazy_ceil20":
        import cap_portfolio_engine as engine
    elif kind.startswith("full_b_fast_maximum_fixed085_dual5_cutoff_v1"):
        import full_b_maximum_cover_engine as engine
    elif kind == "adaptive_b_fast_maximum_fixed085_dual5_cutoff_v1":
        import adaptive_b_cover_engine as engine
    else:
        raise ValueError("Unrecognized evaluator family: " + kind)
    engine.install()
    return engine


class MixedCallbacks:
    """Use a separate evaluator for each of the two recorded acceptance targets."""
    def __init__(self, record):
        self.record = record
        self.engine = install(record["evaluator"]["kind"])
        self.factories = {}
        self.calls = {"origin":0,"continuation":0}

    def score(self, role, descriptor, lo, hi, box):
        from flint import ctx
        expected = self.record["origin_evaluator"] if role == "origin" else self.record["evaluator"]
        require(canonical(descriptor) == canonical(expected), "Evaluator role/configuration changed")
        require((lo, hi) == (self.record["Llo"], self.record["Lhi"]), "Evaluator domain changed")
        with ctx.workprec(100):
            if role not in self.factories:
                choose = self.engine.stable_cover.weight_factory(hi, descriptor["N"], descriptor["kind"], descriptor["target"])
                weight = choose(hi)
                require((weight.T, weight.s) == (descriptor["T"], descriptor["s"]), "Smoothing parameters changed")
                require(weight.rectangle_class.__name__ == "CachedMaximumVarianceRectangles" and
                        weight.coupling is False, "Wrong adaptive evaluator")
                self.factories[role] = choose
            value = self.factories[role](box[-1]).box(lo, hi, *box)
            self.calls[role] += 1
            return value

    def origin(self, descriptor, lo, hi, box):
        return self.score("origin",descriptor,lo,hi,box)

    def continuation(self, descriptor, lo, hi, box):
        return self.score("continuation",descriptor,lo,hi,box)


def replay_band(index):
    from flint import ctx
    cover, result = load_inputs()
    require(type(index) is int and 0 <= index < len(cover["bands"]), "Band index is outside the cover")
    band = cover["bands"][index]
    record = band["record"]
    with ctx.workprec(100):
        if band["type"] == "mixed":
            callbacks = MixedCallbacks(record)
            counts = mixed_record(record, numerical=True, callbacks=(callbacks.origin,callbacks.continuation))
            roles = callbacks.calls
        else:
            import stable_cover
            install(record["kind"])
            counts = stable_cover.replay_band(record, numerical=True)
            roles = None
    return {"band_index":index, "accepted_leaf_quadrature_replayed":True, "counts":counts,
            "source_separated_numerical_calls":roles, "global_proof":False,
            "whole_union_numerically_replayed":False, "gaussian_node_audit_completed":False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--band-index",type=int,help="recompute accepted-leaf integrals for one stored L-interval index")
    group.add_argument("--selected",action="store_true",help="recompute the three selected L intervals, serially")
    args = parser.parse_args()
    if args.selected:
        _, result = load_inputs()
        # Isolate evaluator factories and caches between bands.
        for index in result["selected_replay_band_indices"]:
            subprocess.run([sys.executable,str(Path(__file__).resolve()),"--band-index",str(index)],check=True)
    else:
        print(json.dumps(replay_band(args.band_index),indent=2))
