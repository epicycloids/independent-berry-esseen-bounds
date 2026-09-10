"""Reconstruct the closed cover from stored partitions and numerical bounds.

Checks source hashes, Gaussian interpolation partitions and exact comparisons.
The accepted-leaf and Gaussian-node integrals are not recomputed."""
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

from scalars import scalar_bounds

ROOT = Path(__file__).resolve().parent.parent
KERNEL = Path(__file__).resolve().parent / "kernel"
sys.path.insert(0, str(KERNEL))


class CertificateError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise CertificateError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def finite(value):
    require(type(value) in (int, float) and math.isfinite(value), "Expected a finite number")
    return Fraction(value)


def load_inputs(root=ROOT):
    require(__debug__, "Run without Python optimization: the curated kernels use assertions")
    result = json.loads((root / "result.json").read_text())
    raw = (root / "certificates/cover.json").read_bytes()
    require(sha256(raw) == result["cover_sha256"], "Cover file digest mismatch")
    cover = json.loads(raw)
    for name, expected in cover["kernel_sha256"].items():
        require(Path(name).name == name and Path(name).suffix in (".py", ".json"), "Invalid kernel filename")
        require(sha256((root / "verify/kernel" / name).read_bytes()) == expected,
                "Curated mathematical source changed: " + name)
    return cover, result


def mixed_record(record, *, numerical=False, callbacks=None):
    """Verify that the second tree covers the first tree's unresolved boxes."""
    from flint import ctx
    import prefix_reuse as geometry
    require(record["kind"] == "dt_verified_prefix_continuation_v1" and
            record["method"] == "mixed_source_prefix", "Wrong mixed record kind")
    require(record["geometry_precision_bits"] == 100, "Geometry precision changed")
    require(record["dt_priority"] == [1., 1.], "Continuation subdivision policy changed")
    lo, hi = record["Llo"], record["Lhi"]
    origin, continuation = record["origin"], record["continuation"]
    old_descriptor, new_descriptor = record["origin_evaluator"], record["evaluator"]
    require(record["complete"] is True and not record["stack"] and record["unresolved"] == 0, "Mixed tree is unfinished")
    require(origin["complete"] is False and origin["split_axes"] == "DT" and
            origin["b_policy"] == "full_feasible_each_evaluation", "Wrong origin contract")
    require((origin["Llo"], origin["Lhi"]) == (lo, hi), "Origin domain differs")
    require(origin["target"] == old_descriptor["target"] == "0.453" and
            record["target"] == new_descriptor["target"] == "0.454", "Source-specific targets changed")
    expected_parameters = dict(theta_denominator=32, threshold_domain=[-2, 2], coupling=False,
        b_max_depth=6, b_max_evaluations=127, b_depth_policy="depth_for_box",
        fixed_cap_fraction=[17, 20], reference_cache_entries=16,
        cutoff_portfolio="frozen fast_maximum_cover_engine.CUTOFFS")
    for descriptor in (old_descriptor, new_descriptor):
        require(descriptor["kind"] == "adaptive_b_fast_maximum_fixed085_dual5_cutoff_v1" and
                descriptor["N"] == 512 and descriptor["weight_parameter"] == hi and
                descriptor["parameters"] == expected_parameters, "Mixed evaluator parameters changed")
        require(0 < finite(descriptor["s"]) < 1 and finite(descriptor["T"]) > 0, "Invalid smoothing parameters")
    for field in ("kind", "N", "T", "s", "weight_parameter"):
        require(old_descriptor[field] == origin[field] == new_descriptor[field], "Origin evaluator descriptor mismatch")
    require(new_descriptor["dt_priority"] == [1., 1.], "New descriptor priority mismatch")
    require(not numerical or callbacks is not None, "Numerical replay needs two explicit evaluator callbacks")
    with ctx.workprec(100):
        root = tuple(geometry._canonical_root(hi))
        require(list(root) == record["root"] == origin["root"], "Mixed canonical root differs")
        old_tree = geometry._unpack(origin)
        old = geometry._walk(old_tree, [(root, 0)], hi,
            evaluate=callbacks[0] if numerical else None, descriptor=old_descriptor,
            lo=lo, target=Fraction(origin["target"]))
        frontier = geometry._stack(origin["stack"])
        require(old["stack"] == frontier and len(frontier) == origin["unresolved"] > 0, "Origin frontier mismatch")
        new_tree = geometry._unpack(continuation)
        new = geometry._walk(new_tree, frontier, hi,
            evaluate=callbacks[1] if numerical else None, descriptor=new_descriptor,
            lo=lo, target=Fraction(record["target"]))
        require(not new["stack"], "Continuation does not finish the original frontier")
    for saved, replay, target in ((origin, old, "0.453"), (continuation, new, "0.454")):
        require(0 <= finite(saved["upper"]) < Fraction(target), "Saved component upper exceeds its own target")
        require(len(bytes.fromhex(saved["leaf_chain_sha256"])) == 32, "Invalid leaf chain")
        for key in ("accepted", "infeasible"):
            require(type(saved[key]) is int and saved[key] == replay[key], "Mixed leaf count mismatch")
        if numerical:
            require(replay["chain"] == saved["leaf_chain_sha256"] and replay["upper"] == saved["upper"] and
                    canonical(replay["worst_box"]) == canonical(saved["worst_box"]), "Numerical leaf chain/maximum differs")
    require(record["upper"] == max(origin["upper"], continuation["upper"]), "Mixed maximum differs")
    for key in ("nodes", "accepted", "infeasible"):
        require(type(record[key]) is int and record[key] == origin[key] + continuation[key], "Mixed aggregate counts differ")
    return {key: record[key] for key in ("nodes", "accepted", "infeasible")}


def validate_cover(cover, result):
    from flint import ctx
    import stable_cover
    require(cover["schema"] == "independent-be-cover-v2" and result["schema"] == "independent-be-result-v2", "Unsupported certificate schema")
    require(result["status"] == "computational-candidate", "Unexpected mathematical status")
    require(cover["global_proof"] is result["global_proof"] is False and
            result["independent_mathematical_review_completed"] is False and
            result["sharp_conjecture_resolved"] is False and
            result["all_leaf_quadratures_replayed"] is False, "Unsupported proof or replay promotion")
    require(cover["comparison_target"] == result["comparison_target"] == result["upper_bound"] == "0.454", "Wrong current target")
    domain = cover["domain"]
    for side in ("lower", "upper"):
        require(finite(domain[side]) == Fraction(domain[side + "_fraction"]), "Exact endpoint mismatch")
        require(float.fromhex(domain[side + "_hex"]) == domain[side], "Dyadic endpoint mismatch")
    require(domain["lower"] == .0014 and domain["upper"] == 1.21, "Wrong candidate domain")
    require(bool(cover["bands"]), "Empty finite cover")
    counts = Counter(bands=0, nodes=0, accepted=0, infeasible=0)
    methods, families = Counter(), Counter()
    intervals, uppers = [], []
    with ctx.workprec(100):
        for band in cover["bands"]:
            record = band["record"]
            lo, hi, upper = (finite(record[k]) for k in ("Llo", "Lhi", "upper"))
            require(0 < lo < hi and Fraction(domain["lower"]) <= lo < hi <= Fraction(domain["upper"]), "Band leaves finite domain")
            require(record["complete"] is True and 0 <= upper < Fraction(record["target"]) <= Fraction("0.454"), "Incomplete or nonaccepting band")
            if band["type"] == "standard":
                require(band["weights_regenerated_on_numerical_replay"] is True, "Legacy weight scope differs")
                require(record["method"] in ("stability", "uniform", "variance"), "Unknown band method")
                require(record["kind"] in ("quadratic_maxwell_dual5_cutoff",
                    "quadratic_maxwell_dual5_cutoff_lazy_ceil20",
                    "full_b_fast_maximum_fixed085_dual5_cutoff_v1_b8",
                    "adaptive_b_fast_maximum_fixed085_dual5_cutoff_v1"), "Unsupported evaluator family")
                require(record["weight_parameter"] == record["Lhi"] and record["N"] in (512,1024,2048), "Unexpected evaluator grid")
                local = stable_cover.replay_band(record, numerical=False)
                for key, value in local.items():
                    require(type(value) is int and value >= 0 and (key not in record or record[key] == value), "Standard tree counts disagree")
            else:
                require(band["type"] == "mixed", "Unknown band type")
                local = mixed_record(record)
            counts.update(local)
            counts["bands"] += 1
            methods[record["method"]] += 1
            families[record.get("evaluator",record)["kind"]] += 1
            intervals.append((lo, hi))
            uppers.append(upper)
    cursor = Fraction(domain["lower"])
    for lo, hi in sorted(intervals):
        require(lo <= cursor, "Finite closed union has a gap")
        cursor = max(cursor, hi)
    require(cursor == Fraction(domain["upper"]), "Finite union misses the endpoint")
    require(dict(counts) == result["counts"] == {"bands":1683,"nodes":412067,"accepted":200239,"infeasible":6636}, "Aggregate counts disagree")
    require(dict(methods) == result["band_methods"] and dict(families) == result["evaluator_families"], "Method counts disagree")
    require(max(uppers) == Fraction(result["largest_recorded_upper_fraction"]) == Fraction(result["largest_recorded_upper"]) ==
            Fraction(4089268438506339,9007199254740992), "Maximum upper enclosure differs")
    return {"counts":dict(counts), "largest_recorded_upper_fraction":str(max(uppers)),
            "closed_domain":[str(Fraction(domain["lower"])),str(cursor)]}


def validate_gaussian(data, cover):
    import signed_gaussian
    require(data["schema"] == "independent-be-gaussian-v2", "Unknown Gaussian schema")
    counts = Counter(records=0, threshold_nodes=0, intervals=0)
    for key, scalar in data["records"].items():
        require(key == sha256(canonical(scalar)), "Scalar record digest differs")
        checked = signed_gaussian.verify_partition(scalar)
        require(checked["all_passed"] is True, "Gaussian interpolation failed")
        counts["records"] += 1
        counts["threshold_nodes"] += checked["threshold_nodes"]
        counts["intervals"] += checked["intervals"]
    for key, weights in data["weight_sets"].items():
        require(key == sha256(canonical(weights)), "Weight-set digest differs")
        require(set(weights["scalar_ids"]) <= data["records"].keys(), "Missing scalar record")
        full_splits = {data["records"][c["scalar_id"]]["T"]: c["rounded_split"]
                       for c in weights["cutoffs"] if c["k"] == c["N"]}
        for cutoff in weights["cutoffs"]:
            scalar = data["records"][cutoff["scalar_id"]]
            require(cutoff["N"] > 0 and 1 <= cutoff["k"] <= cutoff["N"], "Invalid Gaussian cutoff")
            require(cutoff["rounded_split"] == scalar["s"] and finite(cutoff["upper"]) >= 0 and
                    finite(cutoff["endpoint_correction"]) >= 0, "Scalar cutoff mismatch")
            if cutoff["k"] == cutoff["N"]:
                # The full split also takes minima with two other Gaussian
                # bounds. Their numerical values are regenerated on replay.
                require(cutoff["endpoint_correction"] == 0 and cutoff["upper"] <= scalar["upper"],
                        "Full-split comparison metadata changed")
            else:
                from flint import arb, ctx
                from interval_bounds import two_kernel, au
                require(cutoff["exact_split"] == "s*k/N", "Wrong partial cutoff convention")
                with ctx.workprec(100):
                    exact = arb(full_splits[scalar["T"]]) * cutoff["k"] / cutoff["N"]
                    approximate = arb(cutoff["rounded_split"])
                    correction = au(abs(exact-approximate)*two_kernel(exact.union(approximate)))
                require(cutoff["endpoint_correction"] >= correction, "Cutoff endpoint correction is understated")
                require(finite(cutoff["upper"]) >= finite(scalar["upper"]) + finite(cutoff["endpoint_correction"]),
                        "Cutoff bound understates its correction")
    for band in cover["bands"]:
        for name in (() if band["type"] == "standard" else ("origin_weights","continuation_weights")):
            require(band[name] in data["weight_sets"], "Missing band weights")
    return dict(counts)


def check(root=ROOT):
    cover, result = load_inputs(root)
    summary = validate_cover(cover, result)
    raw_gaussian = (root / cover["gaussian_file"]).read_bytes()
    require(sha256(raw_gaussian) == result["gaussian_file_sha256"], "Mixed Gaussian data digest differs")
    gaussian = json.loads(raw_gaussian)
    summary["stored_gaussian_partitions"] = validate_gaussian(gaussian, cover)
    summary["scalar_enclosures"] = scalar_bounds(cover["domain"]["lower"], cover["domain"]["upper"])
    summary.update(status="exact saved-data checks passed", global_proof=False,
        tree_geometry_replayed=True, accepted_leaf_quadrature_replayed=False,
        legacy_gaussian_partitions_replayed=False,
        gaussian_node_quadrature_replayed=False, independent_mathematical_review_completed=False)
    return summary


if __name__ == "__main__":
    try:
        print(json.dumps(check(), indent=2))
    except (CertificateError, AssertionError, ArithmeticError, KeyError, ValueError) as error:
        print(f"Certificate check failed: {error}", file=sys.stderr)
        raise SystemExit(1)
