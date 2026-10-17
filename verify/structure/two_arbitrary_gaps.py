





from fractions import Fraction as Q
import base64
import hashlib
import json
import math
from pathlib import Path
import time
import zlib
from flint import arb, ctx

BASE_SOURCE_PIN = "337c26011b32e283fce2d58876bf2afab8f10305fb5431b340e70ce6bc665d9a"
CONTINUATION_SOURCE_PIN = "996ef8abdf9913ed951fca4204eb96eef2cf8991d596a87a256e93b1a4682cdd"
BASE_CERTIFICATE_PIN = "462730ebde070589ce3b4b2b0571e6c89c952fb57730c7c0a43bc00f9ae0bdf7"
CONTINUATION_CERTIFICATE_PIN = "6500c112966f07cd2ea4de83f8d0f7c4824f4e99de1dff001f3062748fefe223"
BASE_TRACE_PIN = "847ef4e0891dfde530d446af78cac379bb999403550c655eb667aaeaa0c7b2ce"
CONTINUATION_TRACE_PIN = "701bdfec369b0ec9271337664f131f382838c19f981aba74b682c5637454a2f4"
REASONS = ["outside_conditional_variance_domain", "outside_conditional_resource_domain",
           "ratio_below_37_over_100", "positive_branch_Holder_below_37_over_100",
           "factored_gap_numerator_strict_positive", "factored_gap_numerator_strict_negative",
           "p_gradient_score_strict_positive", "p_gradient_score_strict_negative",
           "q_gradient_score_strict_positive", "q_gradient_score_strict_negative",
           "fully_factored_gap_strict_positive", "fully_factored_gap_strict_negative"]
VARIANCE_FLOOR = Q(1, 2500)
COMPETITIVE = Q(37, 100)
BASE64_LIMIT = 1480000


def rational(x):
    return arb(x.numerator) / x.denominator


def ball(l, u):
    return arb(str((l + u) / 2), str((u - l) / 2))


def split(box, axis):
    a, b = box[axis]
    mid = (a + b) / 2
    result = []
    for bounds in ((a, mid), (mid, b)):
        child = list(box)
        child[axis] = bounds
        result.append(tuple(child))
    return result


def decode(path, initial=None):
    box = ((Q(0), Q(1)),) * 3 if initial is None else initial
    for digit in path:
        if digit not in "012345":
            raise ValueError("invalid exact path")
        axis, direction = divmod(int(digit), 2)
        box = split(box, axis)[direction]
    return box


def forest(paths):
    ordered = sorted(paths)
    for a, b in zip(ordered, ordered[1:]):
        if b.startswith(a):
            raise ValueError("overlapping or duplicate paths")
    active = set(paths)
    for depth in range(max(map(len, paths), default=0), 0, -1):
        for path in list(active):
            if len(path) != depth or int(path[-1]) % 2:
                continue
            sibling = path[:-1] + str(int(path[-1]) + 1)
            if sibling in active:
                parent = path[:-1]
                if parent in active:
                    raise ValueError("duplicate cube produced by different split axes")
                active.remove(path)
                active.remove(sibling)
                active.add(parent)
    if active != {""}:
        raise ValueError("exact cube forest incomplete")


def var_range(l, u):
    low = min(l * (1 - l), u * (1 - u))
    high = Q(1, 4) if l <= Q(1, 2) <= u else max(l * (1 - l), u * (1 - u))
    return low, high


def independent_bounds(box, owner):
    (pl, pu), (ql, qu), (tl, tu) = box
    vlp, vup = var_range(pl, pu)
    vlq, vuq = var_range(ql, qu)
    vlo, vhi = vlp + tl ** 2 * vlq, vup + tu ** 2 * vuq
    if vhi < VARIANCE_FLOOR:
        return {"outside_conditional_variance_domain": str(vhi)}
    vlo = max(vlo, VARIANCE_FLOOR)
                                                              
    budget_low = max(vlp - 2 * vlp * vlp + tl ** 3 * (vlq - 2 * vlq * vlq), vlo / 75)
    budget_high = vup - 2 * vup * vup + tu ** 3 * (vuq - 2 * vuq * vuq)
    if budget_high < budget_low:
        return {"outside_conditional_resource_domain": str(budget_high)}
    p, q, t = ball(pl, pu), ball(ql, qu), ball(tl, tu)
    variance, budget = ball(vlo, vhi), ball(budget_low, budget_high)
    deviation = (-p - t * q if owner == 0 else
                 (t * (1 - q) - p if owner == 1 else 1 - p - t * q))
    cumulative = ((1 - p) * (1 - q) if owner == 0 else
                  (1 - p if owner == 1 else 1 - p * q))
    rootvar = variance.sqrt()
    standardized = deviation / rootvar
    gaussian = (-standardized / arb(2).sqrt()).erfc() / 2
    density = (-deviation * deviation / (2 * variance)).exp() / (2 * arb.pi()).sqrt()
    distance = cumulative - gaussian
    penalty = 3 * distance + standardized * density
    kp, kq = 1 - 2 * p, 1 - 2 * q
    dfp = -(1 - q) if owner == 0 else (arb(-1) if owner == 1 else -q)
    dfq = -(1 - p) if owner == 0 else (arb(0) if owner == 1 else -p)
                                                                                  
    gp = dfp + density / rootvar + kp * penalty / (2 * variance) - kp ** 3 * distance / budget
    gq = dfq + t * density / rootvar + t ** 2 * kq * penalty / (2 * variance) - t ** 3 * kq ** 3 * distance / budget
    aq = ball(vlq, vuq)
    rq = 1 - 2 * aq
    weight, sign = (q, -1) if owner == 1 else (1 - q, 1)
    raw_gap = (sign * density / rootvar + weight * t * standardized * density / variance
               + 3 * distance * weight * t / variance - 3 * distance * weight * t ** 2 * rq / budget)
    rp = 1 - 2 * ball(vlp, vup)
    ap = ball(vlp, vup)
                                                                
                                        
                                                    
                                                                     
    first = (p * ((1 - p) - t * (1 - q)) if owner == 0 else
             (-p * ((1 - p) + t * q) if owner == 1 else
              (1 - p) * (p + t * (1 - q))))
    gap = density * budget * first + 3 * distance * weight * t * ap * rootvar * (rp - t * rq)
    if owner == 0:
        full = (density * budget * ((1 - p) - t * (1 - q))
                + 3 * distance * (1 - q) * t * (1 - p) * rootvar * (rp - t * rq))
    elif owner == 2:
        full = (density * budget * (p + t * (1 - q))
                + 3 * distance * (1 - q) * t * p * rootvar * (rp - t * rq))
    else:
        full = None
    reciprocal_shape_upper = vup / (1 - 2 * vup) ** 2 + vuq / (1 - 2 * vuq) ** 2
    return {"ratio_below_37_over_100": distance * variance * rootvar / budget,
            "positive_branch_Holder_below_37_over_100": distance * rational(reciprocal_shape_upper).sqrt(),
            "p_gradient_score": gp, "q_gradient_score": gq,
            "factored_gap_numerator": gap, "fully_factored_gap": full,
            "unmultiplied_gap": raw_gap, "positive_gap_multiplier": variance * budget * rootvar}


def check_reason(ex, owner, reason):
    if reason not in REASONS or (reason in REASONS[10:] and owner not in (0, 2)):
        raise ValueError("unknown or inapplicable pinned producer leaf reason")
    if reason in ("outside_conditional_variance_domain", "outside_conditional_resource_domain"):
        return (True, ex[reason]) if reason in ex else (False, None)
    if any(key.startswith("outside_conditional_") for key in ex):
                                                                           
                                                                         
        key = next(iter(ex))
        return True, "empty-domain:" + key + ":" + ex[key]
    if reason in ("ratio_below_37_over_100", "positive_branch_Holder_below_37_over_100"):
        value = ex[reason]
        return (True, value.upper().str(35)) if value < rational(COMPETITIVE) else (False, None)
    if reason.endswith("_strict_positive"):
        key = reason[:-len("_strict_positive")]
        value = ex[key]
        return (True, value.lower().str(35)) if value > 0 else (False, None)
    if reason.endswith("_strict_negative"):
        key = reason[:-len("_strict_negative")]
        value = ex[key]
        return (True, value.upper().str(35)) if value < 0 else (False, None)
    raise ValueError("unknown pinned producer leaf reason")


def recheck(box, owner, reason):
                                                                     
    if reason not in REASONS or (reason in REASONS[10:] and owner not in (0, 2)):
        raise ValueError("unknown or inapplicable pinned producer leaf reason")
    return check_reason(independent_bounds(box, owner), owner, reason)


def classify(box, owner, original_reason):
    ex = independent_bounds(box, owner)
                                                                        
                                                                    
    for empty in REASONS[:2]:
        if empty in ex:
            return True, ex[empty], empty
    valid, witness = check_reason(ex, owner, original_reason)
    if valid:
        return True, witness, original_reason
                                                                           
                                                                            
                                                                           
    for reason in REASONS[2:]:
        if reason == original_reason or (reason in REASONS[10:] and owner == 1):
            continue
        valid, witness = check_reason(ex, owner, reason)
        if valid:
            return True, witness, reason
    return False, None, None


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def packed(value):
    raw = canonical(value)
    zipped = zlib.compress(raw, 9)
    encoded = base64.b64encode(zipped).decode("ascii")
    return raw, zipped, encoded


def unpacked(artifact, digest):
    if artifact["sha256_json"] != digest:
        raise ValueError("source certificate does not match immutable digest")
    blob = base64.b64decode(artifact["data_b64"], validate=True)
    if len(artifact["data_b64"]) > BASE64_LIMIT or len(blob) != artifact["zlib_bytes"]:
        raise ValueError("compressed source certificate bounds mismatch")
    decompressor = zlib.decompressobj()
    raw = decompressor.decompress(blob, 32 * 1024 * 1024 + 1)
    if len(raw) > 32 * 1024 * 1024 or not decompressor.eof or decompressor.unused_data:
        raise ValueError("source certificate exceeds bounded JSON or has trailing data")
    if len(raw) != artifact["json_bytes"] or hashlib.sha256(raw).hexdigest() != artifact["sha256_json"]:
        raise ValueError("source certificate digest mismatch")
    return json.loads(raw)


def validate_path(owner, path):
    if type(owner) is not int or owner not in (0, 1, 2) or not isinstance(path, str):
        raise ValueError("invalid owner/path type")
    if len(path) > 160 or any(c not in "012345" for c in path):
        raise ValueError("invalid bounded exact path")


def load_compound(payload):
    base, continuation = payload["base"], payload["continuation"]
    for part, source, trace in ((base, BASE_SOURCE_PIN, BASE_TRACE_PIN),
                               (continuation, CONTINUATION_SOURCE_PIN, CONTINUATION_TRACE_PIN)):
        if part["source_sha256"] != source or part["accepted_trace_sha256"] != trace:
            raise ValueError("source/trace generation does not match immutable replay pin")
    first = unpacked(base["lossless_certificate"], BASE_CERTIFICATE_PIN)
    second = unpacked(continuation["lossless_certificate"], CONTINUATION_CERTIFICATE_PIN)
    if first["schema"] != "n2-binary-conditional-critical-exact-forest.v1":
        raise ValueError("source certificate schema mismatch")
    if second["schema"] != "n2-binary-critical-continuation.compact.v1":
        raise ValueError("continuation certificate schema mismatch")
    if second["base_certificate_sha256"] != BASE_CERTIFICATE_PIN or second["reason_codes"] != REASONS:
        raise ValueError("continuation acceptance dictionary/base reference mismatch")
    rows, pending = first["leaves"], first["frontier"]
    if len(rows) != 32256 or len(pending) != 24 or len(second["leaves"]) != 5031:
        raise ValueError("immutable compound leaf counts mismatch")
    if second["initial_roots"] != pending:
        raise ValueError("continuation initial roots are not exactly the base frontier")
    trace = hashlib.sha256()
    for row in rows:
        validate_path(row["owner"], row["path"])
        if row["reason"] not in REASONS[:10]:
            raise ValueError("unknown base acceptance condition")
        trace.update((json.dumps(row, sort_keys=True) + "\n").encode())
    if trace.hexdigest() != BASE_TRACE_PIN:
        raise ValueError("accepted source trace mismatch")
    for r in pending:
        validate_path(r["owner"], r["path"])
    for owner in range(3):
        forest([r["path"] for r in rows if r["owner"] == owner]
               + [r["path"] for r in pending if r["owner"] == owner])
    extra, extra_pending = [], []
    trace = hashlib.sha256()
    for row in second["leaves"]:
        if not isinstance(row, list) or len(row) != 3:
            raise ValueError("invalid compact continuation leaf")
        owner, path, code = row
        validate_path(owner, path)
        if owner not in (0, 2) or type(code) is not int or not 0 <= code < len(REASONS):
            raise ValueError("invalid continuation owner/reason code")
        trace.update(canonical(row) + b"\n")
        extra.append({"owner": owner, "path": path, "reason": REASONS[code]})
    if trace.hexdigest() != CONTINUATION_TRACE_PIN:
        raise ValueError("accepted continuation trace mismatch")
    for row in second["frontier"]:
        if not isinstance(row, list) or len(row) != 2:
            raise ValueError("invalid compact continuation frontier")
        owner, path = row
        validate_path(owner, path)
        extra_pending.append({"owner": owner, "path": path})
                                                                         
    relative_groups = {(r["owner"], r["path"]): [] for r in pending}
    for row in extra + extra_pending:
        ancestors = [(owner, prefix) for owner, prefix in relative_groups
                     if owner == row["owner"] and row["path"].startswith(prefix)]
        if len(ancestors) != 1:
            raise ValueError("continuation escaped or overlaps initial roots")
        key = ancestors[0]
        relative_groups[key].append(row["path"][len(key[1]):])
    for paths in relative_groups.values():
        forest(paths)
    combined = [{"owner": r["owner"], "path": r["path"], "reason": r["reason"]} for r in rows] + extra
    for owner in range(3):
        forest([r["path"] for r in combined if r["owner"] == owner]
               + [r["path"] for r in extra_pending if r["owner"] == owner])
    return combined, extra_pending


def float_value(point, owner):
    p, q, t = point
    V = p * (1 - p) + t * t * q * (1 - q)
    R = p * (1 - p) * (p * p + (1 - p) ** 2) + t ** 3 * q * (1 - q) * (q * q + (1 - q) ** 2)
    y = [-p - t * q, -p + t * (1 - q), 1 - p - t * q][owner]
    F = [(1 - p) * (1 - q), 1 - p, 1 - p * q][owner]
    D = F - .5 * math.erfc(-y / math.sqrt(2 * V))
    return D * V ** 1.5 / R, V ** 1.5 / R


def fixed_controls():
    pointq = [Q(41, 100), Q(63, 100), Q(57, 100)]
    point = [float(v) for v in pointq]
    answers = []
    for owner in range(3):
        ex = independent_bounds(tuple((v, v) for v in pointq), owner)
        value, pref = float_value(point, owner)
        residuals = []
        for axis in range(3):
            xp, xm = list(point), list(point)
            xp[axis] += 1e-5
            xm[axis] -= 1e-5
            residuals.append((float_value(xp, owner)[0] - float_value(xm, owner)[0]) / (2e-5 * pref))
        ppass = abs(ex["p_gradient_score"] - arb(str(residuals[0]))) < arb("1e-7")
        qpass = abs(ex["q_gradient_score"] - arb(str(residuals[1]))) < arb("1e-7")
        factor = 1 - point[1] if owner == 1 else point[1]
        gpass = abs(ex["unmultiplied_gap"] - arb(str(residuals[2] / factor))) < arb("1e-7")
        cancellation = abs(ex["factored_gap_numerator"]
                           - ex["positive_gap_multiplier"] * ex["unmultiplied_gap"]) < arb("1e-40")
        fullpass = (owner == 1 or abs(ex["factored_gap_numerator"]
                    - (rational(pointq[0]) if owner == 0 else 1 - rational(pointq[0]))
                    * ex["fully_factored_gap"]) < arb("1e-40"))
        vpass = abs(ex["ratio_below_37_over_100"] - arb(str(value))) < arb("1e-12")
        answers.append({"owner": owner, "value": bool(vpass), "p": bool(ppass), "q": bool(qpass),
                        "gap": bool(gpass), "cancelled_gap_identity": bool(cancellation),
                        "additional_factor": bool(fullpass)})
    forest(["0", "1"])
    invalid = False
    try:
        forest(["0"])
    except ValueError:
        invalid = True
    overlap = False
    try:
        forest(["0", "1", "2", "3"])
    except ValueError:
        overlap = True
    unknown = False
    try:
        recheck(((Q(0), Q(0)),) * 3, 0, "invented_condition")
    except ValueError:
        unknown = True
    sample = {"test": ["0", "1"]}
    raw, zipped, encoded = packed(sample)
    serialization = zlib.decompress(base64.b64decode(encoded)) == raw
    return {"points": answers, "incomplete_forest_rejected": invalid,
            "different_axes_overlap_rejected": overlap, "unknown_condition_rejected": unknown,
            "serialization_roundtrip": serialization,
            "all_passed": invalid and overlap and unknown and serialization
                          and all(all(v for k, v in r.items() if k != "owner") for r in answers)}


def run(payload):
    begun = time.monotonic()
    ctx.prec = min(max(int(payload.get("precision_bits", 256)), 192), 512)
    controls = fixed_controls()
    if not controls["all_passed"]:
        return {"status": "independent_controls_failed", "controls": controls, "complete": False}
    source, producer_pending = load_compound(payload)
    seconds = min(float(payload.get("max_seconds", 140)), 140)
    node_cap = min(int(payload.get("max_nodes", 250000)), 300000)
    extra_depth = min(int(payload.get("max_added_depth", 8)), 12)
    descendants, unresolved, samples = [], [], []
    fallback_reasons, acceptance_reasons = {}, {}
    relative_groups = {}
    trace = hashlib.sha256()
    processed = verified_sources = 0
    stop_reason = None
    next_unstarted = len(source)
    started_indices = []
    for index, original in enumerate(source):
        if time.monotonic() >= begun + seconds - 5 or processed >= node_cap or stop_reason == "certificate_storage_guard":
            next_unstarted = index
            break
        started_indices.append(index)
        relative_groups[index] = []
        rootbox = decode(original["path"])
        stack = [("", rootbox)]
        accepted_paths, pending_paths = [], []
        while stack:
            if time.monotonic() >= begun + seconds - 5 or processed >= node_cap or stop_reason == "certificate_storage_guard":
                pending_paths.extend(path for path, _ in stack)
                stop_reason = stop_reason or "global_budget"
                break
            path, box = stack.pop()
            processed += 1
            valid, witness, accepted_reason = classify(box, original["owner"], original["reason"])
            if valid:
                row = [index, path, REASONS.index(accepted_reason)]
                descendants.append(row)
                accepted_paths.append(path)
                relative_groups[index].append(path)
                trace.update(canonical(row) + b"\n")
                acceptance_reasons[accepted_reason] = acceptance_reasons.get(accepted_reason, 0) + 1
                if accepted_reason != original["reason"]:
                    fallback_reasons[accepted_reason] = fallback_reasons.get(accepted_reason, 0) + 1
                if len(samples) < 16:
                    samples.append({"source_leaf": index, "relative_path": path,
                                    "owner": original["owner"], "original_reason": original["reason"],
                                    "accepted_reason": accepted_reason, "witness": witness})
                if len(descendants) % 2048 == 0:
                    _, _, encoded = packed({"descendants": descendants})
                    if len(encoded) > 1048576:
                        stop_reason = "certificate_storage_guard"
                continue
            if len(path) >= extra_depth:
                pending_paths.append(path)
                continue
            axis = max(range(3), key=lambda i: box[i][1] - box[i][0])
            stack.extend((path + str(2 * axis + side), child) for side, child in enumerate(split(box, axis)))
        forest(accepted_paths + pending_paths)
        if not pending_paths:
            verified_sources += 1
        unresolved.extend([index, path] for path in pending_paths)
        relative_groups[index].extend(pending_paths)
        if stop_reason in ("certificate_storage_guard", "global_budget"):
            next_unstarted = index + 1
            break
    replay_certificate = {"schema": "n2-binary-independent-compound-replay.compact.v2",
                          "base_certificate_sha256": BASE_CERTIFICATE_PIN,
                          "continuation_certificate_sha256": CONTINUATION_CERTIFICATE_PIN,
                          "reason_codes": REASONS,
                          "descendants": descendants, "unresolved": unresolved,
                          "unstarted_source_index_range": [next_unstarted, len(source)]}
    raw, zipped, encoded = packed(replay_certificate)
    if len(encoded) > BASE64_LIMIT:
        raise ValueError("independent lossless certificate output cap exceeded")
    roundtrip = json.loads(zlib.decompress(base64.b64decode(encoded)))
    assert roundtrip == replay_certificate
    decoded_groups = {index: [] for index in started_indices}
    for index, path, code in roundtrip["descendants"]:
        if type(code) is not int or not 0 <= code < len(REASONS):
            raise ValueError("decoded acceptance reason invalid")
        decoded_groups[index].append(path)
    for index, path in roundtrip["unresolved"]:
        decoded_groups[index].append(path)
    if decoded_groups != relative_groups:
        raise ValueError("decoded relative partitions differ from replay trace")
    for paths in decoded_groups.values():
        forest(paths)
    complete = not producer_pending and not unresolved and next_unstarted == len(source)
    return {"status": "bounded independent n2 binary critical-condition replay",
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "base_source_sha256": BASE_SOURCE_PIN, "continuation_source_sha256": CONTINUATION_SOURCE_PIN,
            "base_certificate_sha256": BASE_CERTIFICATE_PIN, "continuation_certificate_sha256": CONTINUATION_CERTIFICATE_PIN,
            "controls": controls,
            "precision_bits": ctx.prec, "source_leaf_count": len(source),
            "verified_source_leaf_count": verified_sources, "processed_nodes": processed,
            "descendant_leaf_count": len(descendants), "unresolved_descendant_count": len(unresolved),
            "unstarted_source_index_range": [next_unstarted, len(source)],
            "producer_frontier_count": len(producer_pending),
            "continuation_roots_equal_base_frontier": True, "both_accepted_source_traces_verified": True,
            "original_exact_forests_verified": True, "every_started_source_subforest_verified": True,
            "independent_trace_sha256": trace.hexdigest(), "complete_independent_arithmetic_cover": complete,
            "acceptance_reasons": acceptance_reasons, "separately_proved_fallback_reasons": fallback_reasons,
            "theorem_complete": False, "stop_reason": stop_reason,
            "elapsed_seconds": time.monotonic() - begun,
            "resources_requested": {"cpu": 1, "memory_gib": 1},
            "witness_samples": samples,
            "lossless_certificate": {"encoding": "zlib9+base64+canonical-json", "data_b64": encoded,
                "json_bytes": len(raw), "zlib_bytes": len(zipped), "base64_bytes": len(encoded),
                "sha256_json": hashlib.sha256(raw).hexdigest(), "remote_roundtrip_verified": True},
            "scope": "independent raw normalization/derivative arithmetic and exact forest replay of the pinned producer, conditional on analytic binary contact-domain proof; no unrestricted support reduction"}


if __name__ == "__main__":
    import sys
    print(json.dumps(run(json.loads(sys.stdin.read())), allow_nan=False))
