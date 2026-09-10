"""Check the saved numerical reevaluation of three intervals against the supplied cover."""
import json
from check_saved import ROOT, canonical, load_inputs, require, sha256


def validate_selected(replay, cover, result):
    require(replay["schema"] == "independent-be-selected-replay-v2", "Unsupported selected replay schema")
    require(replay["cover_sha256"] == result["cover_sha256"], "Selected replay belongs to another cover")
    require(replay["whole_union_numerically_replayed"] is False and
            replay["gaussian_node_replay_claimed"] is False, "Unsupported full replay promotion")
    rows = replay["bands"]
    require(len(rows) == 3 and [r["band_index"] for r in rows] == result["selected_replay_band_indices"], "Wrong selected bands")
    total_nodes = total_leaves = 0
    require(len({row["band_index"] for row in rows}) == 3, "Duplicate selected band")
    for row in rows:
        band = cover["bands"][row["band_index"]]
        require(band["type"] == "mixed" and row["band_sha256"] == sha256(canonical(band)), "Selected mathematical record differs")
        record = band["record"]
        require(row["recorded_numerical_replay_passed"] is True, "Failed selected replay")
        for role in ("origin", "continuation"):
            component = record[role]
            require(row[role + "_nodes"] == component["nodes"], "Selected node counts differ")
            require(row[role + "_leaf_chain_sha256"] == component["leaf_chain_sha256"], "Selected numerical chain differs")
            require(row["numerical_calls"][role] == component["accepted"], "Selected numerical leaf counts differ")
        total_nodes += row["origin_nodes"] + row["continuation_nodes"]
        total_leaves += sum(row["numerical_calls"].values())
    require(total_nodes == replay["total_nodes"] == 2689 and
            total_leaves == replay["total_accepted_leaves"] == 1346, "Selected aggregate counts differ")
    return {"recorded_selected_bands":3, "recorded_selected_nodes":total_nodes,
        "recorded_selected_accepted_leaves":total_leaves,
        "fresh_numerical_evaluation_performed_by_this_check":False,
        "whole_union_numerically_replayed":False, "gaussian_node_replay_claimed":False}


def check(root=ROOT):
    cover, result = load_inputs(root)
    replay = json.loads((root / "certificates/selected-replay.json").read_text())
    return validate_selected(replay, cover, result)


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))
