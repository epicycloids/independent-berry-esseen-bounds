"""Check the stated scope of the saved leaf reevaluations."""
from collections import Counter
import json
import zipfile
from check_saved import ROOT,load_inputs,load_band,require


def validate_replay_scope(cover,result,root=ROOT):
    counts=Counter(individual_leaf_values=0,containing_band_bounds=0)
    with zipfile.ZipFile(root/cover['band_archive']) as archive:
        for row in cover['bands']:
            band=load_band(archive,row)
            for leaf in band['leaves']:
                if leaf['bound_kind']=='leaf':
                    require(leaf['arithmetic_replayed'] is True,'Individual leaf lacks recorded reevaluation')
                    counts['individual_leaf_values']+=1
                else:
                    require(leaf['bound_kind']=='containing_band' and leaf['arithmetic_replayed'] is False,'Containing-band cap relabelled as a reevaluated leaf')
                    counts['containing_band_bounds']+=1
    require(counts['individual_leaf_values']==result['arithmetic_replayed_leaves']==141533,'Reevaluated leaf count differs')
    require(counts['containing_band_bounds']==result['containing_band_bound_leaves'],'Containing-band count differs')
    return dict(counts,reevaluation='separate recorded computations through the frozen evaluator',fresh_computation_in_this_check=False)

def check(root=ROOT):
    cover,result=load_inputs(root);return validate_replay_scope(cover,result,root)

if __name__=='__main__':print(json.dumps(check(),indent=2))
