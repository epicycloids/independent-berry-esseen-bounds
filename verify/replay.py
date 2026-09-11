"""Recompute a specified interval, or a bounded range of its accepted leaves."""
import argparse
from fractions import Fraction
import json
import os
import zipfile
for variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[variable]='1'
from check_saved import ROOT,canonical,load_inputs,load_band,require,validate_band

class Evaluator:
    def __init__(self):self.factories={}
    def __call__(self,descriptor,lo,hi,box,method):
        from flint import arb,ctx
        from interval_bounds import al,au
        import stable_cover
        key=(lo,hi,canonical(descriptor));kind=descriptor['kind']
        with ctx.workprec(100):
            if method=='variance':return au(arb('0.54093655')/arb(lo))
            if key not in self.factories:
                if kind=='quadratic':import quadratic_cover_engine as engine
                elif kind=='quadratic-ceil':import cap_portfolio_engine as engine
                else:import adaptive_b_cover_engine as engine
                engine.install();choose=stable_cover.weight_factory(hi,descriptor['N'],engine.KIND,descriptor['stop_target'])
                weight=choose(hi)
                if kind.startswith('signed'):
                    require(weight.rectangle_class.__name__=='CachedMaximumVarianceRectangles' and not weight.coupling and
                            weight.target==al(arb(descriptor['stop_target'])),'Wrong signed evaluator')
                    if kind=='signed-tail':
                        from cantelli_gaussian_tail import rectangle_class,constant_upper
                        constant_upper();weight.rectangle_class=rectangle_class();weight.rectangles=None
                self.factories[key]=(choose,weight)
            choose,weight=self.factories[key]
            if method=='uniform':return choose(hi).uniform(lo,hi)
            if kind.startswith('quadratic'):return choose(box[-1]).box(lo,hi,*box)
            from full_b_consumer import full_b_box
            from adaptive_b_consumer import partition_upper
            full=full_b_box(hi,box);require(full is not None,'Feasible box lost its full b interval')
            def score(piece):
                if kind=='signed-tail':
                    from cantelli_gaussian_tail import score_box
                    return score_box(weight,lo,hi,piece,full)
                return weight._score_box(lo,hi,piece,full,coupling=False)
            return partition_upper(full,score,weight.target,max_depth=6,max_evaluations=127)[0]

def replay_band(index,first=0,count=None):
    cover,result=load_inputs();require(type(index) is int and 0<=index<len(cover['bands']),'Band index outside the cover')
    with zipfile.ZipFile(ROOT/cover['band_archive']) as archive:band=load_band(archive,cover['bands'][index])
    validate_band(band,cover['evaluators']);leaves=band['leaves'];last=len(leaves) if count is None else min(len(leaves),first+count)
    require(0<=first<last<=len(leaves),'Empty or invalid leaf range');evaluate=Evaluator()
    for leaf in leaves[first:last]:
        value=evaluate(cover['evaluators'][leaf['evaluator']],band['Llo'],band['Lhi'],tuple(leaf['box']),band['method'])
        if leaf['bound_kind']=='leaf':require(value==leaf['upper'],'Fresh numerical leaf value differs')
        else:require(Fraction(value)<=Fraction(leaf['upper']),'Fresh value exceeds containing-band cap')
    return dict(band_index=index,first_leaf=first,next_leaf=last,numerical_leaf_integrals_recomputed=last-first,total_leaves=len(leaves))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--band-index',required=True,type=int)
    parser.add_argument('--first-leaf',default=0,type=int);parser.add_argument('--max-leaves',type=int)
    args=parser.parse_args();print(json.dumps(replay_band(args.band_index,args.first_leaf,args.max_leaves),indent=2))
