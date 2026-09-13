"""Replay fixed moment-domain bounds without choosing new search parameters."""
from fractions import Fraction as Q
import importlib
import math
from moment_proofs import require


class Evaluator:
    def __init__(self):
        self.weights={};self.consumers={}

    def prepared(self,recipe,lo,hi):
        import stable_cover
        flavor=recipe['flavor'];key=(hi,flavor,recipe['stop_target'],recipe['cap_bucket'])
        if key not in self.weights:
            if flavor=='legacy':import quadratic_cover_engine as engine
            else:import adaptive_b_cover_engine as engine
            engine.install()
            choose=stable_cover.weight_factory(hi,512,engine.KIND,recipe['stop_target'])
            # The unrounded coordinate fixes exactly the recorded bucket.
            tau=hi*recipe['cap_bucket']/20
            weight=choose(tau)
            if flavor=='tail':
                from cantelli_gaussian_tail import rectangle_class
                weight.rectangle_class=rectangle_class();weight.rectangles=None
            self.weights[key]=weight
        return self.weights[key]

    def __call__(self,proof,lo,hi):
        from flint import ctx
        if proof['rule']=='original_evaluation':
            from replay import Evaluator as Original
            return Original()(proof['evaluator'],lo,hi,tuple(proof['box']),'partition')
        recipe=proof['recipe'];box=tuple(proof['box']);reference=tuple(recipe['reference_box'])
        with ctx.workprec(100):
            weight=self.prepared(recipe,lo,hi);method=recipe['method'];flavor=recipe['flavor']
            if method=='native_box':
                if flavor=='legacy':return weight.box(lo,hi,*box)
                if flavor=='tail':
                    from cantelli_gaussian_tail import score_box
                    return score_box(weight,lo,hi,box,reference)
                return weight._score_box(lo,hi,box,reference,coupling=False)
            if flavor=='legacy':envelopes=weight.envelopes(lo,hi,box)
            else:
                from full_b_consumer import FullBMaximumMixin
                envelopes=super(FullBMaximumMixin,weight).envelopes(lo,hi,box)
            require(envelopes is not None,'A recorded nonempty score became empty')
            if method=='positive':return float(weight.integrate(*envelopes[:2],lo,envelopes[2]))
            edges=tuple(map(Q,recipe['theta_edges']));key=(id(weight),edges,recipe.get('consumer','cached'))
            if key not in self.consumers:
                from cached_maximum_signed import CachedMaximumVarianceRectangles
                cls=CachedMaximumVarianceRectangles
                if recipe.get('consumer')=='cantelli':
                    from cantelli_gaussian_tail import rectangle_class
                    cls=rectangle_class()
                self.consumers[key]=cls(weight,theta_edges=edges)
            consumer=self.consumers[key]
            for name in ('reference_cache','prepared_reference','score_cache'):
                getattr(consumer,name).clear()
            if 'adapter' in recipe:
                adapter=importlib.import_module(recipe['adapter'])
                parameter=recipe['mode']
                if recipe['adapter'].startswith('cauchy_'):
                    parameter='1' if parameter=='baseline' else recipe['eta']
                answer,_=adapter.score(consumer,lo,hi,box,reference,envelopes,parameter)
                return answer['upper']
            detail=consumer.score(lo,hi,reference,envelopes,detail=True)
            if method=='signed':return detail['upper']
            require(method=='disk_signed','Unknown signed recipe')
            return disk_score(consumer,lo,hi,reference,envelopes,detail,recipe['contexts'])


def signed_arrays(rectangles,lo,hi,reference,envelopes,details):
    """The signed consumer's outward operation order, with cell values exposed."""
    import numpy as np
    import stable_bounds
    from interval_signed import up,normalize,directed_prefix,positive_dot,cell_contributions
    w=rectangles.weights;clipped=stable_bounds.feasible_clip(reference,hi)
    require(clipped is not None,'Empty full reference');dlo,dhi,_,_,taulo,tauhi=clipped
    B,F,lowcf=envelopes
    old_cells=cell_contributions(w,B,lo,lowcf)
    old=float(np.min(w.cutoff_integrals(B,F,lo,lowcf)))
    proposed,eligible=rectangles._proposal(lo,hi,dlo,dhi,taulo,tauhi)
    cells=np.where(eligible,np.minimum(old_cells,proposed),old_cells)
    prefix=directed_prefix(cells);cf=up(w.low*lowcf)
    suffix=np.zeros((rectangles.N+1,rectangles.columns))
    for i in range(rectangles.N-1,-1,-1):suffix[i,:]=up(suffix[i+1,:]+cf[i,:])
    high=np.array([positive_dot(w.high[:,j],F[:,j]) for j in range(rectangles.columns)])
    scores=np.empty_like(rectangles.gamma)
    for r,k in enumerate(w.cutoff_indices):
        common=up(up(suffix[k,:]+high)+rectangles.gamma[:,r,:])
        scores[:,r,:]=up(prefix[:,k,:]+normalize(common,lo,hi))
    require(old==details['old_upper'] and np.min(scores,axis=(1,2)).tolist()==details['per_interval'],
            'The reconstructed signed sum differs')
    return cells,scores


def normalized(value,lo,hi):
    value=Q(value);return value/Q(lo if value>=0 else hi)


def disk_score(consumer,lo,hi,reference,envelopes,detail,contexts):
    """Fixed rational coefficient splits on entire frequency/threshold cells."""
    from flint import arb
    from interval_bounds import au
    weight=consumer.weights;B,_,lowcf=envelopes
    cells,scores=signed_arrays(consumer,lo,hi,reference,envelopes,detail)
    per=list(map(Q,detail['per_interval']))
    def ball(x):
        x=Q(x);return arb(x.numerator)/x.denominator
    def primitive(t):return (t*(t*t/2).exp()-(arb.pi()/2).sqrt()*(t/arb(2).sqrt()).erfi())/2
    for context in contexts:
        h,col,r=context['threshold_index'],context['column'],context['r']
        native=context['cells'];nodes=consumer.nodes[col];T=ball(Q(weight.columns[col].T))
        require(int(weight.cutoff_indices[r])==context['cutoff'] and
                [float(cells[h,i,col]) for i in native]==context['old_cells'] and
                float(scores[h,r,col])==context['old_threshold_upper'],'Changed original cell contributions')
        prefixes={i:sum((arb(float(B[j,col]))*arb(float(weight.inner[j,col])) for j in range(i)),arb(0)) for i in native}
        primitives={i:primitive(ball(nodes[i])) for i in native}
        def rectangle(row):
            i=row['cell'];left,right=map(Q,row['frequency_interval']);a,b=map(Q,row['threshold_interval'])
            require(nodes[i]<=left<right<=nodes[i+1] and consumer.theta[h]<=a<b<=consumer.theta[h+1],
                    'A two-disk rectangle left its original cell')
            t,x=ball(left).union(ball(right)),ball(a).union(ball(b));g=(-t*t/2).exp()
            radius=au(g*(prefixes[i]+arb(float(B[i,col]))*(primitive(t)-primitives[i])))
            u,phase=t/T,x*t;angle=arb.pi()*u
            kr=(1-u)/(2*T);ki=((1-u)*angle.cos()/angle.sin()+1/arb.pi())/(2*T)
            real,imag=kr*phase.cos()+ki*phase.sin(),ki*phase.cos()-kr*phase.sin()
            cf=float(lowcf[i,col]);hr,hj=map(ball,row['coefficient_split'])
            hnorm,cnorm=(hr**2+hj**2).sqrt(),(real**2+imag**2).sqrt()
            first=au(arb(radius)*((real-hr)**2+(imag-hj)**2).sqrt())
            numerator=au(arb(cf)*hnorm-g*hr);ode=au(arb(radius)*cnorm)
            cf_numerator=au(arb(cf)*cnorm-g*real)
            require([cf,radius,first,numerator,ode,cf_numerator]==[row[k] for k in
                ('CF_upper','normalized_ODE_radius_upper','first_term_upper','signed_numerator_upper',
                 'separate_ODE_upper','separate_CF_numerator_upper')],'Interval two-disk arithmetic differs')
            width=2*(right-left);separate=width*min(Q(ode),normalized(cf_numerator,lo,hi))
            joint=min(separate,width*(Q(first)+normalized(numerator,lo,hi)))
            require(separate==Q(row['separate_upper']) and joint==Q(row['joint_upper']),'Two-disk normalization differs')
            return joint
        old=list(map(Q,context['old_cells']))
        require(len(native)==len(old)==len(context['coarse_cells'])==16,'Incomplete native-cell comparison')
        coarse=[min(old[i],rectangle(row)) for i,row in enumerate(context['coarse_cells'])]
        alternatives=[sum(coarse)]
        quarters=context['threshold_pieces']
        if quarters:
            require(len(quarters)==4,'Incomplete threshold refinement');cursor=consumer.theta[h];totals=[]
            for quarter in quarters:
                a,b=map(Q,quarter['threshold_interval']);require(a==cursor and a<b,'Threshold refinement gap');cursor=b
                require(len(quarter['cells'])==16,'Missing refined original cell');total=Q(0)
                for j,cell in enumerate(quarter['cells']):
                    require(cell['cell']==native[j] and len(cell['pieces'])==8,'Incomplete frequency refinement')
                    frequency=nodes[native[j]];subtotal=Q(0)
                    for piece in cell['pieces']:
                        left,right=map(Q,piece['frequency_interval'])
                        require(left==frequency and piece['threshold_interval']==quarter['threshold_interval'],'Frequency refinement gap')
                        frequency=right;subtotal+=rectangle(piece)
                    require(frequency==nodes[native[j]+1],'Frequency refinement is incomplete')
                    total+=min(coarse[j],subtotal)
                totals.append(total)
            require(cursor==consumer.theta[h+1],'Threshold refinement is incomplete');alternatives.append(max(totals))
        updated=Q(context['old_threshold_upper'])-sum(old)+min(alternatives)
        per[h]=min(per[h],updated)
    return min(Q(detail['old_upper']),max(Q(detail['outside_upper']),max(per)))
