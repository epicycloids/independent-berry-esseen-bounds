"""Closed moment partitions and sufficient fixed arithmetic evaluations."""
from fractions import Fraction as Q
import math


def require(condition,message):
    if not condition:raise ValueError(message)


def contains(outer,inner):
    return len(outer)==len(inner)==6 and all(Q(outer[i])<=Q(inner[i])<=Q(inner[i+1])<=Q(outer[i+1]) for i in (0,2,4))


def exact(value):
    if isinstance(value,Q):return value
    require(type(value) in (int,float,str),'Invalid exact upper bound')
    if isinstance(value,float):require(math.isfinite(value),'Nonfinite bound')
    return Q(value)


def signed_edges(count):
    if count==129:return [Q(i,32) for i in range(-64,65)]
    old=[Q(i,32) for i in range(-69,70)]
    if count==139:return old
    require(count==145,'Unknown threshold partition')
    return sorted(set(old+[old[i]+Q(k,4)*(old[i+1]-old[i]) for i in (54,55) for k in (1,2,3)]))


def disk_contexts(recipe,lo,hi):
    """Validate saved frequency/threshold rectangles and coefficient-split sums
    using the recorded function bounds.
    """
    edges=signed_edges(len(recipe['theta_edges']));seen=set()
    cutoffs=[384,400,416,432,448,464,480,496,512]
    def normalized(v):
        v=exact(v);return v/Q(lo if v>=0 else hi)
    for context in recipe['contexts']:
        h,col,r=(context[k] for k in ('threshold_index','column','r'))
        require(type(h)is type(col)is type(r)is int and 0<=h<len(edges)-1 and 0<=col<3 and 0<=r<9,
                'Invalid fixed smoothing cell')
        require((h,col,r) not in seen,'Repeated two-disk context');seen.add((h,col,r))
        scale,shift=((1.,0.),(1.02,-.02),(1.04,-.02))[col]
        T=round(scale*2*math.pi/(hi+hi*recipe['cap_bucket']/20),8)
        s=min(.14+.43*hi+shift,5/T);nodes=[Q(T)*Q(s)*i/512 for i in range(513)]
        native=list(range(cutoffs[r]-16,cutoffs[r]));theta=list(map(str,edges[h:h+2]))
        require(context['cutoff']==cutoffs[r] and context['cells']==native
                and context['threshold_interval']==theta and len(context['old_cells'])==16,
                'Changed complete original interval')
        def rectangle(row,cell,frequency,threshold):
            require(row['cell']==cell and list(map(Q,row['frequency_interval']))==frequency
                    and list(map(Q,row['threshold_interval']))==threshold
                    and len(row['coefficient_split'])==2,'A two-disk rectangle is incomplete')
            tuple(map(Q,row['coefficient_split']))
            require(exact(row['CF_upper'])>=0 and exact(row['normalized_ODE_radius_upper'])>=0
                    and exact(row['first_term_upper'])>=0 and exact(row['separate_ODE_upper'])>=0,
                    'Invalid separate norm enclosure')
            width=2*(frequency[1]-frequency[0])
            separate=width*min(exact(row['separate_ODE_upper']),normalized(row['separate_CF_numerator_upper']))
            joint=min(separate,width*(exact(row['first_term_upper'])+normalized(row['signed_numerator_upper'])))
            require(exact(row['separate_upper'])==separate and exact(row['joint_upper'])==joint,
                    'Two-disk exact signed normalization differs')
        require(len(context['coarse_cells'])==16,'Missing coarse original cell')
        for cell,row in zip(native,context['coarse_cells']):rectangle(row,cell,nodes[cell:cell+2],edges[h:h+2])
        quarters=context['threshold_pieces']
        if quarters:
            require(len(quarters)==4,'Incomplete threshold subdivision')
            for k,quarter in enumerate(quarters):
                threshold=[edges[h]+Q(j,4)*(edges[h+1]-edges[h]) for j in (k,k+1)]
                require(list(map(Q,quarter['threshold_interval']))==threshold and len(quarter['cells'])==16,
                        'Threshold quarter is incomplete')
                for cell,row in zip(native,quarter['cells']):
                    require(row['cell']==cell and len(row['pieces'])==8,'Incomplete frequency subdivision')
                    for j,piece in enumerate(row['pieces']):
                        frequency=[nodes[cell]+Q(k,8)*(nodes[cell+1]-nodes[cell]) for k in (j,j+1)]
                        rectangle(piece,cell,frequency,threshold)


def validate(proof,lo,hi,callback=None):
    """Check every closed child, containing alternative and empty-domain claim.

    Numerical evaluation nodes retain the exact input used by the interval
    kernel.  An optional callback recomputes each such fixed input.
    """
    box=proof['box'];upper=exact(proof['upper']);rule=proof['rule']
    require(len(box)==6 and all(type(v) in (float,int) and math.isfinite(v) for v in box)
            and all(Q(box[i])<=Q(box[i+1]) for i in (0,2,4)) and upper>=0,'Invalid closed moment box')
    if rule=='empty':
        d0,d1,b0,b1,t0,t1=map(Q,box)
        if proof['condition']=='moment_constraints':
            from flint import ctx
            import stable_cover
            with ctx.workprec(100):require(stable_cover.feasible_clip(tuple(box),hi) is None,'Excluded feasible child')
        elif proof['condition']=='b_squared':require(d0**3>b1**2,'The third moment does not exclude this child')
        elif proof['condition']=='remainder_resource':
            excess=t0-Q(hi)+b0
            require(excess>0 and excess**2>d1**3,'The remaining-resource inequality does not exclude this child')
        else:raise ValueError('Unknown empty-domain proof')
        require(upper==0,'Empty domain contributes a positive bound')
        return 0
    if rule=='cover':
        axis=proof['axis'];rows=proof['regions'];require(axis in (0,2) and rows,'Unknown or empty partition')
        cursor=Q(box[axis]);bounds=[];count=0
        for row in rows:
            child=row['box'];require(len(child)==6 and all(child[i]==box[i] for i in range(6) if i not in (axis,axis+1))
                and Q(child[axis])==cursor and Q(child[axis])<=Q(child[axis+1]),'Closed partition has a gap or changed coordinate')
            require(contains(row['proof']['box'],child),'Containing proof does not cover its assigned child')
            count+=validate(row['proof'],lo,hi,callback);bounds.append(exact(row['proof']['upper']));cursor=Q(child[axis+1])
        require(cursor==Q(box[axis+1]) and upper==max(bounds),'Incomplete closed partition or wrong maximum')
        return count
    if rule=='feasible_b':
        from flint import ctx
        from full_b_consumer import full_b_box
        with ctx.workprec(100):full=full_b_box(hi,tuple(box))
        require(full is not None and list(full)==proof['full_box']==proof['proof']['box'],
                'The complete feasible third-moment interval changed')
        count=validate(proof['proof'],lo,hi,callback)
        require(upper==exact(proof['proof']['upper']),'Restricted-domain cap differs')
        return count
    if rule=='original_evaluation':
        from check_saved import canonical,sha256,validate_evaluators
        descriptor=proof['evaluator'];validate_evaluators({sha256(canonical(descriptor)):descriptor})
        require(proof['bound_kind'] in ('leaf','containing_band'),'Invalid original evaluation scope')
    elif rule=='evaluation':
        recipe=proof['recipe']
        require(recipe['flavor'] in ('modern','legacy','tail') and recipe['N']==512
                and recipe['method'] in ('native_box','positive','signed','disk_signed')
                and type(recipe['cap_bucket']) is int and 1<=recipe['cap_bucket']<=20,'Unknown fixed consumer')
        require(recipe['stop_target'] is None if recipe['flavor']=='legacy' else recipe['stop_target'] in ('0.450','0.4499','0.44989',0.44999999999999996),
                'Changed fixed consumer target')
        require(contains(recipe['reference_box'],box),'The reference box does not contain the scoring domain')
        if recipe['method'] in ('signed','disk_signed'):
            edges=list(map(Q,recipe['theta_edges']))
            require(edges==signed_edges(len(edges)),'Incomplete signed threshold domain')
            require(recipe.get('consumer','cached')==('cantelli' if recipe['flavor']=='tail' and len(edges)==129 else 'cached'),
                    'Changed original outside-tail consumer')
        if 'adapter' in recipe:
            require(recipe['method']=='signed' and recipe['adapter'] in (
                'cauchy_reference_adapter','cauchy_endpoint_adapter','packed_weight_adapter',
                'packed_weight_general_adapter','packed_frontier_adapter','packed_polynomial_regional_adapter'),
                'Unknown source inequality')
            require(recipe['mode'] in ('baseline','contracted','polynomial'),'Unknown source comparison')
            if recipe['adapter'].startswith('cauchy_'):require(0<Q(recipe['eta'])<=1,'Invalid radius factor')
        if recipe['method']=='disk_signed':
            require(recipe['contexts'],'Empty disk-intersection portfolio');disk_contexts(recipe,lo,hi)
    else:raise ValueError('Unknown moment proof rule')
    if callback is not None:
        actual=exact(callback(proof,lo,hi))
        if rule=='original_evaluation' and proof['bound_kind']=='containing_band':require(actual<=upper,'Containing numerical cap failed')
        else:require(actual==upper,'Fixed numerical evaluation differs')
    return 1
