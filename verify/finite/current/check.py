"""Verify the current finite certificates without discovering or subdividing cells.

The default checks every frozen decision tree and source hash. --arithmetic
recomputes the accepted inequalities, sequentially, on those same cells.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, io, itertools, json, sys, zipfile
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RECORD=ROOT/'certificates/finite/current-results.json'
ARCHIVE=ROOT/'certificates/finite/current-finite.zip'

def require(condition,message):
    if not condition:raise ValueError(message)
def load(name):
    path=HERE/(name+'.py')
    spec=importlib.util.spec_from_file_location('finite_current_'+name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
    spec.loader.exec_module(module);return module

def binary_box(path,dimensions,widths=None):
    box=[[0,0] for _ in range(dimensions)]
    for position,bit in enumerate(path):
        require(bit in '01','Invalid binary path')
        axis=position%dimensions if widths is None else max(range(dimensions),key=lambda j:widths[j]/(1<<box[j][1]))
        box[axis][0]=2*box[axis][0]+int(bit);box[axis][1]+=1
    return box

def binary_cover(paths):
    require(bool(paths),'Empty original domain')
    previous=None;depth=max(map(len,paths));end=0
    for path in sorted(paths):
        require(previous is None or not path.startswith(previous),'Repeated or overlapping leaf')
        lo=int(path or '0',2)<<(depth-len(path));hi=(int(path or '0',2)+1)<<(depth-len(path))
        require(lo==end,'Gap in original domain');end=hi;previous=path
    require(end==1<<depth,'Incomplete original domain')

def axis_box(case,path):
    box=[(Q(0),Q(99,100))]*2 if case=='endpoint' else [(Q(0),Q(1))]*3
    for token in path:
        require(token in '012345','Invalid axis-side token')
        axis,side=divmod(int(token),2);require(axis<len(box),'Axis outside root')
        lo,hi=box[axis];mid=(lo+hi)/2;box[axis]=(mid,hi) if side else (lo,mid)
    return [[str(a),str(b)] for a,b in box]

def axis_cover(paths):
    stack=[]
    for path in sorted(paths):
        require(not stack or not path.startswith(stack[-1]),'Repeated or overlapping axis path')
        stack.append(path)
        while len(stack)>1:
            a,b=stack[-2:]
            if not a or len(a)!=len(b) or a[:-1]!=b[:-1]:break
            digit=int(a[-1])
            if digit%2 or int(b[-1])!=digit+1:break
            stack[-2:]=[a[:-1]]
    require(stack==[''],'Incomplete axis-side forest')

def recipe_check(recipe,expected):
    # Fixed-axis masks suffice to check exact bisection/projection topology.
    stack=[0];cursor=nodes=leaves=projections=0
    while cursor<len(recipe):
        require(bool(stack),'Recipe extends past closed root');fixed=stack.pop()
        while recipe[cursor]<10:
            axis=recipe[cursor]//2;require(not fixed&(1<<axis),'Repeated face projection')
            fixed|=1<<axis;projections+=1;cursor+=1
            require(cursor<len(recipe),'Missing terminal operation')
        terminal=recipe[cursor];cursor+=1;nodes+=1
        if terminal==10:leaves+=1
        else:
            axis=terminal-11;require(0<=axis<5 and not fixed&(1<<axis),'Invalid exact bisection')
            stack.extend((fixed,fixed))
    require(not stack,'Unfinished recipe frontier')
    require((nodes,leaves,projections)==(expected['processed'],expected['leaves'],expected['projections']),'Recipe count differs')

def replay_recipe(recipe,name):
    from flint import ctx
    worker=load('bernoulli_three');ctx.prec=320
    stack=[(worker.initial_box(),False)];cursor=0
    while cursor<len(recipe):
        box,known=stack.pop()
        while recipe[cursor]<10:
            axis,upper=divmod(recipe[cursor],2);cursor+=1
            if not known:
                info=worker.numerical_bounds(name,box,gradients=True)
                if info['safe']:known=True
                else:
                    derivative=info['gradient'][axis]
                    require(derivative.upper()<=0 if upper else derivative.lower()>=0,'Projection sign failed')
            value=box.bounds[axis][int(upper)];bounds=list(box.bounds);bounds[axis]=(value,value)
            box=worker.Box(tuple(bounds),box.depth)
        terminal=recipe[cursor];cursor+=1
        if terminal==10:
            require(known or worker.numerical_bounds(name,box,gradients=True)['safe'],'Bernoulli leaf failed')
        else:
            axis=terminal-11;lo,hi=box.bounds[axis];mid=(lo+hi)/2
            left=list(box.bounds);right=list(box.bounds);left[axis]=(lo,mid);right[axis]=(mid,hi)
            stack.extend(((worker.Box(tuple(left),box.depth+1),known),(worker.Box(tuple(right),box.depth+1),known)))

def arithmetic_leaf(group,row,module):
    from flint import arb,ctx
    ctx.prec=384 if group in ('lower110','collision221') else 320
    box=tuple(map(tuple,row['box']));case=row['case']
    if group=='margin':
        box=tuple(tuple(Q(x) for x in pair) for pair in box)
        require(module.lower(case,box)>arb(3)/80,'Bernoulli deficit leaf failed');return
    if group in ('collision210','tt221_negative','tt221_positive'):
        result=module.evaluate(case,box,sweeps=row.get('sweeps',10))
        require(bool(result.get('reason')),f'Raw contact leaf failed: {result}');return
    if group=='low210':
        require(module.independent_rejection(box) is not None,'Algebraic contact leaf failed');return
    args=('221',case,box) if group=='collision221' else (case,box)
    values=module.independently_evaluate(*args)
    for name,kind,value in values:
        # At C=49/120 the lower-floor room may vanish. Its rejection must
        # be strictly negative; this is the accepted equality-floor check.
        if name in ('global_ratio_lower_room','ratio_lower_room'):kind='nonnegative'
        if ((kind=='positive' and value<=0) or (kind=='nonnegative' and value<0)
            or (kind=='nonpositive' and value>0) or (kind=='zero' and (value>0 or value<0))):return
    raise ValueError('No valid closed-floor rejection')

def scalar_checks(archive,arithmetic):
    from flint import arb,arb_series,ctx
    ctx.prec,ctx.cap=320,2
    two=load('two_owner321');three=load('three_owner321')
    def enclosure(fn,raw,sqrt_axis=None):
        ends=[(two.a(lo),two.a(hi)) for lo,hi in raw]
        if sqrt_axis is not None:ends[sqrt_axis]=tuple(t.sqrt() for t in ends[sqrt_axis])
        boxes=[lo.union(hi) for lo,hi in ends];centers=[(lo+hi)/2 for lo,hi in ends];radii=[(hi-lo)/2 for lo,hi in ends]
        value=fn(centers)
        for axis,rad in enumerate(radii):
            value+=fn([arb_series([x,int(i==axis)]) for i,x in enumerate(boxes)])[1]*rad.union(-rad)
        require(min(value.upper(),fn(boxes).upper())<=0,'Scalar contact room failed')
    for group in ('three_owner321','endpoint_owner321','adjacent_owner321'):
        rows=json.loads(archive.read(group+'.json'))['accepted'];observed=set()
        if group=='adjacent_owner321':
            completion=json.loads(archive.read('adjacent_owner321_completion.json'))['accepted']
            for row in completion:
                key=tuple(row[k] for k in ('i','j','k'));require(key not in observed,'Duplicate Cartesian owner');observed.add(key)
                require({tuple(c['child']) for c in row['children']}==set(itertools.product(range(2),repeat=3)),'Missing fixed subcell')
                if arithmetic:
                    raw=[two.mesh(Q(1,6),Q(1),key[0],16),two.mesh(two.V0,Q(1,2),key[1],16),two.mesh(Q(0),Q(3,5),key[2],16)]
                    for c in row['children']:
                        small=[two.mesh(lo,hi,j,2) for (lo,hi),j in zip(raw,c['child'])]
                        enclosure(lambda xs:two.room([xs[0],xs[1].sqrt(),xs[2]],'adjacent',row['witness']),small)
        for row in rows:
            if group=='adjacent_owner321':
                key=tuple(row[k] for k in ('i','j','k'));raw=[two.mesh(Q(1,6),Q(1),key[0],16),two.mesh(two.V0,Q(1,2),key[1],16),two.mesh(Q(0),Q(3,5),key[2],16)]
                fn=lambda xs:two.room(xs,'adjacent',row['witness']);sqrt_axis=1
            elif group=='endpoint_owner321':
                key=(row['case'],row['i'],row['j']);coefficient,slack,z0,z1,nv,nz=two.CASES[row['case']]
                raw=[two.mesh(two.V0,Q(1,2),row['i'],nv),two.mesh(z0,z1,row['j'],nz)]
                fn=lambda xs:two.room(xs,'endpoint',row['witness'],coefficient,slack);sqrt_axis=0
            else:
                key=(row['case'],row['i'],row['j']);coefficient,slack={'small_T':(40,Q(0)),'large_T':(24,Q(3,100))}[row['case']]
                raw=[two.mesh(Q(2401,24964),Q(11,40),row['i'],32),two.mesh(Q(-5,4),Q(3,5),row['j'],64)]
                fn=lambda xs:three.original(*xs,coefficient,slack,row['witness']);sqrt_axis=0
            require(key not in observed,'Duplicate Cartesian cell');observed.add(key)
            if row['witness']=='cantelli':require(raw[-1][1]<=0,'Cantelli outside negative thresholds')
            if arithmetic:enclosure(fn,raw,sqrt_axis)
        if group=='adjacent_owner321':expected=set(itertools.product(range(16),repeat=3))
        elif group=='endpoint_owner321':expected={(name,i,j) for name,(_,_,_,_,nv,nz) in two.CASES.items() for i in range(nv) for j in range(nz)}
        else:expected={(name,i,j) for name in ('small_T','large_T') for i in range(32) for j in range(64)}
        require(observed==expected,'Incomplete Cartesian cover')
    if arithmetic:
        load('rectangle').verify()
        mixture=load('mixture')
        for i in range(256):
            for kind,cap in [('power',Q(113,2000)),('log',Q(1,70))]:
                require(mixture.upper(Q(i,768),Q(i+1,768),Q(1,3),kind)<two.a(cap),'Mixture envelope failed')

def verify(arithmetic=False,selected=None):
    record=json.loads(RECORD.read_text())
    if selected is not None:
        choices={name for name,item in record['domains'].items() if name.startswith('n3_') or item['file'].endswith('.jsonl')}|{'scalars'}
        require(arithmetic and selected in choices,'Select an arithmetic family from current-results.json, or scalars')
    require(hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()==record['archive_sha256'],'Certificate archive differs')
    for path,source in record['sources'].items():require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==source['sha256'],'Arithmetic source differs: '+path)
    with zipfile.ZipFile(ARCHIVE) as archive:
        for group,domain in record['domains'].items():
            data=archive.read(domain['file']);require(hashlib.sha256(data).hexdigest()==domain['sha256'],'Certificate member differs')
            if group.startswith('n3_'):
                recipe_check(data,domain)
                if arithmetic and (selected is None or group==selected):
                    print('Replaying '+group,flush=True);replay_recipe(data,group[3:])
            elif domain['file'].endswith('.jsonl'):
                module=None
                if arithmetic and (selected is None or group==selected):module=load({'collision210':'collision210_raw_replay_evaluator_v1'}.get(group,group))
                paths={};count=0
                for line in io.BytesIO(data):
                    row=json.loads(line);case=row['case'];path=row['path'];count+=1
                    if group=='margin':expected=axis_box(case,path)
                    else:
                        dimension=4 if group in ('collision221','tt221_negative','tt221_positive') else 3
                        widths=(Q(1,4),Q(1),Q(1,2)) if group=='low210' else ((Q(1),Q(1),Q(2)) if group=='lower110' else None)
                        expected=binary_box(path,dimension,widths)
                    require(row['box']==expected,f'Leaf box differs from exact path: {group}, {path}, {row["box"]}, {expected}')
                    if 'source_leaf_path' in row:require(path.startswith(row['source_leaf_path']),'Leaf outside recorded source owner')
                    paths.setdefault(case,[]).append(path)
                    if module:arithmetic_leaf(group,row,module)
                require(count==domain['leaves'],'Leaf count differs')
                for case,leaves in paths.items():
                    require(len(leaves)==domain['cases'][case],'Case count differs')
                    (axis_cover if group=='margin' else binary_cover)(leaves)
        scalar_checks(archive,arithmetic and selected in (None,'scalars'))
    return record

def controls():
    from math import comb
    coefficients=(4,-48,154,-58,0,-8,0,2)
    values=[]
    for i in range(32):
        power=[sum(Q(coefficients[j]*comb(j,k)*i**(j-k),64**j) for j in range(k,8)) for k in range(8)]
        values.extend(sum(power[k]*Q(comb(l,k),comb(7,k)) for k in range(l+1)) for l in range(8))
    require(min(values)==Q(11964025693,3848290697216),'Joint curvature polynomial differs')
    for bad in ([],['0'],['0','0','1'],['00','1']):
        try:binary_cover(bad)
        except ValueError:pass
        else:raise AssertionError('Malformed tree was accepted')
    binary_cover(['0','10','11']);axis_cover(['0','12','13'])
    try:axis_cover(['0','14','13'])
    except ValueError:pass
    else:raise AssertionError('Mismatched siblings were accepted')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--arithmetic',action='store_true');parser.add_argument('--case');args=parser.parse_args()
    controls();record=verify(args.arithmetic,args.case)
    print('Current finite certificates: hashes and complete exact partitions passed.'+(' Requested frozen arithmetic passed.' if args.arithmetic else ' Arithmetic was not rerun.'))
