"""Check the current closed moment cover and its saved interval coefficients."""
from collections import Counter, OrderedDict
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import zipfile

from scalars import scalar_bounds

ROOT=Path(__file__).resolve().parent.parent
KERNEL=ROOT/'verify/kernel'
sys.path.insert(0,str(KERNEL))
UPPER=Fraction(8104460571768603,18014398509481984)
CHOICES=((1.,0.),(1.02,-.02),(1.04,-.02))
CUTOFFS=(.75,.78125,.8125,.84375,.875,.90625,.9375,.96875,1.)

class CertificateError(ValueError):pass

def require(condition,message):
    if not condition:raise CertificateError(message)

def sha256(data):return hashlib.sha256(data).hexdigest()
def file_sha256(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def finite(value):
    require(type(value) in (int,float) and math.isfinite(value),'Expected a finite real number')
    return Fraction(value)

def load_inputs(root=ROOT):
    require(__debug__,'Run without Python optimization')
    result=json.loads((root/'result.json').read_text());raw=(root/result['cover_file']).read_bytes()
    require(sha256(raw)==result['cover_sha256'],'Cover index digest mismatch');cover=json.loads(raw)
    for name,pin in cover['kernel_sha256'].items():
        require(Path(name).name==name and Path(name).suffix in ('.py','.json'),'Invalid numerical source name')
        require(file_sha256(root/'verify/kernel'/name)==pin,'Numerical source changed: '+name)
    require(file_sha256(root/cover['band_archive'])==cover['band_archive_sha256'],'Changed band archive')
    for row in cover['gaussian_archives']:
        require(file_sha256(root/row['file'])==row['sha256'],'Changed Gaussian archive')
    return cover,result

def load_band(bundle,row):
    require(row['file'].startswith('bands/') and Path(row['file']).suffix=='.json','Invalid interval filename')
    raw=bundle.read(row['file']);require(sha256(raw)==row['sha256'],'Interval record changed')
    return json.loads(raw)

def validate_evaluators(evaluators):
    require(evaluators,'Empty evaluator collection')
    for ident,value in evaluators.items():
        require(ident==sha256(canonical(value)),'Evaluator descriptor changed')
        require(value['kind'] in ('quadratic','quadratic-ceil','signed','signed-tail') and value['N'] in (512,1024,2048),'Unknown evaluator')
        if value['kind'].startswith('signed'):
            require(value==dict(kind=value['kind'],N=512,stop_target=value['stop_target'],b_max_depth=6,b_max_evaluations=127,
                               coupling=False,cap_bucket=17,theta_denominator=32,threshold_domain=[-2,2],reference_cache_entries=16),
                    'Signed evaluator settings changed')
            require(value['stop_target'] in ('0.450','0.44995','0.4499'),'Unspecified signed early-return level')
        elif value['kind']=='quadratic':require(value['stop_target'] is None,'Quadratic formula has no stopping parameter')
        else:require(value['stop_target']=='0.448','Ceiling portfolio early-return level changed')

def validate_band(band,evaluators,callback=None):
    """Replay each closed split and outward clipping operation."""
    from flint import arb,ctx
    import stable_cover as core
    from interval_bounds import au
    lo,hi=band['Llo'],band['Lhi'];require(0<finite(lo)<finite(hi),'Invalid L interval')
    with ctx.workprec(100):
        root=[0.,min(1.,au(arb(hi)**(arb(2)/3))),0.,hi,0.,min(1.,hi)]
        require(band['root']==root,'The canonical moment root was reduced')
        require(band['method'] in ('partition','uniform','variance'),'Unknown interval method')
        extra=band['extra_clips'];require(extra==sorted(set(extra)),'Repeated clipping location')
        clip_counts=band.get('additional_clip_counts',{p:1 for p in extra})
        require(set(clip_counts)==set(extra) and all(type(n)is int and 1<=n<=3 for n in clip_counts.values()),
                'Invalid repeated outward clipping count')
        stack=[(tuple(root),'')];cursor=excluded=0;seen_extra=set();uppers=[];replayed=0;recipes=set();proof_count=0
        for code in band['tree']:
            require(stack,'Tree extends beyond its root');raw,path=stack.pop()
            box=core.feasible_clip(raw,hi) if band['method']=='partition' else raw
            if path in extra:
                for _ in range(clip_counts[path]):
                    require(box is not None,'Redundant clipping at an empty node');box=core.feasible_clip(box,hi)
                seen_extra.add(path)
            if code=='X':require(box is None,'Excluded feasible moment box');excluded+=1;continue
            require(box is not None,'Accepted empty moment box')
            if code=='A':
                require(cursor<len(band['leaves']),'Missing accepted leaf');leaf=band['leaves'][cursor]
                require(leaf['path']==path and leaf['box']==list(box),'Accepted leaf differs from closed partition')
                from moment_proofs import exact,validate as validate_proof
                upper=exact(leaf['upper']);require(0<=upper<=UPPER,'Leaf exceeds current global upper bound')
                ident=leaf['evaluator'];require(ident in evaluators,'Missing evaluator');recipes.add(ident)
                require(leaf['bound_kind'] in ('leaf','containing_band'),'Unknown upper-bound meaning')
                require(type(leaf['arithmetic_replayed']) is bool and leaf['arithmetic_replayed']==(leaf['bound_kind']=='leaf'),'Reevaluation scope differs')
                if 'proof' in leaf:
                    proof=leaf['proof'];require(proof['box']==leaf['box'] and exact(proof['upper'])==upper,
                                               'Nested moment proof differs from its original domain')
                    proof_callback=None if callback is None else getattr(callback,'proof',None)
                    require(callback is None or callable(proof_callback),'Missing fixed-recipe callback')
                    validate_proof(proof,lo,hi,proof_callback);proof_count+=1
                elif callback is not None:
                    value=callback(evaluators[ident],lo,hi,tuple(box),band['method'])
                    if leaf['bound_kind']=='leaf':require(value==leaf['upper'],'Numerical leaf value differs')
                    else:require(Fraction(value)<=upper,'Numerical leaf exceeds containing-band bound')
                replayed+=leaf['arithmetic_replayed'];uppers.append(upper);cursor+=1
            else:
                require(band['method']=='partition' and code in 'DBT' and len(path)<60,'Invalid moment split')
                left,right=core.children(box,code);stack.extend(((right,path+'1'),(left,path+'0')))
        require(not stack and cursor==len(band['leaves']) and seen_extra==set(extra),'Partition is incomplete')
        if band['method']!='partition':require(band['tree']=='A' and not extra,'Uniform bound has a subdivision')
        require(cursor>0,'Empty complete interval')
        require(band['weights']==sorted(set(band['weights'])) and band['weights'],'Missing interval coefficients')
    return dict(nodes=len(band['tree']),accepted=cursor,infeasible=excluded,replayed_leaves=replayed,
                complete_owner_proofs=proof_count,upper=max(uppers),evaluators=recipes)

def validate_cover(cover,result,root=ROOT):
    require(cover['schema']=='independent-be-cover-v5' and result['schema']=='independent-be-result-v5','Unsupported certificate schema')
    require(result['status']=='interval-bound' and result['sharp_conjecture_resolved'] is False,'Wrong result class')
    require(cover['comparison_target']==result['comparison_target']==result['upper_bound']=='0.44988794','Wrong comparison target')
    require(Fraction(cover['upper_fraction'])==Fraction(result['largest_recorded_upper_fraction'])==finite(result['largest_recorded_upper'])==UPPER,'Exact upper endpoint differs')
    require(UPPER<Fraction('0.44988794') and cover['geometry_precision_bits']==100,'Invalid upper or geometry precision')
    domain=cover['domain']
    for side in ('lower','upper'):
        require(finite(domain[side])==Fraction(domain[side+'_fraction']) and float.fromhex(domain[side+'_hex'])==domain[side],'Exact endpoint differs')
    require(domain['lower']==.0014 and domain['upper']==1.21,'Wrong finite L domain')
    validate_evaluators(cover['evaluators']);counts=Counter(bands=0,nodes=0,accepted=0,infeasible=0);uppers=[];intervals=[];replayed=0;proofs=0
    with zipfile.ZipFile(root/cover['band_archive']) as bundle:
        require(set(bundle.namelist())=={r['file']for r in cover['bands']} and len(bundle.namelist())==len(cover['bands']),'Wrong interval archive members')
        for number,row in enumerate(cover['bands']):
            band=load_band(bundle,row);require(band['index']==number,'Interval order changed')
            local=validate_band(band,cover['evaluators']);counts['bands']+=1
            for key in ('nodes','accepted','infeasible'):require(local[key]==row[key],'Interval count differs');counts[key]+=local[key]
            require(local['upper']==Fraction(row['upper']) and row['interval']==[str(finite(band['Llo'])),str(finite(band['Lhi']))],'Interval bound or domain differs')
            require(set(band['weights'])<=cover['weight_sha256'].keys(),'Missing interval weight set')
            intervals.append(tuple(map(Fraction,row['interval'])));uppers.append(local['upper']);replayed+=local['replayed_leaves'];proofs+=local['complete_owner_proofs']
    cursor=Fraction(domain['lower'])
    for lo,hi in sorted(intervals):
        require(Fraction(domain['lower'])<=lo<hi<=Fraction(domain['upper']) and lo<=cursor,'Finite closed union has a gap')
        cursor=max(cursor,hi)
    require(cursor==Fraction(domain['upper']) and max(uppers)==UPPER,'Incomplete global interval or wrong maximum')
    require(dict(counts)==result['counts'] and counts['bands']==1670,'Aggregate counts differ')
    require(replayed==result['arithmetic_replayed_leaves'] and counts['accepted']-replayed==result['containing_band_bound_leaves']
            and proofs==result['complete_owner_proofs'],'Bound meanings differ')
    analytic=cover['analytic_complements']
    require(analytic['small_L']['endpoint']==domain['lower_fraction'] and analytic['large_L']['endpoint']==domain['upper_fraction'],'Analytic complements do not meet the cover')
    require(Fraction(analytic['small_L']['upper'])<UPPER and Fraction(analytic['large_L']['ratio_upper'])<UPPER,'Analytic complement exceeds global bound')
    return dict(counts=dict(counts),largest_recorded_upper_fraction=str(UPPER),closed_domain=[domain['lower_fraction'],str(cursor)],
                arithmetic_replayed_leaves=replayed,complete_owner_proofs=proofs)

class CoefficientArchive:
    """A bounded cache of the indexed Gaussian archive shards."""
    def __init__(self,cover,root):self.cover=cover;self.root=root;self.opened=OrderedDict()
    def __enter__(self):return self
    def __exit__(self,*args):
        for archive in self.opened.values():archive.close()
    def read(self,name):
        require(name in self.cover['coefficient_members'],'Unindexed coefficient member')
        filename=self.cover['coefficient_members'][name]
        if filename not in self.opened:
            if len(self.opened)>=3:self.opened.popitem(last=False)[1].close()
            self.opened[filename]=zipfile.ZipFile(self.root/'certificates'/filename)
        self.opened.move_to_end(filename);return self.opened[filename].read(name)
    def names(self,prefix):
        found=set()
        for row in self.cover['gaussian_archives']:
            with zipfile.ZipFile(self.root/row['file']) as archive:
                names=archive.namelist();require(len(names)==len(set(names)),'Duplicate coefficient member')
                for name in names:
                    require(self.cover['coefficient_members'].get(name)==Path(row['file']).name,'Coefficient belongs to another shard')
                    if name.startswith(prefix+'/'):
                        require(name not in found,'Repeated coefficient member');found.add(name);yield name
        require(found=={n for n in self.cover['coefficient_members']if n.startswith(prefix+'/')},'Missing coefficient member')

def validate_gaussian(cover,root=ROOT):
    import signed_gaussian
    from flint import arb,ctx
    from interval_bounds import au,two_kernel
    counts=Counter(records=0,threshold_nodes=0,intervals=0,cutoffs=0);used=set();weight_geometry={}
    def member(bundle,prefix,ident,mapping):
        raw=bundle.read(prefix+'/'+ident+'.json');require(sha256(raw)==mapping[ident],'Gaussian file changed')
        value=json.loads(raw);require(sha256(canonical(value))==ident,'Gaussian mathematical identity differs');return value
    with CoefficientArchive(cover,root) as bundle:
        require(len(cover['coefficient_members'])==len(cover['gaussian_sha256'])+len(cover['weight_sha256']),'Unexpected Gaussian member count')
        for name in bundle.names('gaussian'):
            ident=Path(name).stem
            scalar=member(bundle,'gaussian',ident,cover['gaussian_sha256']);checked=signed_gaussian.verify_partition(scalar)
            require(checked['all_passed'],'Gaussian interpolation failed');counts['records']+=1
            counts['threshold_nodes']+=checked['threshold_nodes'];counts['intervals']+=checked['intervals']
        for name in bundle.names('weights'):
            ident=Path(name).stem
            weight=member(bundle,'weights',ident,cover['weight_sha256']);N=weight['N'];hi=weight['Lhi'];expected={}
            require(N in (512,1024,2048) and 0<finite(hi),'Invalid weight geometry')
            weight_geometry[ident]=(hi,N,set(weight['cap_values']))
            for cap in weight['cap_values']:
                require(0<finite(cap),'Invalid smoothing cap')
                for scale,shift in CHOICES:
                    T=float(round(scale*2*math.pi/(hi+cap),8));s0=float(min(.14+.43*hi+shift,5/T))
                    for k in sorted(set([N]+[min(N,max(1,int(round(N*x))))for x in CUTOFFS])):
                        expected[(T,s0,k)]=float(float(s0)*k/N)
            require(len(weight['cutoffs'])==len(expected),'Incomplete smoothing portfolio')
            seen=set()
            with ctx.workprec(100):
                for row in weight['cutoffs']:
                    key=(row['T'],row['original_s'],row['k']);require(key in expected and key not in seen,'Wrong or repeated cutoff');seen.add(key)
                    require(row['N']==N and row['rounded_split']==expected[key],'Changed original cutoff split')
                    scalar=member(bundle,'gaussian',row['scalar_id'],cover['gaussian_sha256']);used.add(row['scalar_id'])
                    require(scalar['T']==row['T'] and scalar['s']==row['rounded_split'],'Scalar belongs to another cutoff')
                    if row['k']==N:require(row['endpoint_correction']==0 and 0<=finite(row['upper'])<=finite(scalar['upper']),'Full Gaussian bound differs')
                    else:
                        exact=arb(row['original_s'])*row['k']/N;rounded=arb(row['rounded_split'])
                        correction=au(abs(exact-rounded)*two_kernel(exact.union(rounded)))
                        require(row['exact_split']=='s*k/N' and finite(row['endpoint_correction'])>=Fraction(correction),'Understated endpoint correction')
                        require(finite(row['upper'])>=finite(scalar['upper'])+finite(row['endpoint_correction']),'Understated Gaussian cutoff')
                    counts['cutoffs']+=1
        require(used==set(cover['gaussian_sha256']),'Unreferenced scalar partition')
    with zipfile.ZipFile(root/cover['band_archive']) as bands:
        for row in cover['bands']:
            band=load_band(bands,row);available={}
            for ident in band['weights']:
                hi,N,caps=weight_geometry[ident];require(hi==band['Lhi'],'Weights belong to another moment interval')
                available.setdefault(N,set()).update(caps)
            for leaf in band['leaves']:
                recipe=cover['evaluators'][leaf['evaluator']];hi=band['Lhi'];N=recipe['N']
                bucket=17 if recipe['kind'].startswith('signed') else (20 if band['method']!='partition' else min(20,max(1,round(20*leaf['box'][-1]/hi))))
                require(N in available and hi*bucket/20 in available[N],'An active leaf lacks its smoothing portfolio')
                def proof_weights(proof):
                    rule=proof['rule']
                    if rule=='evaluation':
                        r=proof['recipe'];require(r['N'] in available and hi*r['cap_bucket']/20 in available[r['N']],
                                                 'A complete owner recipe lacks its smoothing portfolio')
                    elif rule=='original_evaluation':
                        r=proof['evaluator'];b=17 if r['kind'].startswith('signed') else min(20,max(1,round(20*proof['box'][-1]/hi)))
                        require(r['N'] in available and hi*b/20 in available[r['N']],'An inherited proof lacks its smoothing portfolio')
                    elif rule=='feasible_b':proof_weights(proof['proof'])
                    elif rule=='cover':
                        for child in proof['regions']:proof_weights(child['proof'])
                if 'proof' in leaf:proof_weights(leaf['proof'])
    return dict(counts)

def check(root=ROOT):
    cover,result=load_inputs(root);summary=validate_cover(cover,result,root)
    summary['stored_gaussian_partitions']=validate_gaussian(cover,root)
    summary['scalar_enclosures']=scalar_bounds(cover['domain']['lower'],cover['domain']['upper'])
    summary.update(status='current closed cover and stored coefficients verified',numerical_leaf_integrals_recomputed=False,gaussian_node_integrals_recomputed=False)
    return summary

if __name__=='__main__':
    try:print(json.dumps(check(),indent=2))
    except (CertificateError,AssertionError,ArithmeticError,KeyError,ValueError) as error:
        print('Certificate check failed: '+str(error),file=sys.stderr);raise SystemExit(1)
