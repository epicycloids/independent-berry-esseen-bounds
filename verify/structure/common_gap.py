






from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.metadata
from pathlib import Path
import sys
import time
import types

from flint import arb, ctx


ORIGINALS = {
    "singleton": "3a44b54786d29e7c591f33936f060d803f516a2a0fd537f34166cb9b87c8a41e",
    "two_summand": "03f90c832b818b16c376571690dbcedbff312abf80bbf47143a2236c59897644",
}
PRICE_LINE = "    ce = (arb(10).sqrt() + 3) / (6 * (2 * arb.pi()).sqrt())"
ALLOWED_PRICES = {Fraction(199,500),Fraction(2,5)}


def load_original(kind, source, price):
    price=Fraction(price)
    if price not in ALLOWED_PRICES:
        raise ValueError("Unreviewed rational price")
    digest=hashlib.sha256(source.encode()).hexdigest()
    if digest!=ORIGINALS[kind]:
        raise ValueError(f"Original {kind} source hash mismatch")
    count=source.count(PRICE_LINE)
    if count!={"singleton":2,"two_summand":1}[kind]:
        raise ValueError(f"Unexpected {kind} price-line count")
                                                                        
                                                                           
    replacement=f"    ce = arb({price.numerator}) / {price.denominator}"
    patched=source.replace(PRICE_LINE,replacement,1)
    name="bernoulli_price_"+kind
    module=types.ModuleType(name)
    module.__file__="<frozen_original_"+kind+">"
    sys.modules[name]=module
    exec(compile(patched,module.__file__,"exec"),module.__dict__)
    return module,dict(original_sha256=digest,
                       modified_formula_sha256=hashlib.sha256(patched.encode()).hexdigest(),
                       modification=f"First gap-evaluator price assignment only: C_E -> {replacement.strip()[5:]}")


def fixed_controls(singleton,two,price):
    price=Fraction(price)
    C=arb(price.numerator)/price.denominator
    floor=4/(arb(1).exp()*arb.pi().sqrt())
    small=8*C/arb(3).sqrt()-floor-1
    left=8*C*arb(3).sqrt()-10/arb.pi()-floor
    right=8*C*(arb(7)/2).sqrt()-arb(35)/(3*arb.pi())-floor
    assert small>arb(1)/247 and left>1 and right>arb(44)/45
    sb=singleton.Box(singleton.CASE_HIGH_P,(
        (Fraction(4633,10000),)*2,(Fraction(4186,10000),)*2,
        (Fraction(0),)*2,(Fraction(1),)*2))
    sg=singleton.gap_enclosure(sb)
    sf=singleton.gap_enclosure(sb,penalty_multiplier=arb(5)/4)
    assert sg is not None and sg>0 and sf is not None and sf<0
    tb=two.Box(1,(Fraction(131,250),)*2,(Fraction(131,250),)*2)
    tg=two.gap_lower(tb)
    tf=two.gap_lower(tb,constant_shift=-arb(3)/50)
    assert tg is not None and tg>0 and tf is not None and tf<0
    return dict(passed=True,price=str(price),small_variance_margin=str(small),
                analytic_tail_left=str(left),analytic_tail_right=str(right),
                singleton_anchor=str(sg),singleton_inflated_penalty=str(sf),
                two_summand_anchor=str(tg),two_summand_lowered_price=str(tf))


def box_data(box):
    if hasattr(box,"case"):
        return dict(case=box.case,depth=box.depth,
                    bounds=[[str(lo),str(hi)] for lo,hi in box.bounds])
    return dict(threshold=box.threshold,depth=box.depth,
                x_bounds=[str(v) for v in box.x_bounds],
                p_bounds=[str(v) for v in box.p_bounds])


def physical_point_check(module,box,price):







    if not hasattr(box,"case"):
        return None
    midbounds=tuple(((lo+hi)/2,)*2 for lo,hi in box.bounds)
    point=module.Box(box.case,midbounds,box.depth)
    p,x,W,f,_,_=module.physical_parameters(point)
    q=1-p;v=p*q;V=x*(1-x)+W+v
    A0=W-(1-f)*q*q;B0=p*p+(1-f)*(W+2*v)
    denom=(2-f)*B0+A0
    if not(V>0 and 1-q*f>0 and denom>0):
        return dict(evaluated=False,reason="singular midpoint")
    y=f*(W+v)/(1-q*f);z=(y-x)/V.sqrt()
    phi=(-z*z/2).exp()/(2*arb.pi()).sqrt()
    price=Fraction(price)
    C=arb(price.numerator)/price.denominator
    gap=2*C*(4*W+2+(1-2*x)*(1-2*p))/V.sqrt()+z*phi*(2+z*z)-5*V*f*B0/denom
    return dict(evaluated=True,actual_scalar_gap=str(gap),
                negative_certified=bool(gap<0),
                scope="Point in the scalar relaxation, not an integer-binomial array")


def cover(kind,module,deadline,max_boxes,max_depth,price):
    if kind=="singleton":
        stack=[module.unit_box(c) for c in module.CASES]
        evaluator=module.gap_enclosure
    else:
        stack=[module.Box(k,(Fraction(0),Fraction(1)),(Fraction(0),Fraction(1))) for k in range(3)]
        evaluator=module.gap_lower
    accepted=processed=deepest=0
    weakest=None
    reason=None
    box=None
    gap=None
    while stack:
        if processed>=max_boxes:
            reason="max_boxes";break
        if processed%128==0 and time.monotonic()>=deadline:
            reason="wall_time";break
        box=stack.pop();processed+=1
        gap=evaluator(box)
        if gap is not None and gap>0:
            accepted+=1;deepest=max(deepest,box.depth)
            lower=gap.lower()
            if weakest is None or lower<weakest[1]:
                weakest=(box,lower)
            continue
        if box.depth>=max_depth:
            reason="max_depth";break
        stack.extend(module.split_box(box))
    complete=reason is None and not stack
    out=dict(certified=complete,reason=reason or "complete",processed=processed,
             accepted=accepted,pending=len(stack),deepest_accepted=deepest,
             max_boxes=max_boxes,max_depth=max_depth,
             scope=("Compact singleton curvature: p(1-p)<=W<=3,A0>=0,all three charts"
                    if kind=="singleton" else
                    "All closed two-bit owner deficits k0,1,2,p,x in[0,1],V>0; deterministic corners by deletion"))
    if weakest:
        out.update(weakest_accepted_box=box_data(weakest[0]),
                   weakest_accepted_lower=str(weakest[1]),
                   margin_scope="Subdivision lower-endpoint diagnostic; not the physical global minimum")
    if not complete and box is not None:
        out.update(last_box=box_data(box),last_gap=None if gap is None else str(gap))
        out["actual_scalar_midpoint"]=physical_point_check(module,box,price)
    return out


def run(payload=None,**overrides):






    params=dict(price="199/500",arb_precision=192,max_seconds=240,
                singleton_max_boxes=2_000_000,singleton_max_depth=42,
                two_summand_max_boxes=1_000_000,two_summand_max_depth=40)
    if payload is not None:params.update(payload)
    params.update(overrides)
    bundle=params.pop("source_bundle")
    price=Fraction(params["price"])
    if price not in ALLOWED_PRICES:
        raise ValueError("This worker supports only the reviewed prices199/500 and2/5")
    params["price"]=str(price)
    if not(64<=int(params["arb_precision"])<=512 and 0<float(params["max_seconds"])<=600):
        raise ValueError("Invalid precision or wall limit")
    if not all(0<int(params[k])<=4_000_000 for k in ["singleton_max_boxes","two_summand_max_boxes"]):
        raise ValueError("Invalid box budget")
    ctx.prec=int(params["arb_precision"])
    start=time.monotonic();deadline=start+float(params["max_seconds"])
    sm,sh=load_original("singleton",bundle["singleton"],price)
    tm,th=load_original("two_summand",bundle["two_summand"],price)
    controls=fixed_controls(sm,tm,price)
                                                                          
                                                                          
    two=cover("two_summand",tm,deadline,int(params["two_summand_max_boxes"]),int(params["two_summand_max_depth"]),price)
    single=cover("singleton",sm,deadline,int(params["singleton_max_boxes"]),int(params["singleton_max_depth"]),price)
    return dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        original_formulas=dict(singleton=sh,two_summand=th),parameters=params,
        environment=dict(python=sys.version.split()[0],python_flint=importlib.metadata.version("python-flint")),
        controls=controls,two_summand=two,singleton=single,
        all_covers_certified=bool(two["certified"] and single["certified"]),
        elapsed_seconds=time.monotonic()-start,
        theorem_scope="A passing price199/500 return plus the separately reviewed finite n3/4 reduction gives common-gap IID-plus-zero selection at n3/4; price2/5 has its existing n>=5 theorem. The universal16/5 third-difference reduction requires price>=2/5 and is NOT proved at199/500. No unrestricted independent rigidity or winning-parameter uniqueness.")
