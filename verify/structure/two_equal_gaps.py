





from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
import importlib.metadata
from pathlib import Path
import sys
import time
from flint import arb,ctx


PRICE=Q(37,100)


@dataclass(frozen=True)
class Box:
    owner:int
    bounds:tuple[tuple[Q,Q],...]
    depth:int=0


def ball(lo,hi=None):
    if hi is None:hi=lo
    mid=(lo+hi)/2;rad=(hi-lo)/2
    return arb(f"{mid.numerator}/{mid.denominator}",f"{rad.numerator}/{rad.denominator}")


def moments(lo,hi):

    vlo=min(lo*(1-lo),hi*(1-hi))
    vhi=Q(1,4) if lo<=Q(1,2)<=hi else max(lo*(1-lo),hi*(1-hi))
    return (vlo,vhi),(vlo*(1-2*vlo),vhi*(1-2*vhi))


def data(box):
    p=ball(*box.bounds[0])
    x=p if len(box.bounds)==1 else ball(*box.bounds[1])
    pb=box.bounds[0];xb=pb if len(box.bounds)==1 else box.bounds[1]
    pv,pr=moments(*pb);xv,xr=moments(*xb)
    v_bounds=(pv[0]+xv[0],pv[1]+xv[1])
    r_bounds=(pr[0]+xr[0],pr[1]+xr[1])
    V=ball(*v_bounds);R=ball(*r_bounds)
    if 4*v_bounds[1]<PRICE*PRICE:
        return dict(small_variance=True,variance_upper=str(v_bounds[1]))
    if v_bounds[0]<=0:return None
    root=V.sqrt();s=p+x;w=p*x
    z=(box.owner-s)/root
    phi=(-z*z/2).exp()/(2*arb.pi()).sqrt()
    normal=(1+(z/arb(2).sqrt()).erf())/2
    F=[(1-p)*(1-x),1-w,arb(1)][box.owner]
    D=F-normal
    if D<0:
        upper=arb(0)
    else:
        upper=ball(v_bounds[1]).sqrt()*ball(v_bounds[1])*D.upper()/ball(r_bounds[0])
    Fw=[1,-1,0][box.owner];Fs=[-1,0,0][box.owner]
    Vs=1-2*s
    Rw=6-12*s+8*s*s-8*w
    Rs=(1-2*s)**3+(-12+16*s)*w
    Gw=R*(V*Fw+z*phi+3*D)-D*V*Rw
    Gs=2*V*R*Fs+2*root*R*phi+R*z*phi*Vs+3*D*R*Vs-2*D*V*Rs
                                                  
    N=R*(6*Fw*V+z*phi*(3+z*z)+3*D)+8*D*V*V
    return dict(small_variance=False,ratio_upper=upper,Gw=Gw,Gs=Gs,
                stationary_curvature_numerator=N,V=V,R=R,D=D)


def classify(box,diagonal=False):
    d=data(box)
    if d is None:return None,None
    if d["small_variance"]:return "coercive_ratio",d
    if d["ratio_upper"]<ball(PRICE):return "ratio",d
    if diagonal:return None,d
    if d["Gw"]>0 or d["Gw"]<0:return "w_gradient",d
    if d["Gs"]>0 or d["Gs"]<0:return "s_gradient",d
    if d["stationary_curvature_numerator"]>0:return "stationary_curvature",d
    return None,d


def split(box):
    index=max(range(len(box.bounds)),key=lambda i:box.bounds[i][1]-box.bounds[i][0])
    lo,hi=box.bounds[index];mid=(lo+hi)/2
    a=list(box.bounds);b=list(box.bounds)
    a[index]=(lo,mid);b[index]=(mid,hi)
    return Box(box.owner,tuple(a),box.depth+1),Box(box.owner,tuple(b),box.depth+1)


def compact_box(box):
    return dict(owner=box.owner,depth=box.depth,bounds=[[str(v) for v in pair] for pair in box.bounds])


def controls():
                                                                           
    lo=Box(0,((Q(0),)*2,(Q(39,100),)*2))
    high=Box(1,((Q(1),)*2,(Q(39,100),)*2))
    values=[]
    for box in [lo,high]:
        d=data(box)
        actual=d["D"]*d["V"]*d["V"].sqrt()/d["R"]
        assert actual>ball(PRICE)
        assert d["Gw"]<0
        values.append(dict(owner=box.owner,actual_ratio=str(actual),Gw=str(d["Gw"])))
    for k in range(3):
        reason,d=classify(Box(k,((Q(9,20),)*2,)),diagonal=True)
        assert reason=="ratio"
    assert classify(Box(0,((Q(0),Q(1,100)),)*2))[0]=="coercive_ratio"
                                                                               
    p,x=Q(7,20),Q(3,5);s=p+x;w=p*x
    V=p*(1-p)+x*(1-x)
    R=p*(1-p)*(p*p+(1-p)**2)+x*(1-x)*(x*x+(1-x)**2)
    poly=s-3*s*s+4*s**3-2*s**4+(6-12*s+8*s*s)*w-4*w*w
    assert V==s-s*s+2*w and R==poly
    assert 6-12*s+8*s*s-8*w==2*((1-2*p)**2+(1-2*p)*(1-2*x)+(1-2*x)**2)
    return dict(passed=True,price=str(PRICE),boundary_controls=values,
                direct_owner1_upper_bound_invalid=True,
                exact_moment_identities=True,diagonal_anchor=True,coercive_anchor=True)


def cover(diagonal,deadline,max_boxes,max_depth):
    dim=1 if diagonal else 2
    stack=[Box(k,((Q(0),Q(1)),)*dim) for k in range(3)]
    processed=accepted=deepest=0
    counts={k:0 for k in ["coercive_ratio","ratio","w_gradient","s_gradient","stationary_curvature"]}
    reason=None;last=None;lastdata=None
    while stack:
        if processed>=max_boxes:reason="max_boxes";break
        if processed%128==0 and time.monotonic()>=deadline:reason="wall_time";break
        box=stack.pop();processed+=1;last=box
        why,lastdata=classify(box,diagonal)
        if why is not None:
            accepted+=1;counts[why]+=1;deepest=max(deepest,box.depth)
            continue
        if box.depth>=max_depth:reason="max_depth";break
        stack.extend(split(box))
    complete=reason is None and not stack
    out=dict(certified=complete,reason=reason or "complete",processed=processed,
             accepted=accepted,pending=len(stack),deepest_accepted=deepest,criteria_counts=counts,
             max_boxes=max_boxes,max_depth=max_depth,
             scope=("Direct ratio<37/100 on all diagonal probabilities and closed owners0,1,2"
                    if diagonal else
                    "Every heterogeneous interior closed-owner competitive local maximum excluded; square edges only test limits"))
    if not complete and last is not None:
        out["last_box"]=compact_box(last)
        if lastdata is not None:
            out["last_enclosures"]={k:str(v) for k,v in lastdata.items()}
    return out


def run(payload=None,**overrides):
    params=dict(price=str(PRICE),arb_precision=192,max_seconds=145,
                diagonal_max_boxes=100_000,diagonal_max_depth=40,
                full_max_boxes=500_000,full_max_depth=42)
    if payload is not None:params.update(payload)
    params.update(overrides)
    if Q(params["price"])!=PRICE:raise ValueError("Only reviewed price37/100 supported")
    if not(64<=int(params["arb_precision"])<=512 and 0<float(params["max_seconds"])<=600):
        raise ValueError("Invalid resource bound")
    if not all(0<int(params[k])<=2_000_000 for k in ["diagonal_max_boxes","full_max_boxes"]):
        raise ValueError("Invalid box bound")
    if not all(0<int(params[k])<=64 for k in ["diagonal_max_depth","full_max_depth"]):
        raise ValueError("Invalid depth bound")
    ctx.prec=int(params["arb_precision"])
    start=time.monotonic();deadline=start+float(params["max_seconds"])
    checks=controls()
    diagonal=cover(True,deadline,int(params["diagonal_max_boxes"]),int(params["diagonal_max_depth"]))
    full=cover(False,deadline,int(params["full_max_boxes"]),int(params["full_max_depth"]))
    return dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                parameters=params,controls=checks,diagonal=diagonal,full=full,
                all_covers_certified=bool(diagonal["certified"] and full["certified"]),
                elapsed_seconds=time.monotonic()-start,
                environment=dict(python=sys.version.split()[0],python_flint=importlib.metadata.version("python-flint")),
                theorem_scope="Passing covers exclude every two-active-bit common-gap global winner at ratio>=37/100. Classical single-bit optimum exceeds37/100, so all common-gap n2 maximizers have one active bit. No nonbinary or unequal-gap conclusion.")
