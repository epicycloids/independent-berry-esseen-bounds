"""Interval bounds for three Bernoulli summands of arbitrary support diameters."""
import base64

from dataclasses import dataclass

from fractions import Fraction as Q

import hashlib

import itertools

import json

from pathlib import Path

from flint import arb, arb_series, ctx

import flint

PROB_CUT=Q(1,125000)

VAR_CUT=Q(1,25)

GAP_CUT=Q(3,500)

RESOURCE_CAP=Q(4,3)

DISCREPANCY_CAP=Q(271,500)

NORMALIZED_GAP_MIN=Q(2,5)

NORMALIZED_GAP_MAX=Q(200,3)

NAMES=("first","second","upper_pair","last","triangle_lower",
       "triangle_upper","separated_lower","separated_upper")

SPECS={
    "first":("O",(1,0,0)),
    "second":("O",(0,1,0)),
    "upper_pair":("O",(1,0,1)),
    "last":("O",(0,1,1)),
    "triangle_lower":("T",(0,0,1)),
    "triangle_upper":("T",(1,1,0)),
    "separated_lower":("U",(1,1,0)),
    "separated_upper":("U",(0,0,1)),
}

@dataclass(frozen=True)
class Box:
    bounds: tuple
    depth: int=0


def number(x):
    if isinstance(x,Q):
        return arb(x.numerator)/x.denominator
    return arb(x)


def interval(lo,hi):
    if isinstance(lo,Q) and isinstance(hi,Q):
        midpoint=(lo+hi)/2
        radius=(hi-lo)/2
        return arb(f"{midpoint.numerator}/{midpoint.denominator}",
                   f"{radius.numerator}/{radius.denominator}")
    return number(lo).union(number(hi))


def cdf(x):
    return (1+(x/arb(2).sqrt()).erf())/2


def density(x):
    return (-x*x/2).exp()/(2*arb.pi()).sqrt()


def gap_map(chart,s,t):
    if chart=="O":
        return (s*t,s,Q(1))
    if chart=="T":
        return ((1-s)/2+s*t,(1+s)/2,Q(1))
    return (s/2,s/2+(1-s)*t,Q(1))


def event_value(name,b):
    b1,b2,b3=b
    if name=="first": return (1-b2)*(1-b3)
    if name=="second": return (1-b3)*(1-b1*b2)
    if name=="upper_pair": return 1-b2*b3
    if name=="last": return 1-b1*b2*b3
    if name=="triangle_lower": return int(b1+b2+b3<=1)
    if name=="separated_lower": return 1-b3
    return (1-b3)+b3*(1-b1)*(1-b2)


def event_probability(name,p):
    p1,p2,p3=p
    q1,q2,q3=1-p1,1-p2,1-p3
    if name=="first": return q2*q3
    if name=="second": return q3*(1-p1*p2)
    if name=="upper_pair": return 1-p2*p3
    if name=="last": return 1-p1*p2*p3
    if name=="triangle_lower": return q1*q2*q3+p1*q2*q3+q1*p2*q3+q1*q2*p3
    if name=="separated_lower": return q3
    return q3+p3*q1*q2


def influences(name,p):
    p1,p2,p3=p
    q1,q2,q3=1-p1,1-p2,1-p3
    if name=="first": return (Q(0),q3,q2)
    if name=="second": return (p2*q3,p1*q3,1-p1*p2)
    if name=="upper_pair": return (Q(0),p3,p2)
    if name=="last": return (p2*p3,p1*p3,p1*p2)
    if name=="triangle_lower": return (p2+p3-2*p2*p3,p1+p3-2*p1*p3,p1+p2-2*p1*p2)
    if name=="separated_lower": return (Q(0),Q(0),Q(1))
    return (p3*q2,p3*q1,1-q1*q2)


def variance_bounds(lo,hi):
    ends=(lo*(1-lo),hi*(1-hi))
    return min(ends),Q(1,4) if lo<=Q(1,2)<=hi else max(ends)


def raw_data(name,box):
    chart,owner=SPECS[name]
    pb=box.bounds[:3]
    sb,tb=box.bounds[3:]
    corners=[gap_map(chart,s,t) for s,t in itertools.product(sb,tb)]
    hb=[(min(h[i] for h in corners),max(h[i] for h in corners)) for i in range(3)]
    ab=[variance_bounds(*b) for b in pb]
    bb=[(lo*(1-2*lo),hi*(1-2*hi)) for lo,hi in ab]
    vl=sum(a[0]*h[0]**2 for a,h in zip(ab,hb))
    vu=sum(a[1]*h[1]**2 for a,h in zip(ab,hb))
    rl=sum(b[0]*h[0]**3 for b,h in zip(bb,hb))
    ru=sum(b[1]*h[1]**3 for b,h in zip(bb,hb))
    yl=yu=Q(0)
    for bit,p,h in zip(owner,pb,hb):
        if bit:
            yl+=(1-p[1])*h[0]
            yu+=(1-p[0])*h[1]
        else:
            yl-=p[1]*h[1]
            yu-=p[0]*h[0]
    return dict(pb=pb,hb=hb,ab=ab,bb=bb,vl=vl,vu=vu,rl=rl,ru=ru,yl=yl,yu=yu)


def numerical_bounds(name,box,*,gradients=False,audit_gradient=False):
    d=raw_data(name,box)
    pb,hb,ab=d["pb"],d["hb"],d["ab"]
    if hb[0]==(Q(1),Q(1)) and not audit_gradient:
        return dict(safe=True,reason="common_diameter",score=1.)
    # Each collar is a proved exclusion of negative values, not a claimed
    # uniform lower enclosure outside the negative set.
    for lo,hi in pb:
        if min(hi,1-lo)<=PROB_CUT:
            return dict(safe=True,reason="probability_collar",score=1.)
    if hb[0][1]<=GAP_CUT:
        return dict(safe=True,reason="gap_collar",score=1.)
    vl,vu=d["vl"],d["vu"]
    if not vl>0:
        raise AssertionError("clipped probabilities must keep variance positive")
    if any(a[1]*h[1]**2<=VAR_CUT*vl for a,h in zip(ab,hb)):
        return dict(safe=True,reason="variance_collar",score=1.)
    vlr,vur=number(vl),number(vu)
    sqrt_vl,sqrt_vu=vlr.sqrt(),vur.sqrt()
    direct_lower=number(d["rl"])/(vur*sqrt_vu)
    reciprocal_upper=sum(a[1]/(1-2*a[1])**2 for a in ab)
    holder=1/number(reciprocal_upper).sqrt()
    resource_lower=max(direct_lower.lower(),holder.lower())
    resource_upper=(number(d["ru"])/(vlr*sqrt_vl)).upper()
    ce=(3+arb(10).sqrt())/(6*(2*arb.pi()).sqrt())
    if ce.lower()*resource_lower>=number(DISCREPANCY_CAP):
        return dict(safe=True,reason="resource",
                    score=float((ce.lower()*resource_lower-number(DISCREPANCY_CAP)).lower()))
    yl,yu=d["yl"],d["yu"]
    xl=(number(yl)/(sqrt_vu if yl>=0 else sqrt_vl)).lower()
    xu=(number(yu)/(sqrt_vl if yu>=0 else sqrt_vu)).upper()
    if xl>=2 or xu<=-2:
        return dict(safe=True,reason="phase",score=.03)
    normal_lower=max(arb(0),cdf(xl).lower())
    Fupper=number(event_probability(name,[b[0] for b in pb])).upper()
    lower=(ce.lower()*resource_lower+normal_lower-Fupper).lower()
    if lower>=0 and not audit_gradient:
        return dict(safe=True,reason="direct",score=float(lower),lower=lower.str(35))
    answer=dict(safe=False,reason="unresolved",score=float(lower),lower=lower.str(35))
    if not gradients:
        return answer
    # Valid only on f<0; the projection lemma permits these intersections.
    xl=max(xl,arb(-2));xu=min(xu,arb(2))
    resource_upper=min(resource_upper,number(RESOURCE_CAP))
    if resource_lower>resource_upper or xl>xu:
        return dict(safe=True,reason="empty_negative_enclosure",score=1.)
    rho=interval(resource_lower,resource_upper)
    phis=[density(xl),density(xu)]
    phi_lower=min(v.lower() for v in phis)
    phi_upper=max(v.upper() for v in phis)
    if xl<=0<=xu:
        phi_upper=density(arb(0)).upper()
    phi=interval(phi_lower,phi_upper)
    xphi_values=[xl*density(xl),xu*density(xu)]
    for crit in (-1,1):
        if xl<=crit<=xu:
            xphi_values.append(arb(crit)*density(arb(crit)))
    xphi=interval(min(v.lower() for v in xphi_values),max(v.upper() for v in xphi_values))
    lam=(3*ce*rho+xphi)/2
    # On a negative point, all normalized variance shares exceed VAR_CUT.
    # Since h_3(raw)=1 and R<4/3, sqrt(V)>3*VAR_CUT/8. Also
    # R>=h_min(raw)/(2sqrt(V)), hence sqrt(V)>3*h_min(raw)/8.
    # These joint restrictions prevent a spurious enormous gap gradient
    # when a broad probability box touches a deterministic corner.
    gradient_scale_lower=max(sqrt_vl.lower(),number(3*VAR_CUT/8),number(3*hb[0][0]/8))
    H=[]
    for h in hb:
        lo=max((number(h[0])/sqrt_vu).lower(),number(NORMALIZED_GAP_MIN))
        hi=min((number(h[1])/gradient_scale_lower).upper(),number(NORMALIZED_GAP_MAX))
        if lo>hi:
            return dict(safe=True,reason="empty_gap_enclosure",score=1.)
        H.append(interval(lo,hi))
    ico=[influences(name,corner) for corner in itertools.product(*pb)]
    I=[interval(min(v[i] for v in ico),max(v[i] for v in ico)) for i in range(3)]
    p=[interval(*b) for b in pb]
    a=[interval(*b) for b in ab]
    b=[interval(*b) for b in d["bb"]]
    t=[interval(1-2*hi,1-2*lo) for lo,hi in pb]
    t3=[interval((1-2*hi)**3,(1-2*lo)**3) for lo,hi in pb]
    owner=SPECS[name][1]
    gp=[ce*t3[i]*H[i]**3-phi*H[i]+I[i]-lam*t[i]*H[i]**2 for i in range(3)]
    invscale=interval((1/sqrt_vu).lower(),(1/gradient_scale_lower).upper())
    gh=[(3*ce*b[i]*H[i]**2+phi*(owner[i]-p[i])-2*lam*a[i]*H[i])*invscale for i in range(3)]
    s=interval(*box.bounds[3]);tt=interval(*box.bounds[4])
    chart=SPECS[name][0]
    if chart=="O":
        gs=tt*gh[0]+gh[1];gt=s*gh[0]
    elif chart=="T":
        gs=(tt-arb(1)/2)*gh[0]+gh[1]/2;gt=s*gh[0]
    else:
        gs=gh[0]/2+(arb(1)/2-tt)*gh[1];gt=(1-s)*gh[1]
    if not audit_gradient and vu<=4*vl:
        centered=centered_lower(name,box)
        if centered>=0:
            return dict(safe=True,reason="full_box_mean_value",score=float(centered),lower=centered.str(35))
    answer["gradient"]=gp+[gs,gt]
    return answer


def full_gradient(name,box):
    """Enclose the five derivatives on the entire box, without clipping."""
    d=raw_data(name,box)
    pb,hb,ab=d["pb"],d["hb"],d["ab"]
    vl,vu=number(d["vl"]),number(d["vu"])
    sl,su=vl.sqrt(),vu.sqrt()
    rho=interval((number(d["rl"])/(vu*su)).lower(),
                 (number(d["ru"])/(vl*sl)).upper())
    yl,yu=d["yl"],d["yu"]
    xl=(number(yl)/(su if yl>=0 else sl)).lower()
    xu=(number(yu)/(sl if yu>=0 else su)).upper()
    values=[density(xl),density(xu)]
    phi_hi=max(v.upper() for v in values)
    if xl<=0<=xu:
        phi_hi=density(arb(0)).upper()
    phi=interval(min(v.lower() for v in values),phi_hi)
    products=[xl*density(xl),xu*density(xu)]
    for critical in (-1,1):
        if xl<=critical<=xu:
            products.append(arb(critical)*density(arb(critical)))
    xphi=interval(min(v.lower() for v in products),max(v.upper() for v in products))
    ce=(3+arb(10).sqrt())/(6*(2*arb.pi()).sqrt())
    lam=(3*ce*rho+xphi)/2
    H=[interval((number(lo)/su).lower(),(number(hi)/sl).upper()) for lo,hi in hb]
    corners=[influences(name,c) for c in itertools.product(*pb)]
    I=[interval(min(v[i] for v in corners),max(v[i] for v in corners)) for i in range(3)]
    p=[interval(*b) for b in pb]
    a=[interval(*b) for b in ab]
    b=[interval(*b) for b in d["bb"]]
    t=[interval(1-2*hi,1-2*lo) for lo,hi in pb]
    t3=[interval((1-2*hi)**3,(1-2*lo)**3) for lo,hi in pb]
    owner=SPECS[name][1]
    gp=[ce*t3[i]*H[i]**3-phi*H[i]+I[i]-lam*t[i]*H[i]**2 for i in range(3)]
    inv=interval((1/su).lower(),(1/sl).upper())
    gh=[(3*ce*b[i]*H[i]**2+phi*(owner[i]-p[i])-2*lam*a[i]*H[i])*inv for i in range(3)]
    s=interval(*box.bounds[3]);tt=interval(*box.bounds[4])
    chart=SPECS[name][0]
    if chart=="O":
        gs=tt*gh[0]+gh[1];gt=s*gh[0]
    elif chart=="T":
        gs=(tt-arb(1)/2)*gh[0]+gh[1]/2;gt=s*gh[0]
    else:
        gs=gh[0]/2+(arb(1)/2-tt)*gh[1];gt=(1-s)*gh[1]
    return gp+[gs,gt]


def centered_lower(name,box):
    midpoint=tuple((lo+hi)/2 for lo,hi in box.bounds)
    value=exact_point_value(name,midpoint)
    for gradient,(lo,hi) in zip(full_gradient(name,box),box.bounds):
        magnitude=max(abs(gradient.lower()),abs(gradient.upper()))
        value-=number((hi-lo)/2)*magnitude
    return value.lower()


def exact_point_value(name,point):
    p=point[:3]
    h=gap_map(SPECS[name][0],*point[3:])
    a=[z*(1-z) for z in p]
    V=sum(z*g*g for z,g in zip(a,h))
    R=sum(z*(1-2*z)*g**3 for z,g in zip(a,h))
    y=sum((bit-z)*g for bit,z,g in zip(SPECS[name][1],p,h))
    ce=(3+arb(10).sqrt())/(6*(2*arb.pi()).sqrt())
    v=number(V)
    return ce*number(R)/(v*v.sqrt())+cdf(number(y)/v.sqrt())-number(event_probability(name,p))


def self_check():
    ctx.prec=192
    controls=0
    if not (Q(25,24)**3<Q(15,13)**2 and Q(25,24)**3<Q(1153,1084)**2
            and Q(64,9)*PROB_CUT<VAR_CUT**3
            and DISCREPANCY_CAP<Q(49,120)*RESOURCE_CAP
            and NORMALIZED_GAP_MIN**2==4*VAR_CUT
            and NORMALIZED_GAP_MAX==2*RESOURCE_CAP/VAR_CUT
            and GAP_CUT==NORMALIZED_GAP_MIN/NORMALIZED_GAP_MAX):
        raise AssertionError("smooth-deletion rational control failed")
    for name in NAMES:
        for s,t in ((Q(1,3),Q(2,3)),(Q(0),Q(1,2)),(Q(1),Q(0)),(Q(1),Q(1))):
            p=(Q(2,7),Q(3,8),Q(4,9))
            h=gap_map(SPECS[name][0],s,t)
            cut=sum(bit*g for bit,g in zip(SPECS[name][1],h))
            physical=Q(0)
            polynomial=Q(0)
            for bits in itertools.product((0,1),repeat=3):
                weight=Q(1)
                for bit,prob in zip(bits,p): weight*=prob if bit else 1-prob
                if sum(bit*g for bit,g in zip(bits,h))<=cut: physical+=weight
                polynomial+=event_value(name,bits)*weight
            if polynomial!=event_probability(name,p) or polynomial>physical:
                raise AssertionError((name,s,t,"CDF ownership mismatch"))
            point=p+(s,t)
            box=Box(tuple((z,z) for z in point))
            out=numerical_bounds(name,box)
            direct=exact_point_value(name,point)
            if out["reason"]=="direct" and not direct>=0:
                raise AssertionError("point enclosure failed")
            controls+=1
    # A sharp-price decrement must produce an actual negative common-gap
    # three-coordinate point, detecting a vacuous always-accept procedure.
    point=(Q(43,100),)*3+(Q(1),Q(1))
    val=exact_point_value("triangle_lower",point)
    p=Q(43,100);a=p*(1-p);rho=number(p*p+(1-p)**2)/number(3*a).sqrt()
    if not (val-arb(1)/20*rho)<0:
        raise AssertionError("lowered-price control failed")
    # Independent formal Taylor differentiation of the raw expression
    # checks all five physical derivatives against the interval formulas.
    # These interior controls satisfy every conditional negative-set bound.
    derivative_controls=0
    for name in NAMES:
        for s,t in ((Q(3,4),Q(1,2)),(Q(4,5),Q(2,3)),(Q(9,10),Q(3,4))):
            point=(Q(2,7),Q(3,8),Q(4,9),s,t)
            out=numerical_bounds(name,Box(tuple((z,z) for z in point)),
                                 gradients=True,audit_gradient=True)
            if "gradient" not in out:
                raise AssertionError((name,point,"invalid derivative control"))
            full=full_gradient(name,Box(tuple((z,z) for z in point)))
            for axis in range(5):
                variables=[arb_series([number(value),int(j==axis)],prec=2)
                           for j,value in enumerate(point)]
                p=variables[:3]
                h=[arb_series([number(g)],prec=2) if isinstance(g,Q) else g
                   for g in gap_map(SPECS[name][0],*variables[3:])]
                a=[v*(1-v) for v in p]
                V=sum(v*g*g for v,g in zip(a,h))
                R=sum(v*(1-2*v)*g**3 for v,g in zip(a,h))
                y=sum((bit-v)*g for bit,v,g in zip(SPECS[name][1],p,h))
                ce=(3+arb(10).sqrt())/(6*(2*arb.pi()).sqrt())
                value=ce*R/(V*V.sqrt())+cdf(y/V.sqrt())-event_probability(name,p)
                if not value[1].overlaps(out["gradient"][axis]):
                    raise AssertionError((name,point,axis,"Taylor derivative mismatch"))
                if not value[1].overlaps(full[axis]):
                    raise AssertionError((name,point,axis,"full-box derivative mismatch"))
                derivative_controls+=1
    return dict(exact_CDF_controls=controls,lowered_price_control=True,
                formal_Taylor_derivative_controls=derivative_controls,
                full_box_derivative_controls=derivative_controls,
                smooth_deletion_rational_controls=True)


def serialize_box(box):
    return dict(bounds=[[str(lo),str(hi)] for lo,hi in box.bounds],depth=box.depth)


def initial_box():
    return Box(((PROB_CUT,1-PROB_CUT),)*3+((Q(0),Q(1)),)*2)


