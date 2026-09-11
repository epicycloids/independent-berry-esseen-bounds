"""Interval bounds for the uniform deficit of a two-Bernoulli sum."""
from fractions import Fraction as Q

import hashlib

import json

from pathlib import Path

from flint import arb,ctx

CASES=('01','10','endpoint')

def number(x):return arb(x.numerator)/x.denominator


def interval(left,right):return number(left).union(number(right))


def normal(z):return (1+(z/arb(2).sqrt()).erf())/2


def normal_tail(z):return (z/arb(2).sqrt()).erfc()/2


def initial(case):return ((Q(0),Q(99,100)),)*2 if case=='endpoint' else ((Q(0),Q(1)),)*3


def serial(box):return [[str(a),str(b)] for a,b in box]


def decode(case,path):
    out=list(initial(case))
    for token in path:
        axis,side=divmod(int(token),2);assert axis<len(out)
        lo,hi=out[axis];mid=(lo+hi)/2
        out[axis]=(mid,hi) if side else (lo,mid)
    return tuple(out)


def tree_check(paths):
    items=[]
    for path in sorted(paths):
        assert not items or not path.startswith(items[-1])
        items.append(path)
        while len(items)>1:
            left,right=items[-2:]
            if not left or len(left)!=len(right) or left[:-1]!=right[:-1]:break
            char=int(left[-1])
            if char%2 or int(right[-1])!=char+1:break
            items[-2:]=[left[:-1]]
    assert items==['']


def variance_range(bounds):
    lo,hi=bounds;bottom=min(lo*(1-lo),hi*(1-hi))
    top=Q(1,4) if lo<=Q(1,2)<=hi else max(lo*(1-lo),hi*(1-hi))
    return bottom,top


def shape_lower(bounds):
    top=variance_range(bounds)[1]
    if not top:return arb(0)
    a=number(top)
    return ((1-2*a)/a.sqrt()).lower()


def lower(case,box):
    ce=(3+arb(10).sqrt())/(6*(2*arb.pi()).sqrt())
    if case=='endpoint':
        variance_cap=Q(0);odds=Q(0);atom=Q(1)
        for lo,hi in box:
            a=variance_range((lo,hi))[1]
            variance_cap+=a/(1-2*a)**2;odds+=hi/(1-hi);atom*=1-lo
        if not variance_cap:return arb(1)
        return ce/number(variance_cap).sqrt()+normal_tail(number(odds).sqrt())-number(atom)
    p1,p2=(interval(*box[i]) for i in (0,1))
    if case=='01':
        t=interval(*box[2]);av=[variance_range(box[i]) for i in (0,1)]
        a1,a2=(interval(*a) for a in av)
        r1,r2=1-2*a1,1-2*a2
        V=a1*t*t+a2;R=a1*r1*t*t*t+a2*r2
        holder=sum(top/(1-2*top)**2 for _,top in av)
        assert holder>0
        resource=(1/number(holder).sqrt()).lower()
        if V.lower()>0:resource=max(resource,(R/(V*V.sqrt())).lower())
        numerator=1-p2-p1*t
        phi=arb(0)
        if V.lower()>0:phi=normal(numerator/V.sqrt()).lower()
        elif numerator.lower()>=0 and V.upper()>0:
            phi=normal(numerator.lower()/V.upper().sqrt()).lower()
        p1hi=box[0][1]
        if p1hi<1:phi=max(phi,normal_tail(number(p1hi/(1-p1hi)).sqrt()).lower())
        return ce.lower()*resource+phi-number(1-box[0][0]*box[1][0])
    assert case=='10'
    angle_lo,angle_hi=box[2]
    # w is the second variance share. Clip only the known trigonometric range.
    wlo=(arb.pi()*number(angle_lo)/2).sin()**2
    whi=(arb.pi()*number(angle_hi)/2).sin()**2
    w_lower=max(arb(0),wlo.lower());w_upper=min(arb(1),whi.upper())
    first_lower=max(arb(0),1-w_upper).sqrt()
    second_lower=w_lower.sqrt();second_upper=w_upper.sqrt()
    resource=first_lower**3*shape_lower(box[0])+second_lower**3*shape_lower(box[1])
    hi1,hi2=box[0][1],box[1][1]
    if hi2==1:phi=arb(0)
    else:
        positive=arb(0) if hi1==1 else first_lower*number((1-hi1)/hi1).sqrt()
        negative=second_upper*number(hi2/(1-hi2)).sqrt()
        phi=normal(positive-negative).lower()
    return ce.lower()*resource+phi-number(1-box[1][0])


