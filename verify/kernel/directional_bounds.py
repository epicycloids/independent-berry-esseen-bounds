"""Outward evaluation of DIRECTIONAL_STABILITY.md."""
import time
import numpy as np
from flint import arb
from interval_bounds import up,au
from stable_bounds import StableWeights,finite_law_checks


def directional_cells(thi,Bhi,dhi,taulo,tauhi):
    tau=min(Bhi,max(0,taulo))
    x=up(thi*au(arb(dhi).sqrt()))
    q=up(1+up(up(x*x)/3))
    alpha=up(up(2*x)/3)
    excess=max(0,float(up(Bhi-tau)))
    other=au(arb(Bhi)+5*arb(tau)/3)
    first=up(up(alpha*tau)+up(q*excess))
    squared=up(up(first*first)+up(up(up(q*q)*excess)*other))
    return np.minimum(Bhi,up(np.sqrt(np.maximum(0,squared))))


class DirectionalWeights(StableWeights):
    @staticmethod
    def prefactor(thi,Bhi,dhi,taulo,tauhi):
        return directional_cells(thi,Bhi,dhi,taulo,tauhi)


def directional_checks():
    rng=np.random.default_rng(260927);count=0
    for _ in range(100):
        B=float(rng.uniform(.001,2));tl=float(rng.uniform(0,B));th=float(rng.uniform(tl,B))
        d=float(rng.uniform(0,1));ts=np.array([0.,.001,.1,.5,1.,2.,4.,10.])
        out=directional_cells(ts,B,d,tl,th)
        for tau in [tl,(tl+th)/2,th]:
            bb,tt,dd=arb(B),arb(tau),arb(d)
            for t,y in zip(ts,out):
                x=arb(float(t))*dd.sqrt();q=1+x*x/3
                ref=((2*x*tt/3+q*(bb-tt))**2+q*q*(bb-tt)*(bb+5*tt/3)).sqrt()
                assert y==B or arb(float(y))>=ref,(B,d,tl,th,tau,t,y,ref)
                count+=1
    return {'directional_enclosures':count,'all_passed':True}


def pilot():
    start=time.monotonic();checks=directional_checks()
    laws=finite_law_checks(DirectionalWeights,True);rows=[]
    points=[(.3,.079524,.026925768,.282),(.4,.195600146724,.092507500581,.376),
            (.5,.275700450734,.159762546520,.47),(.55,.299414822466,.191336229373,.495),
            (.6,.344675095679,.232355427437,.54),(.7,.445185420716,.349537581321,.63)]
    for L,d,b,tau in points:
        for N in [512,2048]:
            w=DirectionalWeights(L,N);e=2e-12
            rows.append({'L':L,'d':d,'b':b,'tau':tau,'N':N,'T':w.T,'s':w.s,
                'parameter_halfwidth':e,'upper':w.box(L-e,L+e,d-e,d+e,b-e,b+e,tau-e,tau+e)})
    return {'status':'point enclosures and finite-law checks only','checks':checks,
            'finite_laws':laws,'rows':rows,'elapsed_seconds':time.monotonic()-start}
