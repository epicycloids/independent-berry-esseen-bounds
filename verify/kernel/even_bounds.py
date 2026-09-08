"""Outward prefactor from EVEN_SUPPORT.md; the proof is in that document."""
import math
import time
import numpy as np
from flint import arb
from interval_bounds import up,au
from refined_integral import RefinedWeights,integration_checks
from stable_bounds import finite_law_checks


def even_cells(thi,Bhi,dhi,taulo,tauhi):
    tau=min(Bhi,max(0,taulo))
    x=up(thi*au(arb(dhi).sqrt()))
    q=up(1+up(up(x*x)/3))
    alpha=up(up(2*x)/3)
    excess=max(0,float(up(Bhi-tau)))
    other=au(2*arb(Bhi)+2*arb(tau)/3)
    first=up(alpha*tau)
    squared=up(up(first*first)+up(up(up(q*q)*excess)*other))
    return np.minimum(Bhi,up(np.sqrt(np.maximum(0,squared))))


class EvenWeights(RefinedWeights):
    @staticmethod
    def prefactor(thi,Bhi,dhi,taulo,tauhi):
        return even_cells(thi,Bhi,dhi,taulo,tauhi)


def even_checks():
    rng=np.random.default_rng(260928);count=0;support_count=0
    for _ in range(100):
        B=float(rng.uniform(.001,2));tl=float(rng.uniform(0,B));th=float(rng.uniform(tl,B))
        d=float(rng.uniform(0,1));ts=np.array([0.,.001,.1,.5,1.,2.,4.,10.])
        out=even_cells(ts,B,d,tl,th)
        for tau in [tl,(tl+th)/2,th]:
            bb,tt,dd=arb(B),arb(tau),arb(d)
            for t,y in zip(ts,out):
                x=arb(float(t))*dd.sqrt();q=1+x*x/3
                ref=((2*x*tt/3)**2+q*q*(bb-tt)*(2*bb+2*tt/3)).sqrt()
                assert y==B or arb(float(y))>=ref,(B,d,tl,th,tau,t,y,ref)
                count+=1
    for xx in [.001,.01,.1,.5,1.,1.5,1.9,2.]:
        x=arb(xx);f1=x.cos()-x.sin()/x
        derivative=-(x+1/x)*x.sin()-x.cos()
        for yy in [0.,.001,.1,.5,.9,.999,1.001,1.1,1.5,2.,3.,10.,100.,1000.]:
            y=arb(yy);z=x*y
            difference=z.cos()-y*z.sin()/x-f1-derivative*(y*y-1)/2
            assert difference>=0,(xx,yy,difference)
            support_count+=1
    return {'even_prefactor_enclosures':count,'supporting_quadratic_samples':support_count,
            'all_passed':True}


def pilot():
    began=time.monotonic();checks=even_checks()
    laws=finite_law_checks(EvenWeights,True);rows=[]
    points=[(.3,.122163285004,.042698346158,.27),(.4,.18337994274,.088528580873,.36),
            (.5,.257461449454,.143137636601,.45),(.55,.299414822466,.191336229373,.495),
            (.6,.397750191359,.280850860934,.54),(.7,.454995299361,.359409523445,.595)]
    for L,d,b,tau in points:
        bucket=min(20,max(1,math.ceil(20*tau/L)))
        for N in [128,512,1024]:
            w=EvenWeights(L,N,L*bucket/20);e=2e-12
            rows.append({'L':L,'d':d,'b':b,'tau':tau,'bucket':bucket,'N':N,
                'T':w.T,'s':w.s,'parameter_halfwidth':e,
                'upper':w.box(L-e,L+e,d-e,d+e,b-e,b+e,tau-e,tau+e)})
    return {'status':'point enclosures and finite-law checks only; no new global certificate',
            'checks':checks,'finite_laws':laws,'integration_checks':integration_checks(),
            'rows':rows,'elapsed_seconds':time.monotonic()-began}


def probe(paired_cf=False):
    # This replaces only the floating exploratory prefactor. The reusable
    # optimizer, grid, CF envelopes and original source remain unchanged.
    import log_probe
    def floating_even(t,B,d,tau,directional=False):
        excess=max(0,B-tau);q=1+t*t*d/3
        return np.minimum(B,np.sqrt((2*t*math.sqrt(d)*tau/3)**2
                                     +q*q*excess*(2*B+2*tau/3)))
    log_probe.stability=floating_even
    if paired_cf:
        preceding_single=log_probe.cf_single
        def paired_single(t,v,b):
            tau=v**1.5;excess=max(0,b-tau);c=np.abs(np.cos(t*math.sqrt(v)))
            bound=np.sqrt(c*c+c*excess*t**3/3+t**6*excess*(2*b+2*tau/3)/36)
            return np.minimum(preceding_single(t,v,b),bound)
        log_probe.cf_single=paired_single
    result=log_probe.run_stable(True)
    result['prefactor']='EVEN_SUPPORT.md E3'
    result['paired_cf']=paired_cf
    return result
