"""Cached Arb upper integration of the regularized Gaussian kernel."""
import math
import time
from functools import lru_cache
import numpy as np
from flint import arb,acb
from scipy.integrate import quad
from scipy.special import exp1
from interval_bounds import au
from batch_weights import BatchWeights


@lru_cache(maxsize=2048)
def gaussian_G(T_value,s_value,M=8):
    T,s=arb(T_value),arb(s_value);p=arb.pi()
    assert T>0 and 0<s<1
    c=au((p*p/6-sum(arb(1)/(n*n) for n in range(1,M+1)))/(1-s*s/(M+1)**2))
    def integrand(u,analytic):
        h=2*u/p*(sum(1/(n*n-u*u) for n in range(1,M+1))+c)
        return (1-u)*(1+h*h).sqrt(analytic=analytic)*(-T*T*u*u/2).exp()
    value=acb.integral(integrand,acb(0),acb(s),rel_tol=1e-12,abs_tol=1e-14,
                       eval_limit=15000,depth_limit=30,use_heap=True)
    assert value.is_finite() and value.imag.contains(0),value
    total=value.real+(T*T*s*s/2).expint(1)/(2*p)
    assert total>0
    return au(total)


class GaussianBatchWeights(BatchWeights):
    def __init__(self,Lhi,N=512,tau_cap=None):
        super().__init__(Lhi,N,tau_cap)
        self.G=np.minimum(self.G,[gaussian_G(w.T,w.s) for w in self.columns])


def install_factory():
    import stable_cover
    if getattr(stable_cover,'_gaussian_installed',False):return
    original=stable_cover.weight_factory
    def factory(hi,N,kind,target=None):
        if kind!='gaussian':return original(hi,N,kind,target)
        cache={}
        def choose(tauhi):
            bucket=min(20,max(1,round(20*tauhi/hi)))
            if bucket not in cache:cache[bucket]=GaussianBatchWeights(hi,N,hi*bucket/20)
            return cache[bucket]
        return choose
    stable_cover.weight_factory=factory;stable_cover._gaussian_installed=True


def pilot():
    began=time.monotonic();checks=[];rows=[]
    for T,s in [(3.,.6),(5.7,.375),(10.,.3),(100.,.05),(400.,.0125)]:
        def integrand(u):
            z=math.pi*u
            h=z/3+z**3/45+2*z**5/945 if z<1e-3 else 1/z-1/math.tan(z)
            return (1-u)*math.sqrt(1+h*h)*math.exp(-T*T*u*u/2)
        ref,err=quad(integrand,0,s,epsabs=1e-13,epsrel=1e-12)
        ref+=float(exp1(T*T*s*s/2))/(2*math.pi)
        upper=gaussian_G(T,s)
        assert upper>=ref+abs(err),(T,s,upper,ref,err)
        checks.append({'T':T,'s':s,'upper':upper,'floating_reference':ref,'error_estimate':err})
    for L,d,b,tau in [(.4,.195600146724,.092507500581,.376),
        (.5945776691627666,.3543870444748948,.2470204829050391,.5203934540269098),
        (.6,.397750191359,.280850860934,.54),(.7,.454995299361,.359409523445,.595)]:
        cap=L*min(20,max(1,round(20*tau/L)))/20
        for N in [512,1024]:
            w=GaussianBatchWeights(L,N,cap);e=2e-12;box=(d-e,d+e,b-e,b+e,tau-e,tau+e)
            r=w.envelopes(L-e,L+e,box);vals=w.integrals(r[0],r[1],L-e,r[2])
            new_G=w.G.copy();w.G=np.array([col.G for col in w.columns])
            old=w.integrals(r[0],r[1],L-e,r[2])
            rows.append({'L':L,'N':N,'old_G':w.G.tolist(),'new_G':new_G.tolist(),
                'old_columns':old.tolist(),'new_columns':vals.tolist(),'upper':float(np.min(vals))})
    return {'status':'Gaussian coefficient certificates and point bounds only',
        'all_passed':True,'integral_checks':checks,'rows':rows,
        'elapsed_seconds':time.monotonic()-began}
