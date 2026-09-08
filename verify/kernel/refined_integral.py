"""Positive Fubini integration and freely chosen smoothing parameters."""
import math
import time
import numpy as np
from flint import arb
from interval_bounds import (PI,U,up,al,au,positive_dot,regularized_kernel,two_kernel)
from directional_bounds import DirectionalWeights,directional_checks


class RefinedWeights(DirectionalWeights):
    def __init__(self,Lhi,N=512,tau_cap=None,scale=1.,shift=0.):
        self.N=N;self.log_correction=True
        cap=Lhi if tau_cap is None else tau_cap
        self.T=float(round(scale*2*math.pi/(Lhi+cap),8))
        self.s=float(min(.14+.43*Lhi+shift,5/self.T))
        assert self.T>0 and 0<self.s<1
        T,s=arb(self.T),arb(self.s)
        ts=[T*s*j/N for j in range(N+1)]
        self.t=(np.array([al(x) for x in ts]),np.array([au(x) for x in ts]))
        self.t[0][0]=self.t[1][0]=0.
        self.tlo=self.t[0][:-1];self.thi=self.t[1][1:]
        self.gauss=np.array([au((-x*x/2).exp()) for x in ts[:-1]])
        primitive=[(x*(x*x/2).exp()-(PI/2).sqrt()*(x/arb(2).sqrt()).erfi())/2 for x in ts]
        self.inner=np.array([au(primitive[j+1]-primitive[j]) for j in range(N)])
        e1=[None]+[(x*x/2).expint(1) for x in ts[1:]]
        self.low=np.zeros(N);self.outer=np.zeros(N);self.local=np.zeros(N)
        for j in range(N):
            a,b=ts[j:j+2]
            A=regularized_kernel((s*j/N).union(s*(j+1)/N))/PI
            if j==0:self.local[j]=au(A*b**3/18)
            else:
                logarithm=(b/a).log()
                self.low[j]=au(A*logarithm)
                self.outer[j]=au(A*(e1[j]-e1[j+1])/2)
                self.local[j]=au(A*((b**3-a**3)/18-a**3*logarithm/6))
        us=[s+(1-s)*j/N for j in range(N+1)]
        high_ts=[T*u for u in us]
        self.ht=(np.array([al(x) for x in high_ts]),np.array([au(x) for x in high_ts]))
        self.high=np.array([au(two_kernel(us[j].union(us[j+1]))*(us[j+1]-us[j]))
                            for j in range(N)])
        z=T*s/arb(2).sqrt();e=(-z*z).exp();er=z.erf()
        i0=(PI/2).sqrt()/T*er;i1=(1-e)/(T*T)
        i2=(PI/2).sqrt()/T**3*er-s*e/(T*T)
        self.G=au(i0-i1+PI*PI/18*i2+(z*z).expint(1)/(2*PI))
        assert all(np.all(np.isfinite(x)) and np.all(x>=0) for x in
                   [self.inner,self.low,self.outer,self.local,self.high])

    def first_term(self,B,Llo,lowcf=None):
        cumul=up(np.cumsum(up(B*self.inner))*(1+4*self.N*U))
        previous=np.r_[0.,cumul[:-1]]
        contributions=up(up(self.outer*previous)+up(self.local*B))
        if lowcf is not None:
            cap=up(up(up(self.low*lowcf)+self.outer)/Llo)
            contributions[1:]=np.minimum(contributions[1:],cap[1:])
        return positive_dot(np.ones(self.N),contributions)

    def integrate(self,B,fcf,Llo,lowcf=None):
        i1=self.first_term(B,Llo,lowcf);i2=positive_dot(self.high,fcf)
        return float(up(i1+up(up(i2+self.G)/Llo)))


def integration_checks():
    from scipy.integrate import quad
    from scipy.special import erfi
    from global_probe import kernel
    rows=[]
    for L in [.1,.4,.8]:
        w=RefinedWeights(L,256)
        for kind in ['one','Gaussian']:
            B=np.ones(w.N) if kind=='one' else w.gauss
            def integrand(u):
                t=w.T*u
                if t==0:return 0.
                R=(t-math.sqrt(math.pi/2)*math.exp(-t*t/2)*erfi(t/math.sqrt(2)))/2 if kind=='one' else t**3*math.exp(-t*t/2)/6
                return 2*abs(kernel(np.array([u]))[0])*R
            ref,err=quad(integrand,0,w.s,epsabs=1e-12,limit=150)
            upper=w.first_term(B,L)
            assert upper>=ref+abs(err), (L,kind,upper,ref,err)
            rows.append({'L':L,'inner_envelope':kind,'upper':upper,
                         'floating_reference':ref,'reference_error_estimate':err})
    return {'status':'numerical cross-check of two analytically known ODE solutions',
            'all_passed':True,'rows':rows}


def pilot():
    began=time.monotonic();checks=directional_checks();integration=integration_checks();rows=[]
    points=[(.3,.122163285004,.042698346158,.27),(.4,.18337994274,.088528580873,.36),
            (.5,.257461449454,.143137636601,.45),(.55,.299414822466,.191336229373,.495),
            (.6,.397750191359,.280850860934,.54),(.7,.454995299361,.359409523445,.595)]
    for L,d,b,tau in points:
        bucket=min(20,max(1,math.ceil(20*tau/L)))
        for N in [128,512,1024]:
            w=RefinedWeights(L,N,L*bucket/20);e=2e-12
            rows.append({'L':L,'d':d,'b':b,'tau':tau,'bucket':bucket,'N':N,
                'T':w.T,'s':w.s,'parameter_halfwidth':e,
                'upper':w.box(L-e,L+e,d-e,d+e,b-e,b+e,tau-e,tau+e)})
    return {'status':'narrow-box enclosures only; whole cover still required',
            'checks':checks,'integration_checks':integration,'rows':rows,
            'elapsed_seconds':time.monotonic()-began}
