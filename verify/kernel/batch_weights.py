"""The same outward cosine envelopes, with smoothing choices as columns."""
import time
import numpy as np
from interval_bounds import up,U,positive_dot
from paired_bounds import CosineWeights
from stable_bounds import finite_law_checks

CHOICES=[(1.,0.),(1.02,-.02),(1.04,-.02)]


class BatchWeights(CosineWeights):
    def __init__(self,Lhi,N=512,tau_cap=None):
        self.columns=[CosineWeights(Lhi,N,tau_cap,scale,shift) for scale,shift in CHOICES]
        self.N=N;self.T=self.columns[0].T;self.s=self.columns[0].s;self.log_correction=True
        for name in ['t','ht']:
            setattr(self,name,tuple(np.column_stack([getattr(w,name)[j] for w in self.columns]) for j in [0,1]))
        for name in ['tlo','thi','gauss','inner','low','outer','local','high']:
            setattr(self,name,np.column_stack([getattr(w,name) for w in self.columns]))
        self.G=np.array([w.G for w in self.columns])

    def first_terms(self,B,Llo,lowcf=None):
        cumul=up(np.cumsum(up(B*self.inner),axis=0)*(1+4*self.N*U))
        previous=np.vstack([np.zeros((1,len(CHOICES))),cumul[:-1,:]])
        contributions=up(up(self.outer*previous)+up(self.local*B))
        if lowcf is not None:
            cap=up(up(up(self.low*lowcf)+self.outer)/Llo)
            contributions[1:,:]=np.minimum(contributions[1:,:],cap[1:,:])
        return np.array([positive_dot(np.ones(self.N),contributions[:,j]) for j in range(len(CHOICES))])

    def integrals(self,B,F,Llo,lowcf=None):
        i1=self.first_terms(B,Llo,lowcf)
        i2=np.array([positive_dot(self.high[:,j],F[:,j]) for j in range(len(CHOICES))])
        return up(i1+up(up(i2+self.G)/Llo))

    def integrate(self,B,F,Llo,lowcf=None):
        return float(np.min(self.integrals(B,F,Llo,lowcf)))


def install_factory():
    import stable_cover
    if getattr(stable_cover,'_batch_installed',False):return
    original=stable_cover.weight_factory
    def factory(hi,N,kind,target=None):
        if kind!='batch':return original(hi,N,kind,target)
        cache={}
        def choose(tauhi):
            bucket=min(20,max(1,round(20*tauhi/hi)))
            if bucket not in cache:cache[bucket]=BatchWeights(hi,N,hi*bucket/20)
            return cache[bucket]
        return choose
    stable_cover.weight_factory=factory;stable_cover._batch_installed=True


def pilot():
    began=time.monotonic();rows=[];count=0;max_difference=0.
    points=[(.4,.195600146724,.092507500581,.376),(.5,.257461449454,.155637636600,.45),
            (.5945776691627666,.3543870444748948,.2470204829050391,.5203934540269098),
            (.6,.397750191359,.280850860934,.54),(.7,.454995299361,.359409523445,.595)]
    for L,d,b,tau in points:
        for N in [512,1024]:
            cap=L*min(20,max(1,round(20*tau/L)))/20;w=BatchWeights(L,N,cap)
            for e in [2e-12,.001]:
                box=(d-e,d+e,b-e,b+e,tau-e,tau+e)
                result=w.envelopes(L-e,L+e,box);assert result is not None
                vals=w.integrals(result[0],result[1],L-e,result[2])
                for j,col in enumerate(w.columns):
                    ref=col.envelopes(L-e,L+e,box);assert ref is not None
                    assert all(np.array_equal(a[:,j],z) for a,z in zip(result,ref))
                    old=col.integrate(ref[0],ref[1],L-e,ref[2])
                    error=abs(float(vals[j])-old);max_difference=max(max_difference,error)
                    assert error<2e-14,(L,N,e,j,error)
                    count+=1
                rows.append({'L':L,'d':d,'b':b,'tau':tau,'N':N,'halfwidth':e,
                             'columns':vals.tolist(),'upper':float(np.min(vals))})
    # The law harness uses one frequency column at a time.
    laws=[]
    for column in range(len(CHOICES)):
        class View:
            def __init__(self,Lhi,N):
                self.w=BatchWeights(Lhi,N);self.tlo=self.w.tlo[:,column];self.thi=self.w.thi[:,column]
                self.ht=(self.w.ht[0][:,column],self.w.ht[1][:,column])
            def envelopes(self,*args):
                r=self.w.envelopes(*args)
                return None if r is None else tuple(a[:,column] for a in r)
        laws.append(finite_law_checks(View,True))
    w=BatchWeights(.6,1024,.54);box=(.395,.4,.278,.283,.538,.542)
    start=time.monotonic()
    for _ in range(20):w.box(.595,.6,*box)
    batch_seconds=time.monotonic()-start;start=time.monotonic()
    for _ in range(20):
        for col in w.columns:col.box(.595,.6,*box)
    separate_seconds=time.monotonic()-start
    return {'status':'column equivalence and finite-law checks only','all_passed':True,
        'column_comparisons':count,'maximum_integral_difference':max_difference,'finite_laws':laws,
        'rows':rows,'timing_20_boxes':{'batch':batch_seconds,'separate':separate_seconds},
        'elapsed_seconds':time.monotonic()-began}
