"""Contact inequalities for a triple/pair lower-gap collision, evaluated from centered marginal laws."""
from fractions import Fraction

from flint import arb


def q(value):
    value=Fraction(value)
    return arb(value.numerator)/value.denominator


class Impossible(Exception):
    def __init__(self,name,low,high):
        self.name,self.low,self.high=name,low,high


class Candidate:
    def __init__(self):
        self.low={}
        self.high={}
        self.steps=0

    def meet(self,name,value=None,lo=None,hi=None):
        lows=[] if name not in self.low else [self.low[name]]
        highs=[] if name not in self.high else [self.high[name]]
        if value is not None:
            if not isinstance(value,arb) or not value.is_finite():
                raise ValueError('finite ordinary interval required')
            lows.append(value.lower())
            highs.append(value.upper())
        if lo is not None:lows.append((lo if isinstance(lo,arb) else q(lo)).lower())
        if hi is not None:highs.append((hi if isinstance(hi,arb) else q(hi)).upper())
        left,right=max(lows),min(highs)
        if left>right:raise Impossible(name,left,right)
        self.low[name],self.high[name]=left,right
        self.steps+=1
        return left.union(right)

    def get(self,name):return self.low[name].union(self.high[name])

    def inv(self,name):
        if not self.low[name]>0:raise ValueError('positive candidate denominator required')
        return (1/self.high[name]).union(1/self.low[name])

    def power(self,name,power):
        if self.low[name]<0:raise ValueError('nonnegative candidate power required')
        if power=='sqrt':return self.low[name].sqrt().union(self.high[name].sqrt())
        if power=='three_halves':return (self.low[name]*self.low[name].sqrt()).union(self.high[name]*self.high[name].sqrt())
        if power=='two_thirds':return (self.low[name].root(3)**2).union(self.high[name].root(3)**2)
        raise ValueError('unknown monotone power')

    def linear(self,terms,right=0):
        for name,coefficient in terms.items():
            value=q(right)
            for other,a in terms.items():
                if other!=name:value-=q(a)*self.get(other)
            self.meet(name,value/q(coefficient))


def gaussian(z):
    return (-z*z/2).exp()/(2*arb.pi()).sqrt(),(1+(z/arb(2).sqrt()).erf())/2


def potential(x,b):
    negative=(-x).max(arb(0))
    return (x-b)**2*(x+2*b-q(Fraction(3,2)))+2*negative**3


def interpolate(t,c,b,h=None,v=None,A=None):
    if h is None:h=t+c
    if v is None:v=b-c
    if A is None:A=t+b
    cubic=t**3/h/A-abs(c)**3/h/v+b**3/A/v
    signed=-t*t/h/A-c*abs(c)/h/v+b*b/A/v
    return cubic,signed


def evaluate(case,box,*,sweeps=10,snapshot=False):
    if case not in ('negative','positive') or not 1<=sweeps<=16:
        raise ValueError('invalid fixed replay configuration')
    x=[q(Fraction(n,1<<d)).union(q(Fraction(n+1,1<<d))) for n,d in box]
    if len(x)!=3:raise ValueError('original three-cube box required')
    state=Candidate()
    put,get=state.meet,state.get
    phi0=1/(2*arb.pi()).sqrt()
    dmin=q(Fraction(49,120))/arb(2).sqrt()
    try:
        b=put('b',1+3*x[0]/2,lo=1,hi=Fraction(5,2))
        t=put('t',7*x[1]/2,lo=Fraction(13,1280),hi=Fraction(7,2))
        c=put('c',-t*(1-x[2]) if case=='negative' else b*x[2])
        put('A',t+b,lo=1,hi=6)
        put('h',t+c,lo=q(Fraction(13,128)).sqrt(),hi=6)
        put('v',b-c,lo=Fraction(25,672),hi=6)
        for name in ('p','q'):
            put(name,lo=Fraction(13,18432),hi=1)
        for name in ('ell','m','r'):
            put(name,lo=0,hi=1)
        for name in ('d','e'):
            put(name,lo=q(Fraction(13,18432))*q(Fraction(13,128)).sqrt(),hi=6)
        put('VT',lo=Fraction(13,512),hi=4)
        put('VB',lo=Fraction(13,512),hi=4)
        put('V',lo=Fraction(2401,4356 if case=='negative' else 5476),hi=4)
        for name in ('RT','RB','R'):put(name,lo=0,hi=11)
        put('G',lo=Fraction(1,2),hi=20)
        put('C',lo=Fraction(49,120),hi=Fraction(227,500))
        put('sigma',lo=Fraction(1,80),hi=48)
        put('D',lo=dmin,hi=Fraction(11,20))
        put('F',lo=dmin,hi=1)
        put('Phi',lo=0,hi=1-dmin)
        put('phi',lo=Fraction(9,50),hi=phi0)
        put('z',lo=Fraction(-5,4),hi=0 if case=='negative' else Fraction(3,5))
        put('o',lo=-4,hi=0 if case=='negative' else 4)
        put('rootV',state.power('V','sqrt'))
        for iteration in range(sweeps):
            for row,right in (
                ({'A':1,'b':-1,'t':-1},0),({'A':1,'h':-1,'v':-1},0),
                ({'h':1,'t':-1,'c':-1},0),({'v':1,'b':-1,'c':1},0),
                ({'h':1,'d':-1,'e':-1},0),({'o':1,'c':-1,'d':1},0),
                ({'o':1,'e':-1,'t':1},0),({'p':1,'q':1},1),
                ({'ell':1,'m':1,'r':1},1),({'V':1,'VT':-1,'VB':-1},0),
                ({'R':1,'RT':-1,'RB':-1},0),({'F':1,'D':-1,'Phi':-1},0)):
                state.linear(row,right)
            b,t,c,A,h,v=(get(key) for key in ('b','t','c','A','h','v'))
            beta=put('beta',b*(b-1),lo=0,hi=Fraction(15,4))
            put('b',(1+(1+4*beta).sqrt())/2)
            G=put('G',potential(-t,b))
            N=put('N',potential(c,b),lo=0,hi=20)
            G=put('G',t**3-3*t*t/2+3*beta*t+b*b*(2*b-q(Fraction(3,2))))
            put('q',N*state.inv('G'))
            put('p',1-get('q'))
            put('N',get('q')*G)
            put('d',get('p')*h)
            put('e',get('q')*h)
            d,e=get('d'),get('e')
            put('VB',d*e)
            put('d',get('VB')*state.inv('e'))
            put('e',get('VB')*state.inv('d'))
            put('p',d*state.inv('h'))
            put('q',e*state.inv('h'))
            kB=put('kB',h-2*get('VB')*state.inv('h'),lo=h/2,hi=h)
            put('RB',get('VB')*kB)
            put('VB',get('RB')*state.inv('kB'))
            Ht=put('Ht',3*(t*t-t+beta),lo=0,hi=48)
            Hm=put('Hm',3*(beta+c-c*abs(c)),lo=0,hi=48)
            gamma,betaDD=interpolate(t,c,b,h,v,A)
            gamma=put('gamma',gamma,lo=A/2,hi=2*A)
            betaDD=put('betaDD',betaDD,lo=-1,hi=1)
            endpointR=t*b*(t*t+b*b)*state.inv('A')
            endpointM=t*b*(b-t)*state.inv('A')
            put('VT',t*b-get('m')*h*v)
            put('ell',(get('VT')+b*c)*state.inv('A')*state.inv('h'))
            put('m',(t*b-get('VT'))*state.inv('h')*state.inv('v'))
            put('r',(get('VT')-t*c)*state.inv('A')*state.inv('v'))
            put('VT',get('ell')*t*t+get('m')*c*c+get('r')*b*b)
            RT=put('RT',endpointR+gamma*(get('VT')-t*b))
            put('RT',get('ell')*t**3+get('m')*abs(c)**3+get('r')*b**3)
            put('VT',t*b+(RT-endpointR)*state.inv('gamma'))
            signedM=endpointM+betaDD*(get('VT')-t*b)
            put('sigma',3*(beta-signedM))
            put('sigma',get('ell')*Ht+get('m')*Hm)
            put('beta',get('sigma')/3+signedM)
            sigma=get('sigma')
            put('F',get('ell')+get('m')*get('q'))
            put('F',b*state.inv('A')+get('m')*(get('q')-v*state.inv('A')))
            j=get('e')-get('d')
            J=j*(j*j-3*h/2)
            put('m',(sigma*h-J)*state.inv('G'))
            put('sigma',(get('m')*G+J)*state.inv('h'))
            gapflow=3*get('p')*get('q')*h*(1-kB)
            put('sigma',(get('ell')*Ht-gapflow)*state.inv('p'))
            put('sigma',(get('m')*Hm+gapflow)*state.inv('q'))
            if state.low['Ht']>0:put('ell',(get('p')*get('sigma')+gapflow)*state.inv('Ht'))
            if state.low['Hm']>0:put('m',(get('q')*get('sigma')-gapflow)*state.inv('Hm'))
            put('m',(b-get('ell')*A)*state.inv('v'))
            put('m',(t-get('r')*A)*state.inv('h'))
            put('VT',get('V')-get('VB'),lo=13*get('V')/200,hi=187*get('V')/200)
            put('VB',get('V')-get('VT'),lo=13*get('V')/200,hi=187*get('V')/200)
            put('V',get('VT')+get('VB'))
            put('RT',lo=state.power('VT','three_halves'))
            put('RB',lo=state.power('VB','three_halves'))
            put('R',get('RT')+get('RB'),lo=state.power('V','three_halves')/arb(2).sqrt(),hi=q(Fraction(66,49))*state.power('V','three_halves'))
            put('C',state.power('V','three_halves')*state.inv('G'))
            put('CG',get('C')*get('G'),lo=0)
            put('V',state.power('CG','two_thirds'))
            put('G',state.power('V','three_halves')*state.inv('C'))
            put('rootV',state.power('V','sqrt'))
            put('o',get('c')-get('d'))
            put('z',get('o')*state.inv('rootV'))
            phi,Phi=gaussian(get('z'))
            put('phi',phi)
            put('Phi',Phi)
            put('D',get('F')-get('Phi'))
            put('D',get('R')*state.inv('G'))
            put('R',get('D')*get('G'))
            put('F',get('D')+get('Phi'))
            put('Phi',get('F')-get('D'))
            put('phi',get('sigma')*get('rootV')*state.inv('G'))
            put('sigma',get('phi')*get('G')*state.inv('rootV'))
            put('rootV',get('G')*get('phi')*state.inv('sigma'))
            put('V',get('rootV')**2)
            put('R',get('V')-get('sigma')*get('o')/3)
            put('V',get('R')+get('sigma')*get('o')/3)
            put('o',3*(get('V')-get('R'))*state.inv('sigma'))
            put('z2',-(2*arb.pi()*get('phi')**2).log(),lo=0)
            put('z2',get('z')**2)
            rootz2=state.power('z2','sqrt')
            if state.high['z']<=0:put('z',-rootz2)
            elif state.low['z']>=0:put('z',rootz2)
            else:put('z',lo=-rootz2,hi=rootz2)
            put('o',get('z')*get('rootV'))
            put('b',lo=get('V')*state.inv('h')/3)
            put('h',lo=get('V')*state.inv('b')/3)
            put('t',lo=get('V')*state.inv('v')/3)
            put('v',lo=get('V')*state.inv('t')/3)
            put('h',lo=2*state.power('VB','sqrt'))
    except Impossible as failure:
        return dict(reason='independent_empty_'+failure.name,lower=failure.low.str(70),
                    upper=failure.high.str(70),method='independent_raw_law_candidate',
                    sweeps_completed=locals().get('iteration',-1),intersections=state.steps)
    return dict(reason=None,method='independent_raw_law_candidate',sweeps_completed=sweeps,
                intersections=state.steps,
                **({'candidate_bounds':{name:[state.low[name].str(70),state.high[name].str(70)] for name in state.low}} if snapshot else {}))
