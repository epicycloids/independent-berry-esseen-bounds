"""Contact inequalities for two triples with an asymmetric collision and negative parameter."""
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import importlib.util
from pathlib import Path
from flint import arb


@lru_cache(maxsize=1)
def arithmetic():
    path=Path(__file__).with_name('collision210_raw_replay_evaluator_v1.py')
    assert hashlib.sha256(path.read_bytes()).hexdigest()=='b619e499eee526a9166afc97ac1c258a1377804fa50256e0cce0e1833354a8fd'
    spec=importlib.util.spec_from_file_location('independent_two_triple_interval_arithmetic',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def low_potential(x,a):
    return (x+a)**2*(2*a-arb(3)/2-x)+2*x.max(arb(0))**3


def original_hessian(case,state):
    """Original fixed-moment Lagrangian Hessian, before flow elimination."""
    g=state.get
    u,h,v,A1,A2=(g(name) for name in ('u','h','v','A1','A2'))
    inv=lambda name:state.inv(name)
    Z,r1,l2,m2=(g(name) for name in ('Z','r1','l2','m2'))
    Hc,Hb,Ht,Hd=(g(name) for name in ('Hc','Hb','Ht','Hd'))
    r_eps=-Z*inv('A1')*inv('h')
    r_eta=-r1*(A1+h)*inv('A1')*inv('h')
    m2_eps=m2*(v-h)*inv('h')*inv('v')
    m2_eta=l2*A2*inv('h')*inv('v')
    l2_eps=-m2*v*inv('h')*inv('A2')
    l2_eta=-l2*(h+A2)*inv('h')*inv('A2')
    if case=='negative':
        first_weighted_curvature=-4*Z
        first_middle_derivative_times_slope=-Z*(h-u)*inv('h')
    else:
        first_weighted_curvature=Z*(6*g('c')-3)*inv('u')
        first_middle_derivative_times_slope=-Z*(h-u)*g('Jc')*inv('u')*inv('h')
    first_eta_middle_times_slope=r1*A1*g('Jc')*inv('h')
    J11=-2*g('G')*r_eps*m2_eps-first_weighted_curvature-m2*(6*abs(g('d'))-3)
    J11+=2*first_middle_derivative_times_slope-2*m2_eps*Hd
    J22=-2*g('G')*r_eta*m2_eta-r1*(6*g('b')-3)-l2*(6*g('t')-3)
    J22+=2*r_eta*Hb-2*l2_eta*Ht
    J12=-g('G')*(r_eps*m2_eta+r_eta*m2_eps)+first_eta_middle_times_slope+r_eps*Hb
    J12-=l2_eps*Ht+m2_eta*Hd
    return J11,J22,J12


def support_hessian(case,state):
    J11,J22,J12=original_hessian(case,state)
    state.meet('minus_J11',-J11,lo=0)
    state.meet('minus_J22',-J22,lo=0)
    state.meet('J12',J12)
    state.meet('support_determinant',state.get('minus_J11')*state.get('minus_J22')-state.get('J12')**2,lo=0)



def negative_geometry(state):
    """Separate implementation of the accepted analytic location lemma."""
    q=arithmetic().q
    put,get=state.meet,state.get
    put('d',hi=Q(-1,2))
    put('v',lo=arb(3)/2)
    put('h',hi=get('v'))
    put('v',lo=get('h'))
    put('h',hi=1+arb(2).sqrt()/2-q(Q(2,3))*get('u'))
    put('u2',get('u')**2,lo=0,hi=Q(9,16))
    put('u2',lo=(2*get('e')-1)**2-2)
    put('e',hi=arb(1)/2+(arb(1)/2+get('u2')/4).sqrt())
    put('u',lo=state.power('u2','sqrt'))

def evaluate(case,box,*,sweeps=10,snapshot=False):
    if case not in ('negative','positive') or not 1<=sweeps<=16:
        raise ValueError('invalid independent TT configuration')
    raw=arithmetic()
    q=raw.q
    x=[q(Q(n,1<<depth)).union(q(Q(n+1,1<<depth))) for n,depth in box]
    if len(x)!=4:raise ValueError('original four-cube box required')
    state=raw.Candidate()
    put,get=state.meet,state.get
    inv=state.inv
    sqrt=lambda name:state.power(name,'sqrt')
    pow3=lambda name:state.power(name,'three_halves')
    phi0=1/(2*arb.pi()).sqrt()
    dmin=q(Q(49,120))/arb(2).sqrt()
    try:
        k=put('k',x[0],lo=0,hi=1)
        if case=='negative':
            a=put('a',(k+2)/4,lo=Q(1,2),hi=Q(3,4))
            c=put('c',(k-1)/2,lo=Q(-1,2),hi=0)
            put('u',3*k/4,lo=0,hi=Q(3,4))
        else:
            a=put('a',3*(1+k)**2/(2*(k**3+3*k+2)),lo=Q(3,4),hi=1)
            c=put('c',k*a,lo=0,hi=1)
            put('u',a+c,lo=Q(3,4),hi=2)
        B=put('B',(1+(1+4*a*(1-a)).sqrt())/2,lo=1,hi=Q(5,4))
        put('w',x[1],lo=0,hi=1)
        put('b',x[1]*B,lo=Q(13,512),hi=Q(5,4))
        put('e',1+3*x[2]/2,lo=1,hi=Q(5,2))
        put('t',7*x[3]/2,lo=Q(13,1280),hi=Q(7,2))
        put('h',get('b')-c,lo=Q(25,672),hi=2)
        put('d',get('h')-get('t'),lo=-4,hi=2)
        put('v',get('e')-get('d'),lo=Q(25,912),hi=6)
        put('A1',a+get('b'),lo=Q(1,2),hi=Q(9,4))
        put('A2',get('t')+get('e'),lo=1,hi=6)
        put('gap',get('A2')-get('A1'),lo=0,hi=6)
        put('o',get('b')-get('t'),lo=-4,hi=Q(5,4))
        for name in ('q1','r1','l2','m2','r2'):put(name,lo=0,hi=1)
        for name in ('Z','Y'):put(name,lo=0,hi=get('u'))
        for name in ('V1','V2'):put(name,lo=Q(13,512),hi=4)
        put('V',lo=Q(2401,5476),hi=4)
        for name in ('R1','R2','R'):put(name,lo=0,hi=11)
        for name in ('M1','M2'):put(name,lo=-4,hi=4)
        put('G',lo=Q(1,2),hi=20)
        put('U',lo=0,hi=20)
        put('g',lo=0,hi=20)
        for name in ('Hb','Ht','Hd','Hc'):put(name,lo=0,hi=48)
        put('Jc',lo=0,hi=3)
        put('C',lo=Q(49,120),hi=Q(113,250))
        put('sigma',lo=Q(1,80),hi=3)
        put('D',lo=dmin,hi=Q(11,20))
        put('F',lo=dmin,hi=1)
        put('Phi',lo=0,hi=1-dmin)
        put('phi',lo=Q(9,50),hi=phi0)
        put('z',lo=Q(-5,4),hi=Q(3,5))
        put('rootV',sqrt('V'))
        for iteration in range(sweeps):
            if case=='negative':negative_geometry(state)
            for terms,right in (
                ({'A1':1,'a':-1,'b':-1},0),({'A1':1,'u':-1,'h':-1},0),
                ({'A2':1,'t':-1,'e':-1},0),({'A2':1,'h':-1,'v':-1},0),
                ({'h':1,'b':-1,'c':1},0),({'h':1,'t':-1,'d':-1},0),
                ({'v':1,'e':-1,'d':1},0),({'u':1,'a':-1,'c':-1},0),
                ({'gap':1,'A2':-1,'A1':1},0),({'gap':1,'v':-1,'u':1},0),
                ({'o':1,'b':-1,'t':1},0),({'o':1,'c':-1,'d':-1},0),
                ({'q1':1,'r1':1},1),({'l2':1,'m2':1,'r2':1},1),
                ({'V':1,'V1':-1,'V2':-1},0),({'R':1,'R1':-1,'R2':-1},0),
                ({'F':1,'D':-1,'Phi':-1},0)):
                state.linear(terms,right)
            if case=='negative':
                for terms,right in (({'a':3,'u':-1},Q(3,2)),({'c':3,'u':-2},Q(-3,2)),
                                    ({'u':4,'k':-3},0)):
                    state.linear(terms,right)
            else:
                put('k',get('c')*inv('a'))
                k=get('k')
                put('a',3*(1+k)**2/(2*(k**3+3*k+2)))
                put('c',k*get('a'))
            a,b,c,t,d,e,u,h,v,A1,A2=(get(name) for name in ('a','b','c','t','d','e','u','h','v','A1','A2'))
            ba=put('beta_a',a*(1-a),lo=0,hi=Q(1,4))
            be=put('beta_e',e*(e-1),lo=0,hi=Q(15,4))
            put('disc_a',1-4*ba,lo=0,hi=1)
            put('a',(1+sqrt('disc_a'))/2)
            put('e',(1+(1+4*be).sqrt())/2)
            B=put('B',(1+(1+4*ba).sqrt())/2)
            put('b',get('w')*B)
            put('w',get('b')*inv('B'))
            put('G',raw.potential(-t,e))
            put('g',raw.potential(d,e))
            put('U',-low_potential(b,a))
            put('q1',get('g')*inv('G'))
            put('m2',get('U')*inv('G'))
            put('g',get('q1')*get('G'))
            put('U',get('m2')*get('G'))
            put('Z',A1*get('q1')-b,lo=0,hi=u)
            put('Y',b-h*get('q1'),lo=0,hi=u)
            put('q1',(get('Z')+b)*inv('A1'))
            put('r1',(a-get('Z'))*inv('A1'))
            put('Z',u*get('q1')-get('Y'))
            put('q1',(b-get('Y'))*inv('h'))
            put('V1',a*b-h*get('Z'))
            put('Z',(a*b-get('V1'))*inv('h'))
            put('l2',(get('V2')+e*d)*inv('A2')*inv('h'))
            put('m2',(t*e-get('V2'))*inv('h')*inv('v'))
            put('r2',(get('V2')-t*d)*inv('A2')*inv('v'))
            put('V2',get('l2')*t*t+get('m2')*d*d+get('r2')*e*e)
            put('m2',(e-get('l2')*A2)*inv('v'))
            put('m2',(t-get('r2')*A2)*inv('h'))
            put('V2',t*e-get('m2')*h*v)
            if case=='negative':
                put('V1',a*a*get('q1')+get('Z')*(c-a)+get('r1')*b*b)
                put('R1',a**3*get('q1')-get('Z')*(a*a-a*c+c*c)+get('r1')*b**3)
                put('M1',-get('V1')+2*get('r1')*b*b)
                gamma1=a-b-c+2*b**3*inv('A1')*inv('h')
                beta1=-1+2*b*b*inv('A1')*inv('h')
                put('Jc',u)
            else:
                ell=put('l1',get('Y')*inv('u'),lo=0,hi=1)
                middle=put('m1',get('Z')*inv('u'),lo=0,hi=1)
                put('R1',ell*a**3+middle*c**3+get('r1')*b**3)
                put('V1',ell*a*a+middle*c*c+get('r1')*b*b)
                put('M1',-ell*a*a+middle*c*c+get('r1')*b*b)
                gamma1=-a+b+c+2*a**3*inv('A1')*inv('u')
                beta1=1-2*a*a*inv('A1')*inv('u')
                put('Jc',3*(ba+c-c*c)*inv('u'))
            gamma1=put('gamma1',gamma1,lo=A1/2,hi=2*A1)
            beta1=put('beta1',beta1,lo=-1,hi=1)
            endpointR1=a*b*(a*a+b*b)*inv('A1')
            endpointM1=a*b*(b-a)*inv('A1')
            put('R1',endpointR1+gamma1*(get('V1')-a*b))
            put('V1',a*b+(get('R1')-endpointR1)*inv('gamma1'))
            put('M1',endpointM1+beta1*(get('V1')-a*b))
            gamma2,beta2=raw.interpolate(t,d,e,h,v,A2)
            gamma2=put('gamma2',gamma2,lo=A2/2,hi=2*A2)
            beta2=put('beta2',beta2,lo=-1,hi=1)
            endpointR2=t*e*(t*t+e*e)*inv('A2')
            endpointM2=t*e*(e-t)*inv('A2')
            put('R2',get('l2')*t**3+get('m2')*abs(d)**3+get('r2')*e**3)
            put('M2',-get('l2')*t*t+get('m2')*d*abs(d)+get('r2')*e*e)
            put('R2',endpointR2+gamma2*(get('V2')-t*e))
            put('V2',t*e+(get('R2')-endpointR2)*inv('gamma2'))
            put('M2',endpointM2+beta2*(get('V2')-t*e))
            put('sigma',3*(ba-get('M1')))
            put('sigma',3*(be-get('M2')))
            put('beta_a',get('sigma')/3+get('M1'))
            put('beta_e',get('sigma')/3+get('M2'))
            put('M1',get('beta_a')-get('sigma')/3)
            put('M2',get('beta_e')-get('sigma')/3)
            put('Hc',3*(ba+c-c*abs(c)))
            put('Hb',3*(ba+b-b*b))
            put('Ht',3*(t*t-t+be))
            put('Hd',3*(be+d-d*abs(d)))
            put('Hc',u*get('Jc'))
            put('flow_middle',get('Z')*get('Jc'),lo=0)
            put('flow_middle',get('m2')*get('Hd'))
            put('flow_endpoint',get('r1')*get('Hb'),lo=0)
            put('flow_endpoint',get('l2')*get('Ht'))
            put('sigma',get('flow_middle')+get('flow_endpoint'))
            for mass,slope,flow in (('Z','Jc','flow_middle'),('m2','Hd','flow_middle'),
                                    ('r1','Hb','flow_endpoint'),('l2','Ht','flow_endpoint')):
                if state.low[slope]>0:put(mass,get(flow)*inv(slope))
            put('F',get('l2')+get('m2')*get('q1'))
            put('F',e*inv('A2')+get('m2')*(get('q1')-v*inv('A2')))
            put('V1',get('V')-get('V2'),lo=13*get('V')/200,hi=187*get('V')/200)
            put('V2',get('V')-get('V1'),lo=13*get('V')/200,hi=187*get('V')/200)
            put('V',get('V1')+get('V2'))
            put('R1',lo=pow3('V1'))
            put('R2',lo=pow3('V2'))
            put('R',get('R1')+get('R2'),lo=pow3('V')/arb(2).sqrt(),hi=q(Q(66,49))*pow3('V'))
            put('C',pow3('V')*inv('G'))
            put('CG',get('C')*get('G'),lo=0)
            put('V',state.power('CG','two_thirds'))
            put('G',pow3('V')*inv('C'))
            put('rootV',sqrt('V'))
            put('z',get('o')*inv('rootV'))
            phi,Phi=raw.gaussian(get('z'))
            put('phi',phi)
            put('Phi',Phi)
            put('D',get('F')-get('Phi'))
            put('D',get('R')*inv('G'))
            put('R',get('D')*get('G'))
            put('F',get('D')+get('Phi'))
            put('Phi',get('F')-get('D'))
            put('phi',get('sigma')*get('rootV')*inv('G'))
            put('sigma',get('phi')*get('G')*inv('rootV'))
            put('rootV',get('G')*get('phi')*inv('sigma'))
            put('V',get('rootV')**2)
            put('sigma',get('V')*get('phi')*inv('C'))
            put('V',get('sigma')*get('C')*inv('phi'))
            put('R',get('V')-get('sigma')*get('o')/3)
            put('V',get('R')+get('sigma')*get('o')/3)
            put('o',3*(get('V')-get('R'))*inv('sigma'))
            put('z2',-(2*arb.pi()*get('phi')**2).log(),lo=0)
            put('z2',get('z')**2)
            rootz=sqrt('z2')
            if state.high['z']<=0:put('z',-rootz)
            elif state.low['z']>=0:put('z',rootz)
            else:put('z',lo=-rootz,hi=rootz)
            put('o',get('z')*get('rootV'))
            if state.high['o']<=0:put('V',lo=Q(2401,4356))
            for span in ('A1','A2'):
                put('V',lo=q(Q(2401,24964))*get(span)**2)
                put(span,hi=q(Q(158,49))*sqrt('V'))
            width1=(b*u).min(a*h)
            width2=(e*h).min(t*v)
            put('V',hi=3*(width1+width2))
            jump1=low_potential(get('c')-get('v'),get('a'))
            put('r2',hi=jump1*inv('G'))
            jump2=raw.potential(get('u')+get('d'),get('e'))
            put('Y',hi=get('u')*jump2*inv('G'))
            support_hessian(case,state)
    except raw.Impossible as failure:
        return dict(reason='independent_TT_empty_'+failure.name,lower=failure.low.str(70),
                    upper=failure.high.str(70),method='independent_two_triple_raw_law',
                    sweeps_completed=locals().get('iteration',-1),intersections=state.steps)
    return dict(reason=None,method='independent_two_triple_raw_law',sweeps_completed=sweeps,
                intersections=state.steps,
                **({'candidate_bounds':{name:[state.low[name].str(70),state.high[name].str(70)] for name in state.low}} if snapshot else {}))
