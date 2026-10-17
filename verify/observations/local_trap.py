#!/usr/bin/env python3
"""Verify the stored unequal-gap four-Bernoulli example by interval arithmetic.

The probability box is fixed: radius 1e-28, precision 320 bits. The analytic
stationary equations, complete probability Hessian, and coordinate cubic
majorants are evaluated without root finding or a candidate search. See
../../docs/numerical-observations.md for the local binary-law implication.
"""

from __future__ import annotations

import argparse
import copy
from decimal import Decimal
import json
import math
from pathlib import Path

from flint import arb, ctx


class Dual:
    def __init__(self,v,d=(0,0)):
        self.v = v if isinstance(v,arb) else arb(v)
        self.d = tuple(q if isinstance(q,arb) else arb(q) for q in d)
    def coerce(self,a):
        return a if isinstance(a,Dual) else Dual(a)
    def __add__(self,a):
        a=self.coerce(a); return Dual(self.v+a.v,[x+y for x,y in zip(self.d,a.d)])
    __radd__=__add__
    def __neg__(self):
        return Dual(-self.v,[-x for x in self.d])
    def __sub__(self,a):
        return self+-self.coerce(a)
    def __rsub__(self,a):
        return self.coerce(a)+-self
    def __mul__(self,a):
        a=self.coerce(a); return Dual(self.v*a.v,[x*a.v+self.v*y for x,y in zip(self.d,a.d)])
    __rmul__=__mul__
    def __truediv__(self,a):
        a=self.coerce(a); return Dual(self.v/a.v,[(x*a.v-self.v*y)/(a.v*a.v) for x,y in zip(self.d,a.d)])
    def __rtruediv__(self,a):
        return self.coerce(a)/self
    def __pow__(self,n):
        if not isinstance(n,int):
            raise TypeError("Dual supports integer powers")
        return Dual(self.v**n,[n*self.v**(n-1)*x for x in self.d]) if n else Dual(1)
    def sqrt(self):
        v=self.v.sqrt(); return Dual(v,[x/(2*v) for x in self.d])
    def exp(self):
        v=self.v.exp(); return Dual(v,[v*x for x in self.d])
    def erf(self):
        factor=2*(-self.v*self.v).exp()/arb.pi().sqrt()
        return Dual(self.v.erf(),[factor*x for x in self.d])


def expression(x,y,sqrt,exp,erf,pi):
    qx,qy=1-x,1-y
    V=3*x*qx+4*y*qy
    R=3*x*qx*(x*x+qx*qx)+8*y*qy*(y*y+qy*qy)
    o=2-3*x-2*y
    z=o/sqrt(V)
    phi=exp(-z*z/2)/sqrt(2*pi)
    w=phi/sqrt(V)
    F=(1-y)*(1-x**3)+y*qx**3
    D=F-(1+erf(z/sqrt(2)))/2
    lam=D/R
    A=(3*D+o*w)/(2*V)
    ga=-((1-y)*x*x+y*qx*qx)
    gb=-3*x*qx
    kx,ky=1-2*x,1-2*y
    G=(ga+w+A*kx-lam*kx**3,gb+2*w+4*A*ky-8*lam*ky**3)
    return dict(x=x,y=y,V=V,R=R,o=o,z=z,phi=phi,w=w,F=F,D=D,lam=lam,A=A,
                C=D*V*sqrt(V)/R,ga=ga,gb=gb,kx=kx,ky=ky,G=G)


def arb_state(x,y):
    return expression(x,y,lambda v:v.sqrt() if hasattr(v,"sqrt") else arb(v).sqrt(),
                      lambda v:v.exp(),lambda v:v.erf(),arb.pi())


def generic_pmf(probabilities,weights):
    one=probabilities[0]*0+1
    out=[one]
    for p,w in zip(probabilities,weights):
        nxt=[one*0 for _ in range(len(out)+w)]
        for j,v in enumerate(out):
            nxt[j]+=v*(1-p); nxt[j+w]+=v*p
        out=nxt
    return out


def generic_cdf(pmf,k):
    return sum(pmf[:max(0,min(len(pmf),k+1))],pmf[0]*0)


def krawczyk(center,inverse,radius):
    thin=[arb(v) for v in center]
    box=[arb(v,radius) for v in center]
    dual=arb_state(Dual(box[0],(1,0)),Dual(box[1],(0,1)))
    J=[[dual["G"][i].d[j] for j in range(2)] for i in range(2)]
    G=arb_state(*thin)["G"]
    B=[[arb(inverse[i][j]) for j in range(2)] for i in range(2)]
    Bdet=B[0][0]*B[1][1]-B[0][1]*B[1][0]
    residual=[-sum((B[i][j]*G[j] for j in range(2)),arb(0)) for i in range(2)]
    M=[[arb(int(i==j))-sum((B[i][k]*J[k][j] for k in range(2)),arb(0)) for j in range(2)] for i in range(2)]
    rad=arb(radius)
    bounds=[abs(residual[i])+sum((abs(M[i][j])*rad for j in range(2)),arb(0)) for i in range(2)]
    contraction=max(float(sum((abs(M[i][j]) for j in range(2)),arb(0))) for i in range(2))
    return box,dict(inclusion=bool(all(v<rad for v in bounds)),
        preconditioner_nonsingular=bool(Bdet>0 or Bdet<0),
        displacement_bound_balls=[str(v) for v in bounds],box_radius=radius,
        derivative_defect_balls=[[str(v) for v in row] for row in M],
        contraction_bound_approx=contraction,contraction_certified=bool(all(sum((abs(v) for v in row),arb(0))<1 for row in M)))


def certificate(box):
    st=arb_state(*box)
    x,y=box; pqx,pqy=x*(1-x),y*(1-y)
    feasible=bool(x>0 and x<1 and y>0 and y<1 and st["V"]>0 and st["R"]>0 and st["D"]>0 and st["lam"]>0)
    shape_x,shape_y=x*x+(1-x)**2,y*y+(1-y)**2
    bx=x+pqx*(2*st["A"]-3*st["lam"]*shape_x)/st["w"]
    by=y+2*pqy*(2*st["A"]-6*st["lam"]*shape_y)/st["w"]
    coupling_ok=bool(bx>0 and by>0)
    coupling=[bx/2,bx/2,bx/2,by]
    h=1/st["V"].sqrt(); gaps=[h,2*h]
    kappas=[1-2*x,1-2*y]; bm=[gaps[i]**2*kappas[i] for i in range(2)]
    g=[st["ga"],st["gb"]]; A=(3*st["D"]+st["z"]*st["phi"])/2
    J=3*st["D"]+(st["z"]**3+3*st["z"])*st["phi"]
    def off(i,j,gij):
        return gij+3*(bm[i]*g[j]+bm[j]*g[i])/2+J*bm[i]*bm[j]/4+(st["z"]**2+2)*st["phi"]*(bm[i]*gaps[j]+bm[j]*gaps[i])/2+st["z"]*st["phi"]*gaps[i]*gaps[j]
    def diag(i):
        return -2*A*gaps[i]**2+J*bm[i]**2/4+(st["z"]**2+2)*st["phi"]*bm[i]*gaps[i]+st["z"]*st["phi"]*gaps[i]**2+3*bm[i]*g[i]+6*st["C"]*gaps[i]**3*kappas[i]**2
    da,db=diag(0),diag(1)
    oo,ab=off(0,0,y-x),off(0,1,2*x-1)
    exchange=da-oo
    B11,B12,B22=3*da+6*oo,3*ab,db
    det=B11*B22-B12*B12
    hessian_ok=bool(exchange<0 and B11<0 and det>0)
    pmf=generic_pmf([x,x,x,y],[1,1,1,2])
    owner_gaps=[]
    owner_ok=True
    for t in range(6):
        z=(t-3*x-2*y)/st["V"].sqrt()
        Phi=(1+(z/arb(2).sqrt()).erf())/2
        plus=generic_cdf(pmf,t)-Phi
        minus=Phi-generic_cdf(pmf,t-1)
        for side,D in (("closed",plus),("strict",minus)):
            if t==2 and side=="closed":
                continue
            gap=st["D"]-D
            owner_gaps.append(dict(owner=t,side=side,gap=str(gap)))
            owner_ok=owner_ok and bool(gap>0)
    duals=[]
    dual_ok=True
    for group,p,a,probs,weights in ((0,x,1,[x,x,y],[1,1,2]),(1,y,2,[x,x,x],[1,1,1])):
        deleted=generic_pmf(probs,weights)
        s=a*a*p*(1-p); rho=a**3*p*(1-p)*(p*p+(1-p)**2)
        M=a*a*p*(1-p)*(1-2*p)
        q0=st["F"]+st["A"]*s-st["lam"]*rho
        B=-(st["w"]+3*st["lam"]*M)
        def gap(u,Fdel):
            return q0+B*u-st["A"]*u*u+st["lam"]*abs(u)**3-Fdel
        checks=[]
        for j in range(len(deleted)):
            if j in (2,2-a):
                continue # Contacts follow exactly from the certified G=0 root.
            u=2-a*p-j
            checks.append(("breakpoint",gap(u,generic_cdf(deleted,j))))
        extrema=[arb(0)]
        ambiguous=False
        for qa in (3*st["lam"],-3*st["lam"]):
            qb=-2*st["A"]; disc=qb*qb-4*qa*B
            if disc<0:
                continue
            if not disc>0:
                ambiguous=True; continue
            for u in ((-qb-disc.sqrt())/(2*qa),(-qb+disc.sqrt())/(2*qa)):
                if (qa>0 and u>0) or (qa<0 and u<0):
                    extrema.append(u)
                elif not (u>0 or u<0):
                    ambiguous=True
        for u in extrema:
            threshold=2-a*p-u
            j=math.floor(float(threshold))
            if threshold<0:
                Fdel=arb(0)
            elif threshold>len(deleted)-1:
                Fdel=generic_cdf(deleted,len(deleted)-1)
            elif threshold>j and threshold<j+1:
                Fdel=generic_cdf(deleted,j)
            else:
                ambiguous=True; continue
            checks.append(("stationary_or_zero",gap(u,Fdel)))
        ok=not ambiguous and all(v>0 for label,v in checks)
        dual_ok=dual_ok and ok
        duals.append(dict(group=group,certified=bool(ok),ambiguous=ambiguous,
            off_contact_gaps=[dict(type=label,gap=str(v)) for label,v in checks],
            exact_contact_basis="Probability stationarity G_i=0 plus zero mean score; original two support masses are positive."))
    return dict(feasibility_certified=feasible,coupling_certified=coupling_ok,coupling_balls=[str(v) for v in coupling],
        hessian_negative_definite_certified=hessian_ok,exchange_eigenvalue_ball=str(exchange),
        constant_group_B11_ball=str(B11),constant_group_determinant_ball=str(det),
        owner_unique_certified=owner_ok,owner_gap_balls=owner_gaps,
        coordinate_duals_certified=bool(dual_ok),coordinate_duals=duals,ratio_ball=str(st["C"]),
        nonglobal_comparison_certified=bool(st["C"]<arb(2).sqrt()/4),
        comparison="Two symmetric Bernoulli coordinates plus two zeros have ratio sqrt(2)/4 at the centered zero owner.")


def midpoint_preconditioner(center):
    """A fixed decimal inverse of the analytic Jacobian at the stored center."""
    x,y=map(arb,center)
    dual=arb_state(Dual(x,(1,0)),Dual(y,(0,1)))
    J=[[dual["G"][i].d[j] for j in range(2)] for i in range(2)]
    det=J[0][0]*J[1][1]-J[0][1]*J[1][0]
    if not (det>0 or det<0):
        raise ValueError("The center Jacobian is not certified nonsingular")
    inverse=[[J[1][1]/det,-J[0][1]/det],[-J[1][0]/det,J[0][0]/det]]
    return [[v.mid().str(90,radius=False) for v in row] for row in inverse]


def readable_root_check(raw):
    result=copy.deepcopy(raw)
    result.pop("contraction_bound_approx")
    # Outward-rounded zero-centered printed balls give rigorous row bounds.
    result["derivative_row_sum_upper_bounds"]=[
        str(sum((Decimal(v.removeprefix("[+/- ").removesuffix("]")) for v in row),Decimal(0)))
        for row in result["derivative_defect_balls"]
    ]
    return result


def readable_local_check(raw):
    result=copy.deepcopy(raw)
    result["attaining_threshold_unique_certified"]=result.pop("owner_unique_certified")
    result["threshold_discrepancy_margins"]=result.pop("owner_gap_balls")
    for item in result["threshold_discrepancy_margins"]:
        item["raw_threshold"]=item.pop("owner")
    result["comparison"]=result["comparison"].replace("centered zero owner","centered zero threshold")
    result.pop("all_local_systems_certified",None)
    return result


def verify(document):
    """Recheck this fixed box and the published concrete interval enclosures."""
    if document["schema"]!="local-trap-fixed-box-v1":
        raise ValueError("Unsupported certificate schema")
    inputs=document["inputs"]
    if (inputs["gaps"]!=[1,1,1,2] or inputs["selected_raw_threshold"]!=2
        or inputs["selected_side"]!="closed" or inputs["arb_precision_bits"]!=320
        or inputs["probability_box_radius"]!="1e-28"):
        raise ValueError("This verifier covers only the stated fixed four-variable box")
    previous_precision=ctx.prec
    try:
        ctx.prec=320
        center=inputs["probability_centers"]
        inverse=midpoint_preconditioner(center)
        if inverse!=inputs["preconditioner"]:
            raise ValueError("The fixed analytic-Jacobian preconditioner changed")
        box,root=krawczyk(center,inverse,inputs["probability_box_radius"])
        local=certificate(box)
        root_flags=("inclusion","preconditioner_nonsingular","contraction_certified")
        local_flags=("feasibility_certified","coupling_certified",
            "hessian_negative_definite_certified","owner_unique_certified",
            "coordinate_duals_certified","nonglobal_comparison_certified")
        if not all(root[k] for k in root_flags) or not all(local[k] for k in local_flags):
            raise ValueError("A required strict interval inequality failed")
        actual_root=readable_root_check(root)
        actual_local=readable_local_check(local)
        if actual_root!=document["stationary_root_check"]:
            raise ValueError("The root-check enclosures differ from the published certificate")
        if actual_local!=document["local_and_comparison_checks"]:
            raise ValueError("The local-check enclosures differ from the published certificate")
        state=arb_state(*box)
        return dict(verified=True,stationary_root_check=actual_root,
            local_and_comparison_checks=actual_local,
            feasibility_balls={k:str(state[k]) for k in ("x","y","V","R","D","lam")},
            comparison_ratio_ball=str(arb(2).sqrt()/4),
            comparison_advantage_ball=str(arb(2).sqrt()/4-state["C"]))
    finally:
        ctx.prec=previous_precision


def main():
    parser=argparse.ArgumentParser(description=(
        "Verify one saved four-Bernoulli probability box with gaps (1,1,1,2), "
        "radius 1e-28, and 320-bit Arb arithmetic. Rechecks a unique stationary "
        "root, the full probability Hessian, positive collision weights, all "
        "CDF thresholds, both coordinate inequalities on the real line, and "
        "the larger sqrt(2)/4 comparison. No optimization or root search."))
    parser.add_argument("--certificate",type=Path,
        default=Path(__file__).resolve().parents[2]/"certificates/observations/local-trap.json",
        help="Fixed-box certificate (default: certificates/observations/local-trap.json)")
    parser.add_argument("--json",action="store_true",help="Print every replayed enclosure")
    args=parser.parse_args()
    try:
        result=verify(json.loads(args.certificate.read_text()))
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(1,f"Verification failed: {exc}\n")
    if args.json:
        print(json.dumps(result,indent=2,sort_keys=True))
    else:
        checks=result["local_and_comparison_checks"]
        print("Verified the fixed unequal-gap probability box and all strict certificate inequalities.")
        print("Ratio:",checks["ratio_ball"])
        print("Larger four-slot comparison:",result["comparison_ratio_ball"])
        print("The binary-law local-stability argument is in docs/numerical-observations.md.")


if __name__=="__main__":
    main()
