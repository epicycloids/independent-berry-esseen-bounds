"""Check scalar certificates and reconstruct their coefficient tables.

The default check verifies partitions, exact correction and tail inequalities,
Bernstein bounds and logarithmic tangent bounds. It uses saved transcendental
leaf bounds. --replay FAMILY recomputes those bounds for the requested family,
serially. --limit N caps the number of complete affine certificates or energy
frequency records checked numerically; the result reports the actual coverage."""
import argparse
from fractions import Fraction as Q
import importlib
import json
from math import comb, factorial, isfinite, nextafter, inf
import os
from pathlib import Path
import sys
import zipfile

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

PUBLIC = Path(__file__).resolve().parent.parent
PRICES = tuple(map(Q, ("5/4", "3/2", "2", "3", "4")))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rational(value):
    require(not isinstance(value, bool), "A boolean is not a scalar coefficient")
    if isinstance(value, float):
        require(isfinite(value), "Nonfinite scalar coefficient")
    return Q(value)


def load(path):
    return json.loads(path.read_bytes())


def modules(root):
    directory = str(root/"certificates/scalar/source")
    if directory not in sys.path:
        sys.path.insert(0, directory)


def check_polygon(polygon):
    eta = Q(997,1000)
    upper = [tuple(map(Q,p)) for p in polygon["upper"]]
    full = [tuple(map(Q,p)) for p in polygon["full"]]
    require(len(upper)==25 and upper[0]==(1,0) and upper[12]==(0,1) and upper[-1]==(-1,0),"Directional polygon axes differ")
    require(len(set(upper))==25 and all(y>0 for _,y in upper[1:-1]),"Bad upper directions")
    require(full==upper+[(x,-y) for x,y in upper[-2:0:-1]] and len(set(full))==48,"Directional polygon reflection differs")
    def cross(a,b): return a[0]*b[1]-a[1]*b[0]
    for i,a in enumerate(full):
        b=full[(i+1)%48]
        edge=(b[0]-a[0],b[1]-a[1])
        determinant=cross(a,b)
        require(determinant>0 and determinant**2>=eta**2*(edge[0]**2+edge[1]**2),"Polygon edge misses required disk")
        require(all(cross(edge,(v[0]-a[0],v[1]-a[1]))>0 for j,v in enumerate(full) if j not in (i,(i+1)%48)),"Polygon supporting half-plane fails")
    require(Q(polygon["inradius_lower"])==eta,"Wrong polygon inradius")
    return upper


def tree_leaves(c, disk=False):
    """Reconstruct each closed rectangle from its exact binary split tree."""
    tree=c["tree"]
    xl,xr=map(Q,c["x_interval"])
    radius=Q(c["tail_radius"])
    require(0<=xl<xr<=2 and radius>=1 and tree,"Invalid scalar root domain")
    pending=[(0,(xl,xr,Q(0) if disk else -radius,radius))]
    seen=set()
    leaves=0
    while pending:
        i,box=pending.pop()
        require(type(i) is int and 0<=i<len(tree) and i not in seen,"Scalar partition is not a finite tree")
        seen.add(i)
        node=tree[i]
        if len(node)==(2 if disk else 1):
            leaves+=1
            yield box,tuple(map(Q,node))
            continue
        require(len(node)==4 and type(node[0]) is int and node[0] in (0,1),"Unknown scalar node format")
        axis,cut,left,right=node
        cut=Q(cut)
        a,b,l,r=box
        if axis==0:
            require(a<cut<b,"Frequency split does not partition its parent")
            boxes=((a,cut,l,r),(cut,b,l,r))
        else:
            require(l<cut<r,"Spatial split does not partition its parent")
            boxes=((a,b,l,cut),(a,b,cut,r))
        pending.extend(zip((left,right),boxes))
    require(len(seen)==len(tree) and leaves==c["leaves"],"Unreachable scalar nodes or wrong leaf count")


def disk_tail(c):
    beta,gamma,c0,c2,lam=map(Q,c["coefficients"])
    delta=Q(c["correction"])
    xl,xr=map(Q,c["x_interval"])
    radius=Q(c["tail_radius"])
    if lam>1:
        a3,a2=lam-1,c2-abs(gamma)
        a1=-(xr*xr/2+abs(beta))
        a0=c0+delta+abs(gamma)
    else:
        require(xl>0 and lam>0,"Invalid small-price centered tail")
        bound=max(abs(gamma-2/xl),abs(gamma-2/xr))
        a3,a2=lam,c2-bound
        a1=-(4/(xl*xl)+abs(beta+2))
        a0=c0+delta+bound-2/xl
    coefficients=(a0+a1*radius+a2*radius**2+a3*radius**3,
        a1+2*a2*radius+3*a3*radius**2,a2+3*a3*radius,a3)
    require(coefficients==tuple(map(Q,c["tail_coefficients"])) and min(coefficients)>=0,"Disk analytic tail algebra fails")


def check_affine(root, index, family, replay=False, limit=None):
    from flint import ctx
    disk=family=="disk"
    module=importlib.import_module("disk" if disk else "directional")
    data=index[family]
    require(tuple(map(Q,data["prices"]))==PRICES,"Affine prices differ")
    upper=None if disk else check_polygon(data["polygon"])
    expected_cells=([(Q(i,64),Q(i+1,64)) for i in range(128)] if disk else
        [(Q(0),Q(1,8))]+[(Q(i,32),Q(i+1,32)) for i in range(4,64)])
    require(len(data["rows"])==len(expected_cells),"Wrong scalar frequency cell count")
    constants=[]
    cumulative=[None]*5
    certificates=leaves=nodes=replayed_certificates=replayed_leaves=0
    expected_entries=set()
    with zipfile.ZipFile(root/"certificates/scalar"/(family+".zip")) as archive:
        for row,interval in zip(data["rows"],expected_cells):
            require(tuple(map(Q,row["x"]))==interval and len(row["certificates"])==5,"Affine frequency coverage differs")
            row_constants=[]
            for price,names in zip(PRICES,row["certificates"]):
                require(len(names)==(1 if disk else 25),"Incomplete scalar directions")
                direction_constants=[]
                for direction_index,name in enumerate(names):
                    require(name not in expected_entries,"Repeated scalar certificate")
                    expected_entries.add(name)
                    c=json.loads(archive.read(name))
                    require(tuple(map(Q,c["x_interval"]))==interval,"Scalar payload frequency differs")
                    correction=Q(c["correction"])
                    require(correction>=0 and type(c["precision_bits"]) is int and c["precision_bits"]>=64,"Invalid scalar correction or precision")
                    actual_replay=replay and (limit is None or replayed_certificates<limit)
                    with ctx.workprec(c["precision_bits"]):
                        if disk:
                            beta,gamma,c0,c2,source_price=map(Q,c["coefficients"])
                            require(0<source_price<=price,"Ineligible direct-disk price")
                            base=Q(c["base_correction"])
                            require(0<=base<=correction,"Invalid disk base correction")
                            extra=correction-base
                            disk_tail(c)
                            if actual_replay:
                                problem=module.DiskProblem.make(c["coefficients"],c["x_interval"],center_moments=True)
                                evaluator=module._Evaluator(problem,base)
                        else:
                            c0,c1,c2,source_price=map(Q,c["coefficients"])
                            require(source_price==price and tuple(map(Q,c["direction"]))==upper[direction_index],"Directional identity differs")
                            problem=module.Problem.make(c["coefficients"],c["x_interval"],c["direction"],center_moments=True)
                            evaluator=module._Evaluator(problem)
                            # This is a short scalar Arb tail calculation, not
                            # a replay of any compact transcendental leaf.
                            tail=evaluator.tail(Q(c["tail_radius"]))
                            stored=c["tail_coefficients"]
                            require(len(stored)==2 and all(len(r)==4 for r in stored),"Incomplete directional tail")
                            require(all(v>=0 for r in tail for v in r),"Directional analytic tail fails")
                            require(all(Q(v)>=0 for r in stored for v in r),"Negative saved directional tail")
                        for rectangle,saved in tree_leaves(c,disk):
                            if disk:
                                low_p,low_h=saved
                                require(low_p+extra>=0 and low_h+2*extra*low_p+extra**2>=0,"Saved disk leaves do not imply corrected positivity")
                                if actual_replay:
                                    fresh=module._cell(evaluator,*rectangle)
                                    fp,fh=fresh["lower_p"],fresh["lower_h"]
                                    require(fp+extra>=0 and fh+2*extra*fp+extra**2>=0,"Disk leaf numerical replay failed")
                            else:
                                require(saved[0]+correction>=0,"Saved directional leaf exceeds its correction")
                                if actual_replay:
                                    require(module._node(evaluator,*rectangle).lower+correction>=0,"Directional leaf numerical replay failed")
                            leaves+=1
                            replayed_leaves+=bool(actual_replay)
                    certificates+=1
                    nodes+=len(c["tree"])
                    replayed_certificates+=bool(actual_replay)
                    intercept=c0+c2+correction+source_price-price
                    direction_constants.append(intercept)
                constant=max(direction_constants)
                require(constant+price>=0,"Scalar norm bound is negative at rho=1")
                row_constants.append(constant)
            if disk:
                cumulative=[value if old is None else max(old,value) for old,value in zip(cumulative,row_constants)]
                constants.append((row_constants,list(cumulative)))
            else:
                constants.append(row_constants)
        require(set(archive.namelist())==expected_entries,"Missing or unused affine payloads")
    require((certificates,nodes,leaves)==tuple(data[k] for k in ("certificate_count","node_count","leaf_count")),"Affine corpus counts differ")
    runtime=load(root/"verify/kernel"/("disk_table.json" if disk else "scalar_table.json"))
    require(tuple(map(Q,runtime["prices"]))==PRICES and len(runtime["rows"])==len(constants),"Runtime affine table shape differs")
    for row,interval,compiled in zip(runtime["rows"],expected_cells,constants):
        require(tuple(map(Q,(row["xlo"],row["xhi"])))==interval,"Runtime affine interval differs")
        if disk:
            require(list(map(Q,row["cell_constants"]))==compiled[0] and list(map(Q,row["constants"]))==compiled[1],"Runtime disk prefix is not compiled from the supplied certificates")
        else:
            require(list(map(Q,row["constants"]))==compiled,"Runtime directional table is not compiled from supplied certificates")
    if not disk:
        require(runtime["direction_polygon"]==data["polygon"] and Q(runtime["inradius_lower"])==Q(997,1000),"Runtime direction geometry differs")
    return {"certificates":certificates,"nodes":nodes,"leaves":leaves,
        "all_partitions_and_table_compilation_checked":True,
        "numerically_replayed_certificates":replayed_certificates,
        "numerically_replayed_leaves":replayed_leaves,
        "complete_numerical_leaf_replay":replayed_certificates==certificates}


def check_odd(root):
    import odd_polynomial as odd
    data=load(root/"certificates/scalar/odd-quadratic.json")
    require(data["order"]==24 and list(map(Q,data["rectangle"]))==[0,4,0,49],"Wrong odd certificate domain/order")
    polynomials=odd.components(24)
    for ij,c in [((0,1),-Q(1,6)),((1,1),-Q(1,36)),((1,0),-Q(1,24)),((2,0),-Q(1,144))]: odd.add(polynomials[0],ij,c)
    for ij,c in [((1,0),-Q(1,6)),((2,0),-Q(1,36))]: odd.add(polynomials[2],ij,c)
    B25=Q(49,1764)*25**2*196**25/factorial(48)
    tail=Q(5,18)*B25
    require(Q(196*26**2,25**2*50*49)<Q(1,10),"Odd analytic series-tail ratio fails")
    minima={}
    for name,p in zip(("F1","F0","G"),polynomials):
        stored=data["components"][name]
        require([[a,b,str(c)] for (a,b),c in sorted(p.items())]==stored["powers"],"Odd exact polynomial coefficients differ")
        minimum=odd.bernstein_lower(p,tuple(map(Q,data["rectangle"])))
        require(minimum==Q(stored["bernstein_minimum"]) and tail==Q(stored["absolute_tail_upper"]),"Odd Bernstein minimum or tail differs")
        require(minimum-tail>Q(1,250),"Odd polynomial positivity fails")
        minima[name]=str(minimum)
    require(Q(data["component_lower"])==Q(1,250) and Q(data["spatial_tail_ratio_upper"])==Q(16,25)==Q(400,5**4)<1,"Odd spatial tail fails")
    return {"order":24,"exact_polynomials_regenerated":True,"bernstein_minima":minima,"analytic_series_and_spatial_tails_checked":True}


def check_direct_log(root):
    from flint import arb,ctx
    data=load(root/"certificates/scalar/direct-log.json")
    def a(q):
        q=Q(q)
        return arb(q.numerator)/q.denominator
    active=list(map(Q,("3/4","1","3/2","2","3")))
    require([Q(r["price"]) for r in data["supports"]]==active,"Wrong active full-log prices")
    require(Q(data["kappa"])==Q(49581,500000),"Wrong cubic envelope constant")
    count=0
    with ctx.workprec(192):
        K=a(data["kappa"])
        A=1/(32*K*K)
        require(0<A<arb(27)/8,"Full-log shape condition fails")
        for group in data["supports"]:
            lam=a(group["price"])
            qstar=lam/(2*(1+lam))
            conjugate=(lam-(1+lam).log())/2
            require(qstar<=a(group["qstar_upper"]) and a(group["conjugate_lower"])<=conjugate,"Full-log cap constants not outward")
            previous_u=Q(-1)
            previous_domain=Q(1)
            for j,row in enumerate(group["rows"]):
                u,alo,bhi,M=map(Q,row)
                require(previous_u<u and u<=M<=previous_domain<=1,"Full-log lookup order or domain fails")
                previous_u,previous_domain=u,M
                if j==0:
                    require((u,alo,bhi,M)==(0,0,0,1),"Missing full-log zero support")
                    count+=1
                    continue
                z=a(u)
                r=A*z*z*(1-z)
                g=-(1-2*r).log()/2-r
                gp=2*r/(1-2*r)
                gpp=2/(1-2*r)**2
                zrp=A*z*z*(2-3*z)
                D=gpp*zrp*zrp-6*r*gp+6*g
                alpha=(gp*zrp-2*g)/z**3
                beta=alpha*z-g/z**2
                require(z<arb(2)/3 and r<=qstar and D>=0 and alpha>=0 and beta>=0,"Full-log tangent shape conditions fail")
                endpoint=a(M)
                re=A*endpoint**2*(1-endpoint)
                phi=-(1-2*re).log()/2-re if re<=qstar else lam*re-conjugate
                require(phi/endpoint**2>=alpha*endpoint-beta,"Full-log endpoint does not certify the support")
                require(0<=alo and a(alo)<=64*K**3*alpha and 0<=bhi and 16*K*K*beta<=a(bhi),"Full-log support coefficients not outward")
                count+=1
    return {"active_prices":[str(p) for p in active],"support_rows":count,"all_tangents_and_endpoints_checked_with_arb":True,"precision_bits":192}


def rectangle_cover(record):
    xl,xr=map(Q,record["x"])
    ymax=Q(record["y_max"])
    boxes=[]
    for row in record["rectangles"]:
        require(len(row)==5,"Wrong energy rectangle shape")
        a,b,l,r,lower=map(Q,row)
        require(xl<=a<=b<=xr and 0<=l<r<=ymax and lower>=0,"Invalid energy rectangle or saved lower bound")
        boxes.append((a,b,l,r,lower))
    require(boxes,"Empty energy certificate")
    xnodes=sorted({xl,xr,*(v for box in boxes for v in box[:2])})
    # Check closed boundaries and a representative of every open slab;
    # the set of covering closed rectangles is constant inside each slab.
    for x in xnodes+[(a+b)/2 for a,b in zip(xnodes[:-1],xnodes[1:])]:
        end=Q(0)
        for l,r in sorted((l,r) for a,b,l,r,_ in boxes if a<=x<=b):
            require(l<=end,"Energy scalar rectangle cover has a gap")
            end=max(end,r)
        require(end==ymax,"Energy scalar rectangle cover misses its endpoint")
    return boxes


def check_energy(root,replay=False,limit=None):
    from flint import arb,ctx
    import energy_majorant
    import energy_small
    data=load(root/"certificates/scalar/energy.json")
    rectangles=replayed_rectangles=replayed_records=0
    counts={}
    with ctx.workprec(128):
        for family,module,start in (("small",energy_small,0),("upper",energy_majorant,64)):
            records=data[family]
            require(len(records)==64,"Incomplete energy frequency family")
            coefficients=module.tail_coefficients()
            require(all(c>0 for c in coefficients),"Energy uniform-tail polynomial fails")
            count=0
            for i,record in enumerate(records,start):
                require(tuple(map(Q,record["x"]))==(Q(i,128),Q(i+1,128)) and Q(record["price"])==Q(9,8),"Energy frequency cover or price differs")
                require(record["y_max"]==(48 if family=="small" else 16),"Wrong energy compact domain")
                require(tuple(map(Q,record["tail_coefficients"]))==coefficients,"Energy exact tail coefficients differ")
                boxes=rectangle_cover(record)
                actual_replay=replay and (limit is None or replayed_records<limit)
                if actual_replay:
                    function=module.normalized_gap if family=="small" else module.gap_quotient
                    for a,b,l,r,lower in boxes:
                        value=function(module.interval(a,b),module.interval(l,r),Q(9,8))
                        require(module.exact(lower)<=value and value>=0,"Energy numerical scalar rectangle replay failed")
                    replayed_records+=1
                    replayed_rectangles+=len(boxes)
                rectangles+=len(boxes)
                count+=len(boxes)
            counts[family]=count
        # The zero-frequency identity is exact, including its nonzero
        # normalized limit. Compare polynomial coefficients, not samples.
        # Insert E_n(0)=1/n! into the supplied entire formula with z=y-1.
        # The resulting coefficients at lambda=9/8 are shown explicitly.
        price=Q(9,8)
        # Substitute c=S=u=a=1, R=14/9, K=-1 and E_n(0)=1/n!
        # into the entire formula; each row is a polynomial contribution.
        terms=((price/16,price/4,price/4),(Q(35,18),Q(7,9),Q(0)),
            (Q(-1),Q(0),Q(0)),(Q(-1),Q(-1),Q(0)),
            (Q(0),Q(0),Q(-1,3)),(Q(0),Q(0),Q(1,12)))
        zero=tuple(sum((row[j] for row in terms),Q(0)) for j in range(3))
        expected=(Q(1,144)+(price-1)/16,Q(1,36)+(price-1)/4,(price-1)/4)
        require(zero==expected,"The entire energy formula has the wrong zero-frequency limit")
        require(all(c>0 for c in zero),"Zero-frequency energy gap fails")
        runtime=load(root/"verify/kernel/energy_table_complete.json")
        rows=runtime["rows"]
        require(len(rows)==65 and runtime["new_upper_interval_complete"] is True,"Incomplete runtime energy table")
        q0=energy_small.q_upper(Q(1,2))
        price0=float(Q(9,8))
        first={"left":"0","right":"1/2","branch":"small-energy","Q_upper":q0,
            "lambda_upper":price0,"prefix_Q_upper":q0,"prefix_lambda_upper":price0}
        require(rows[0]==first,"Runtime small energy prefix differs")
        prefix_q,prefix_price=q0,price0
        for i,row in enumerate(rows[1:],64):
            right=Q(i+1,128)
            x=energy_majorant.exact(right)
            q=energy_majorant.upper(arb(8)*x.sin()**4/(3*x**4*x.cos()**2))
            price=energy_majorant.upper(energy_majorant.exact(Q(9,8)))
            prefix_q,prefix_price=max(prefix_q,q),max(prefix_price,price)
            expected={"left":str(Q(i,128)),"right":str(right),"branch":"new","Q_upper":q,
                "lambda_upper":price,"prefix_Q_upper":prefix_q,"prefix_lambda_upper":prefix_price}
            require(row==expected,"Runtime energy prefix is not compiled from supplied scalar cells")
    return {"scalar_records":128,"rectangles":rectangles,"rectangles_by_family":counts,
        "all_rectangle_covers_and_table_compilation_checked":True,
        "zero_frequency_gap_coefficients":[str(c) for c in zero],
        "numerically_replayed_records":replayed_records,"numerically_replayed_rectangles":replayed_rectangles,
        "complete_numerical_leaf_replay":replayed_records==128}


def check_all(root=PUBLIC, *, replay=(), limit=None):
    root=Path(root)
    modules(root)
    index=load(root/"certificates/scalar/affine-index.json")
    result={"schema":"current-scalar-certificate-check-v1",
        "default_scope":"All exact saved topology, rational leaf consequences, analytic tails, scalar table compilation; odd Bernstein regenerated and logarithmic tangents certified. Saved transcendental rectangle enclosures remain numerical inputs unless explicitly replayed.",
        "odd_quadratic":check_odd(root),"direct_log":check_direct_log(root)}
    for family in ("directional","disk"):
        result[family]=check_affine(root,index,family,family in replay,limit)
    result["energy"]=check_energy(root,"energy" in replay,limit)
    result["complete_numerical_leaf_replay"]=all(result[k]["complete_numerical_leaf_replay"] for k in ("directional","disk","energy"))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay",choices=("directional","disk","energy"),action="append",default=[])
    parser.add_argument("--limit",type=int,help="maximum affine certificates or energy frequency records to recompute per family")
    parser.add_argument("--root",type=Path,default=PUBLIC)
    args=parser.parse_args()
    require(args.limit is None or args.limit>0,"Numerical reevaluation limit must be positive")
    print(json.dumps(check_all(args.root,replay=args.replay,limit=args.limit),indent=2,sort_keys=True))


if __name__=="__main__":
    main()
