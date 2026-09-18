import math, json

def log_w(loga):
    # w*exp(w)=exp(loga), solved without exp(loga).
    lo, hi = -max(50.0, abs(loga)+50.0), max(1.0, loga+1.0)
    for _ in range(200):
        t = (lo+hi)/2
        if t+math.exp(t) < loga: lo=t
        else: hi=t
    return math.exp((lo+hi)/2)

def state(delta, costs, beta, owner, leaders=None, core=None):
    a=abs(beta)
    groups=list(dict.fromkeys(owner))
    aa=[math.exp(d+beta*c) for d,c in zip(delta,costs)]
    af={g:sum(v for v,f in zip(aa,owner) if f==g) for g in groups}
    xx={}
    if leaders is not None:
        for g in groups:
            if g not in leaders: xx[g]=log_w(math.log(af[g])-1)
        reserve=1+sum(xx.values())
        for g in leaders: xx[g]=reserve*log_w(math.log(af[g])-1-math.log(reserve))
        hh={g:1+xx[g]/reserve if g in leaders else 1+xx[g] for g in groups}
    else:
        xx={g:log_w(math.log(af[g])-1) if g in core else af[g]/math.e for g in groups}
        hh={g:1+xx[g] if g in core else 1 for g in groups}
    p=[c-hh[f]/beta for c,f in zip(costs,owner)]
    q0=1/(1+sum(xx.values()))
    q=[q0*v*math.exp(-hh[f]) for v,f in zip(aa,owner)]
    return p,q,q0,xx,hh,aa

def inverse(q,delta,beta):
    q0=1-sum(q)
    return [(math.log(z/q0)-d)/beta for z,d in zip(q,delta)]

def profit(q, delta, costs, beta, owner, firm):
    p=inverse(q,delta,beta)
    sign=-1 if beta>0 else 1
    return sign*sum((z-c)*s for z,c,s,f in zip(p,costs,q,owner) if f==firm)

checks=[]
for beta in (-6.0,6.0):
    owner=['a','a','b','c','d','d']
    costs=[1.2,1.4,1.0,1.1,1.3,1.25]
    delta=[-beta*c+d for c,d in zip(costs,(.1,-.4,.2,-.1,-.5,-.7))]
    for model, roles in [('CF',{'a','b'}),('C',set(owner)),('CS',{'a','b'})]:
        p,q,q0,xx,hh,aa=state(delta,costs,beta,owner,
          leaders=roles if model=='CS' else None, core=roles if model!='CS' else None)
        eu=[math.exp(d+beta*z) for d,z in zip(delta,p)]
        demand=[z/(1+sum(eu)) for z in eu]
        error=max(abs(s-t) for s,t in zip(q,demand))
        assert error<2e-14
        assert abs(q0-(1-sum(q)))<2e-14
        for j, firm in enumerate(owner):
            if model!='CS' and firm not in roles: continue
            def objective(step):
                qp=q.copy(); qp[j]+=step
                if model=='CS' and firm in roles:
                    followers={f for f in owner if f not in roles}
                    rem=1-sum(z for z,f in zip(qp,owner) if f in roles)
                    reserve=1+sum(xx[f] for f in followers)
                    for i,f in enumerate(owner):
                        if f in followers: qp[i]=rem/reserve*aa[i]*math.exp(-hh[f])
                return profit(qp,delta,costs,beta,owner,firm)
            eps=1e-6
            grad=(objective(eps)-objective(-eps))/(2*eps)
            assert abs(grad)<2e-8, (model,beta,j,grad)
            checks.append({'model':model,'beta':beta,'product':j,'gradient_error':abs(grad),'share_error':error})
print(json.dumps({'checks':len(checks),'max_gradient_error':max(x['gradient_error'] for x in checks),'max_share_error':max(x['share_error'] for x in checks),'cases':checks},indent=2))
