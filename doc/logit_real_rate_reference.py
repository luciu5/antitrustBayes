"""Independent quantity-profit check for real input Cournot rates.

Pure standard library; no package solver, Stan or posterior refit.
The outside-price shift is absorbed in delta in this reference.
"""
import json
import math


def w_from_log(loga):
    lo = min(-50.0, loga - 1)
    hi = max(1.0, math.log(max(1.0, loga)) + 1)
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid + math.exp(mid) < loga:
            lo = mid
        else:
            hi = mid
    return math.exp((lo + hi) / 2)


def equilibrium(delta, values, beta, owner, roles, stack):
    firms = list(dict.fromkeys(owner))
    aj = [math.exp(d + beta * v) for d, v in zip(delta, values)]
    af = {f: sum(a for a, g in zip(aj, owner) if g == f) for f in firms}
    x = {}
    if stack:
        for f in firms:
            if f not in roles:
                x[f] = w_from_log(math.log(af[f]) - 1)
        reserve = 1 + sum(x.values())
        for f in roles:
            x[f] = reserve * w_from_log(math.log(af[f]) - 1 - math.log(reserve))
        h = {f: 1 + x[f] / reserve if f in roles else 1 + x[f] for f in firms}
    else:
        reserve = None
        x = {f: w_from_log(math.log(af[f]) - 1) if f in roles
             else af[f] / math.e for f in firms}
        h = {f: 1 + x[f] if f in roles else 1 for f in firms}
    s0 = 1 / (1 + sum(x.values()))
    p = [v - h[f] / beta for v, f in zip(values, owner)]
    q = [s0 * a * math.exp(-h[f]) for a, f in zip(aj, owner)]
    return p, q, x, h, aj, reserve


def profit(q, delta, values, beta, owner, firm):
    outside = 1 - sum(q)
    rates = [(math.log(s / outside) - d) / beta for s, d in zip(q, delta)]
    return sum((v - p) * s for v, p, s, f in zip(values, rates, q, owner)
               if f == firm)


def maximize(objective, upper):
    lo, hi = 1e-12, upper - 1e-12
    ratio = (math.sqrt(5) - 1) / 2
    a, b = hi - ratio * (hi - lo), lo + ratio * (hi - lo)
    fa, fb = objective(a), objective(b)
    for _ in range(150):
        if fa < fb:
            lo, a, fa = a, b, fb
            b = lo + ratio * (hi - lo)
            fb = objective(b)
        else:
            hi, b, fb = b, a, fa
            a = hi - ratio * (hi - lo)
            fa = objective(a)
    return (lo + hi) / 2


rows = []
beta = .4
delta = [.3, -.2, .1, -.4, -.1]
base_values = [.4, .5, .2, .8, .3]
for owner in (['a', 'b', 'c', 'd', 'e'], ['a', 'a', 'b', 'c', 'c']):
    firms = list(dict.fromkeys(owner))
    for stack in (False, True):
        for roles in (firms[:1], firms[:2], firms):
            for shock in (0, .1):
                values = [v * (1 + shock) for v in base_values]
                p, q, x, h, aj, reserve = equilibrium(
                    delta, values, beta, owner, roles, stack)
                assert min(p) < 0
                utilities = [d + beta * z for d, z in zip(delta, p)]
                weights = [math.exp(u) for u in utilities]
                demand = [w / (1 + sum(weights)) for w in weights]
                share_error = max(abs(a - b) for a, b in zip(q, demand))
                assert share_error < 2e-14
                for j, f in enumerate(owner):
                    if not stack and f not in roles:
                        continue  # MonCom's perceived derivative differs from Cournot.
                    leader = stack and f in roles
                    fixed = sum(z for i, (z, g) in enumerate(zip(q, owner))
                                if i != j and (not leader or g in roles))
                    upper = 1 - fixed

                    def objective(choice):
                        candidate = q.copy()
                        candidate[j] = choice
                        if leader:
                            remaining = 1 - sum(z for z, g in zip(candidate, owner)
                                                if g in roles)
                            for i, g in enumerate(owner):
                                if g not in roles:
                                    candidate[i] = remaining / reserve * aj[i] * math.exp(-h[g])
                        return profit(candidate, delta, values, beta, owner, f)

                    eps = 1e-6
                    gradient = (objective(q[j] + eps) - objective(q[j] - eps)) / (2 * eps)
                    best_response = maximize(objective, upper)
                    optimum_error = abs(best_response - q[j])
                    assert abs(gradient) < 2e-8
                    assert optimum_error < 2e-8
                    rows.append(dict(model='CS' if stack else 'CF', multiproduct=len(firms) < 5,
                                     leaders_or_core=len(roles), shock=shock, product=j,
                                     min_rate=min(p), share_error=share_error,
                                     gradient_error=abs(gradient), optimum_error=optimum_error))
print(json.dumps(dict(checks=len(rows),
                     max_share_error=max(r['share_error'] for r in rows),
                     max_gradient_error=max(r['gradient_error'] for r in rows),
                     max_optimum_error=max(r['optimum_error'] for r in rows),
                     results=rows), indent=2))
