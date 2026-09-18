#!/usr/bin/env python3
"""Small-fixture oracle only: Decimal dense Cholesky, no external packages."""
import csv, json, sys
from decimal import Decimal, localcontext
from pathlib import Path

PI = '3.14159265358979323846264338327950288419716939937510582097494459230781640628620899'

def dense_log_density(fixture, tau):
    D = Decimal
    rs = list(map(D, fixture['share_residual']))
    rm = list(map(D, fixture['margin_residual']))
    z = rs + rm
    n, q = len(rs), len(rm)
    market, group = fixture['market'], fixture['group']
    idx = [i-1 for i in fixture['margin_idx']]
    ss, sm, rho = map(D, (fixture['sigma_share'], fixture['sigma_margin'], fixture['rho']))
    groups = {m:len({g for g,mm in zip(group,market) if mm==m}) for m in set(market)}
    covariance = [[D(0) for _ in z] for _ in z]
    for i in range(n):
        for j in range(n):
            if market[i] == market[j]:
                covariance[i][j] = tau*tau*(D(int(group[i]==group[j]))-D(1)/groups[market[i]])
        covariance[i][i] += ss*ss
    for j in range(q):
        covariance[n+j][n+j] = sm*sm
        covariance[idx[j]][n+j] = covariance[n+j][idx[j]] = rho*ss*sm
    L = [[D(0) for _ in z] for _ in z]
    for i in range(len(z)):
        for j in range(i+1):
            value = covariance[i][j]-sum((L[i][k]*L[j][k] for k in range(j)),D(0))
            L[i][j] = value.sqrt() if i==j else value/L[j][j]
    solved = []
    for i in range(len(z)):
        solved.append((z[i]-sum((L[i][j]*solved[j] for j in range(i)),D(0)))/L[i][i])
    return -(D(len(z))*(D(2)*D(PI)).ln()+D(2)*sum((L[i][i].ln() for i in range(len(z))),D(0))+
             sum((x*x for x in solved),D(0)))/D(2)

def main(path):
    path = Path(path)
    f = json.loads(path.read_text())
    rows = []
    for precision in (50,80):
        with localcontext() as context:
            context.prec = precision
            tau = Decimal(f['tau'])
            value = dense_log_density(f,tau)
            native = Decimal(f['native_likelihood'])
            rows.append(dict(check='dense_density_'+str(precision)+'_digits',
                reference=str(value), native=str(native), error=str(abs(value-native)),
                tolerance='0.00002',passed=abs(value-native)<=Decimal('.00002')))
            def target(log_shift):
                t = tau*log_shift.exp()
                # All other proper prior constants cancel in this derivative.
                return dense_log_density(f,t)-((t-1)/Decimal(f['prior_tau_sd']))**2/2
            for h in ('0.00002','0.000004'):
                step=Decimal(h)
                fd=(target(step)-target(-step))/(2*step)
                ad=Decimal(f['native_tau_AD'])
                error=abs(fd-ad)/(1+abs(ad))
                rows.append(dict(check='tau_AD_FD_'+str(precision)+'_digits_h'+h,
                    reference=str(fd), native=str(ad), error=str(error),
                    tolerance='0.00002',passed=error<=Decimal('.00002')))
    output=path.with_name('high-precision-validation.csv')
    with output.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for row in rows:print(row)
    if not all(row['passed'] for row in rows):sys.exit(1)

if __name__ == '__main__':main(sys.argv[1])
