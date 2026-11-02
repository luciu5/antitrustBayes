#!/usr/bin/env python3
"""Independent reference checks for the firm-uniform tariff route.

This file deliberately has no dependency on antitrustBayes, coordination, R,
or scipy.  It computes demand, inverse demand, physical-cost tariff profits,
effective-cost profits, finite-difference FOCs, and one-firm best responses
directly.  It is a small mathematical oracle, rather than a production game
solver.

The output CSV is evidence for the restricted route in
doc/tariff_posterior_reuse_plan.md: every post-policy firm must have one
positive retention.  BF/CF rows use a two-firm strategic core with one fringe
follower; BS/CS rows use one leader and two followers.  Responses are
conditional one-firm responses, which keeps this reference independent from
native coordination equilibrium code while testing the payoff and follower
identities used by that code.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Mapping, Sequence, Tuple


EPS = 1.0e-10


def _exp(x: float) -> float:
    # Fixture values stay moderate, but this protects the independent demand
    # implementation from accidental overflow when it is reused.
    if x > 700.0:
        return math.exp(700.0)
    if x < -745.0:
        return 0.0
    return math.exp(x)


def _max_abs(values: Iterable[float]) -> float:
    vals = [abs(float(v)) for v in values]
    return max(vals) if vals else 0.0


def _max_diff(x: Sequence[float], y: Sequence[float]) -> float:
    return _max_abs(a - b for a, b in zip(x, y))


@dataclass(frozen=True)
class Demand:
    """A generic output-market Logit or CES demand system."""

    kind: str
    delta_or_quality: Tuple[float, ...]
    slope: float
    scale: float
    outside_price: float = 1.0

    def quantity(self, prices: Sequence[float]) -> List[float]:
        if self.kind == "logit":
            alpha = self.slope
            z = [d - alpha * p for d, p in zip(self.delta_or_quality, prices)]
            zmax = max([0.0] + z)
            ez = [_exp(v - zmax) for v in z]
            den = _exp(-zmax) + sum(ez)
            return [self.scale * v / den for v in ez]
        if self.kind == "ces":
            gamma = self.slope
            p0 = self.outside_price
            quality = self.delta_or_quality
            den = p0 ** (1.0 - gamma) + sum(
                a * p ** (1.0 - gamma) for a, p in zip(quality, prices)
            )
            return [
                self.scale * a * p ** (-gamma) / den
                for a, p in zip(quality, prices)
            ]
        raise ValueError(f"unknown demand kind {self.kind!r}")

    def inverse(self, quantities: Sequence[float]) -> List[float]:
        if self.kind == "logit":
            total = sum(quantities)
            outside = self.scale - total
            if outside <= 0.0 or any(q <= 0.0 for q in quantities):
                raise ValueError("Logit inverse requires q_j > 0 and sum(q) < M")
            alpha = self.slope
            return [
                (d - math.log(q / outside)) / alpha
                for d, q in zip(self.delta_or_quality, quantities)
            ]
        if self.kind == "ces":
            gamma = self.slope
            if gamma <= 1.0 or any(q <= 0.0 for q in quantities):
                raise ValueError("CES inverse requires gamma > 1 and q_j > 0")
            quality = self.delta_or_quality
            A = self.outside_price ** (1.0 - gamma)
            H = sum(
                a ** (1.0 / gamma) * q ** ((gamma - 1.0) / gamma)
                for a, q in zip(quality, quantities)
            )

            # With z > 0, Marshallian CES inversion reduces to
            # B = A z^gamma + H z and p_j = z (a_j/q_j)^(1/gamma).
            def equation(z: float) -> float:
                return A * z ** gamma + H * z - self.scale

            lo, hi = 0.0, 1.0
            while equation(hi) < 0.0:
                hi *= 2.0
            for _ in range(180):
                mid = 0.5 * (lo + hi)
                if equation(mid) < 0.0:
                    lo = mid
                else:
                    hi = mid
            z = 0.5 * (lo + hi)
            return [z * (a / q) ** (1.0 / gamma) for a, q in zip(quality, quantities)]
        raise ValueError(f"unknown demand kind {self.kind!r}")


def finite_gradient(fn: Callable[[Sequence[float]], float], x: Sequence[float]) -> List[float]:
    """Central finite differences used by the independent FOC checks."""
    out: List[float] = []
    for j, value in enumerate(x):
        h = max(2.0e-6, abs(float(value)) * 2.0e-5)
        xp, xm = list(x), list(x)
        xp[j] += h
        xm[j] -= h
        out.append((fn(xp) - fn(xm)) / (2.0 * h))
    return out


def _firm_products(owner: Sequence[str], firm: str) -> List[int]:
    return [j for j, f in enumerate(owner) if f == firm]


def _firms(owner: Sequence[str]) -> List[str]:
    return list(dict.fromkeys(owner))


def check_uniform(owner: Sequence[str], retention: Sequence[float]) -> Dict[str, float]:
    """Return per-firm rates or raise for the restricted route."""
    if len(owner) != len(retention):
        raise ValueError("owner and retention lengths differ")
    out: Dict[str, float] = {}
    for firm in _firms(owner):
        rates = [retention[j] for j in _firm_products(owner, firm)]
        if any((not math.isfinite(r)) or r <= 0.0 for r in rates):
            raise ValueError(f"firm {firm} has non-positive/non-finite retention")
        if max(rates) - min(rates) > 1.0e-12:
            raise ValueError(f"firm {firm} has heterogeneous post-policy tariffs")
        out[firm] = rates[0]
    return out


def recover_production_costs(kappa_pre: Sequence[float], r_pre: Sequence[float]) -> List[float]:
    if len(kappa_pre) != len(r_pre):
        raise ValueError("kappa_pre and r_pre lengths differ")
    if any(r <= 0.0 for r in r_pre):
        raise ValueError("baseline retention must be positive")
    return [r * k for k, r in zip(kappa_pre, r_pre)]


def post_effective_costs(
    kappa_pre: Sequence[float],
    r_pre: Sequence[float],
    r_post: Sequence[float],
    cost_shock: Sequence[float],
) -> Tuple[List[float], List[float], List[float]]:
    """Return c_pre, c_post, kappa_post from original structural costs."""
    c_pre = recover_production_costs(kappa_pre, r_pre)
    c_post = [c * (1.0 + d) for c, d in zip(c_pre, cost_shock)]
    if any(r <= 0.0 for r in r_post):
        raise ValueError("post-policy retention must be positive")
    kappa_post = [k * (1.0 + d) * rp / rn
                  for k, d, rp, rn in zip(kappa_pre, cost_shock, r_pre, r_post)]
    return c_pre, c_post, kappa_post


def tariff_from_producer_rate(producer_rate: Sequence[float]) -> List[float]:
    return [t / (1.0 + t) for t in producer_rate]


def producer_from_tariff_rate(tau: Sequence[float]) -> List[float]:
    if any(t >= 1.0 for t in tau):
        raise ValueError("consumer-price tariff share must be below one")
    return [t / (1.0 - t) for t in tau]


def net_margin(prices: Sequence[float], c: Sequence[float], r: Sequence[float]) -> List[float]:
    return [(rr * p - cc) / (rr * p) for p, cc, rr in zip(prices, c, r)]


def net_from_gross_margin(gross_margin: Sequence[float], tau: Sequence[float]) -> List[float]:
    return [(m - t) / (1.0 - t) for m, t in zip(gross_margin, tau)]


def net_from_profit_per_gross(profit_per_gross: Sequence[float], r: Sequence[float]) -> List[float]:
    return [a / rr for a, rr in zip(profit_per_gross, r)]


def _prices_quantities(
    x: Sequence[float], demand: Demand, mode: str
) -> Tuple[List[float], List[float]]:
    if mode == "B":
        prices = list(x)
        return prices, demand.quantity(prices)
    if mode == "C":
        quantities = list(x)
        return demand.inverse(quantities), quantities
    raise ValueError("mode must be B or C")


def firm_profit(
    x: Sequence[float],
    demand: Demand,
    mode: str,
    owner: Sequence[str],
    firm: str,
    c: Sequence[float],
    retention: Sequence[float],
    fixed_cost: Mapping[str, float],
    basis: str,
) -> float:
    """One firm's payoff in original or effective-cost units.

    In effective units, c/r is used and the firm's fixed cost is divided by
    its common retention.  Thus, for uniform retention, original payoff is
    exactly r_f times this payoff, including fixed costs.  The ``naive``
    heterogeneous counterexample below intentionally calls this function with
    a non-uniform firm; its unweighted effective objective then ceases to be
    an equivalent game.
    """
    prices, quantities = _prices_quantities(x, demand, mode)
    idx = _firm_products(owner, firm)
    if not idx:
        raise ValueError(f"unknown firm {firm}")
    if basis == "original":
        variable = sum((retention[j] * prices[j] - c[j]) * quantities[j] for j in idx)
        return variable - fixed_cost.get(firm, 0.0)
    if basis == "effective":
        kappa = [c[j] / retention[j] for j in range(len(c))]
        variable = sum((prices[j] - kappa[j]) * quantities[j] for j in idx)
        rf = retention[idx[0]]
        return variable - fixed_cost.get(firm, 0.0) / rf
    if basis == "weighted_effective":
        # Algebraically identical to original variable profit, expressed as
        # r_j(p_j-kappa_j)q_j; useful for diagnosing the heterogeneous case.
        kappa = [c[j] / retention[j] for j in range(len(c))]
        variable = sum(
            retention[j] * (prices[j] - kappa[j]) * quantities[j] for j in idx
        )
        return variable - fixed_cost.get(firm, 0.0)
    raise ValueError(f"unknown profit basis {basis!r}")


def golden_max(fn: Callable[[float], float], lo: float, hi: float, iterations: int = 65) -> float:
    """Deterministic bounded one-dimensional maximizer; stdlib only."""
    if not hi > lo:
        return lo
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = hi - phi * (hi - lo)
    x2 = lo + phi * (hi - lo)
    f1, f2 = fn(x1), fn(x2)
    for _ in range(iterations):
        if f1 < f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + phi * (hi - lo)
            f2 = fn(x2)
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - phi * (hi - lo)
            f1 = fn(x1)
    return x1 if f1 >= f2 else x2


def _action_bounds(x: Sequence[float], j: int, demand: Demand, mode: str) -> Tuple[float, float]:
    if mode == "B":
        return 1.0e-5, 8.0
    if demand.kind == "logit":
        other = sum(x[k] for k in range(len(x)) if k != j)
        # Preserve a positive outside quantity for inverse Logit demand.
        return 1.0e-8, max(1.1e-7, demand.scale - other - 1.0e-8)
    return 1.0e-8, 6.0


def best_response(
    x_start: Sequence[float],
    demand: Demand,
    mode: str,
    owner: Sequence[str],
    firm: str,
    c: Sequence[float],
    retention: Sequence[float],
    fixed_cost: Mapping[str, float],
    basis: str,
    rounds: int = 18,
) -> List[float]:
    """Coordinate best response for one firm, holding rivals fixed."""
    x = list(x_start)
    products = _firm_products(owner, firm)
    for _ in range(rounds):
        before = list(x)
        for j in products:
            lo, hi = _action_bounds(x, j, demand, mode)

            def objective(value: float) -> float:
                trial = list(x)
                trial[j] = value
                try:
                    return firm_profit(
                        trial, demand, mode, owner, firm, c, retention, fixed_cost, basis
                    )
                except (ArithmeticError, ValueError, OverflowError):
                    return -1.0e300

            x[j] = golden_max(objective, lo, hi)
        if _max_diff(before, x) < 2.0e-8:
            break
    return x


def row(
    rows: List[Dict[str, object]],
    test: str,
    demand: str = "-",
    mode: str = "-",
    case: str = "-",
    max_error: float = 0.0,
    secondary_error: float = 0.0,
    tolerance: float = 1.0e-8,
    passed: bool = True,
    detail: str = "",
) -> None:
    rows.append({
        "test": test,
        "demand": demand,
        "mode": mode,
        "case": case,
        "max_error": f"{max_error:.12g}",
        "secondary_error": f"{secondary_error:.12g}",
        "tolerance": f"{tolerance:.12g}",
        "pass": int(bool(passed)),
        "detail": detail,
    })


def run(out_dir: pathlib.Path) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, object]] = []
    failures: List[Dict[str, object]] = []

    # Structural-fit effective costs and a post-merger tariff policy.  A and B
    # have different pre-policy retentions, but their products receive the same
    # post-policy retention after the merger.
    kappa_pre = [1.20, 1.50, 1.05, 1.35]
    owner_pre = ["A", "A", "B", "C"]
    r_pre = [0.78, 0.78, 0.91, 0.86]
    owner_post = ["A", "A", "A", "C"]
    r_post = [0.83, 0.83, 0.83, 0.88]
    shock = [0.10, -0.05, 0.08, 0.00]
    fixed = {"A": 0.35, "B": 0.25, "C": 0.20}
    c_pre, c_post, kappa_post = post_effective_costs(kappa_pre, r_pre, r_post, shock)

    try:
        pre_rates = check_uniform(owner_pre, r_pre)
        post_rates = check_uniform(owner_post, r_post)
        row(rows, "post_policy_uniformity", case="merger_A_B", max_error=0.0,
            detail=f"pre={pre_rates}; post={post_rates}; differing pre rates retained productwise")
    except Exception as exc:
        failures.append({"attempt": "post_policy_uniformity", "status": "failed", "error": repr(exc)})
        row(rows, "post_policy_uniformity", case="merger_A_B", passed=False, detail=repr(exc))

    try:
        bad = list(r_post)
        bad[2] = 0.91
        check_uniform(owner_post, bad)
        row(rows, "heterogeneous_policy_rejection", case="merged_A_B", passed=False,
            detail="expected rejection did not occur")
    except ValueError as exc:
        row(rows, "heterogeneous_policy_rejection", case="merged_A_B", passed=True,
            detail=f"expected rejection: {exc}")

    try:
        expected_cpost = [c * (1.0 + d) for c, d in zip(c_pre, shock)]
        expected_kpost = [k * (1.0 + d) * rp / rn
                          for k, d, rp, rn in zip(kappa_pre, shock, r_pre, r_post)]
        e1, e2 = _max_diff(c_pre, recover_production_costs(kappa_pre, r_pre)), _max_diff(c_post, expected_cpost)
        e3 = _max_diff(kappa_post, expected_kpost)
        row(rows, "cost_recovery_formula", case="pre_post", max_error=max(e1, e2, e3),
            tolerance=2.0e-14, detail="cpre=rpre*kappapre; cpost=cpre*(1+d); kappapost formula")
    except Exception as exc:
        failures.append({"attempt": "cost_recovery_formula", "status": "failed", "error": repr(exc)})
        row(rows, "cost_recovery_formula", case="pre_post", passed=False, detail=repr(exc))

    try:
        producer_rate = [0.0, 0.25, 1.0, 2.5]
        tau = tariff_from_producer_rate(producer_rate)
        recovered = producer_from_tariff_rate(tau)
        row(rows, "tariff_rate_conversion", case="producer_vs_consumer_price",
            max_error=_max_diff(producer_rate, recovered), tolerance=2.0e-14,
            detail="r=1-tau=1/(1+t); tau=t/(1+t)")
    except Exception as exc:
        failures.append({"attempt": "tariff_rate_conversion", "status": "failed", "error": repr(exc)})
        row(rows, "tariff_rate_conversion", case="producer_vs_consumer_price", passed=False, detail=repr(exc))

    logit = Demand("logit", (1.55, 1.35, 1.10, 0.95), 1.15, 10.0)
    ces = Demand("ces", (1.20, 1.05, 0.90, 0.80), 3.0, 18.0)
    base_prices = [2.35, 2.25, 2.20, 2.10]
    demands = {"logit": logit, "ces": ces}
    for name, demand in demands.items():
        try:
            q = demand.quantity(base_prices)
            recovered_p = demand.inverse(q)
            row(rows, "demand_inverse", demand=name, case="price_quantity_roundtrip",
                max_error=_max_diff(base_prices, recovered_p), tolerance=3.0e-9,
                detail="independent direct demand and inverse demand")
        except Exception as exc:
            failures.append({"attempt": f"demand_inverse:{name}", "status": "failed", "error": repr(exc)})
            row(rows, "demand_inverse", demand=name, case="price_quantity_roundtrip", passed=False, detail=repr(exc))

    # Pointwise accounting identities use the post-policy ownership.  Fixed
    # costs are included so the tested objective is not only variable profit.
    try:
        p = [2.5, 2.15, 2.30, 2.05]
        for name, demand in demands.items():
            for firm in _firms(owner_post):
                original = firm_profit(p, demand, "B", owner_post, firm, c_post,
                                        r_post, fixed, "original")
                effective = firm_profit(p, demand, "B", owner_post, firm, c_post,
                                        r_post, fixed, "effective")
                rf = r_post[_firm_products(owner_post, firm)[0]]
                row(rows, "firm_payoff_identity", demand=name, mode="B", case=firm,
                    max_error=abs(original - rf * effective), tolerance=2.0e-10,
                    detail="fixed cost included and divided by r_f in effective units")
                go = finite_gradient(
                    lambda z, f=firm: firm_profit(z, demand, "B", owner_post, f,
                                                   c_post, r_post, fixed, "original"), p
                )
                ge = finite_gradient(
                    lambda z, f=firm: firm_profit(z, demand, "B", owner_post, f,
                                                   c_post, r_post, fixed, "effective"), p
                )
                idx = _firm_products(owner_post, firm)
                grad_error = _max_abs(go[j] - rf * ge[j] for j in idx)
                row(rows, "firm_foc_scaling", demand=name, mode="B", case=firm,
                    max_error=grad_error, tolerance=3.0e-7,
                    detail="d pi_original/d choice = r_f d pi_effective/d choice")
    except Exception as exc:
        failures.append({"attempt": "payoff_and_foc_identity", "status": "failed", "error": repr(exc)})
        row(rows, "firm_payoff_identity", passed=False, detail=repr(exc))

    # Net-revenue margin and the two pointwise noiseless conversions in the
    # plan.  The generated math review records why these do not transform an
    # already-fitted noisy likelihood automatically.
    try:
        p = [2.5, 2.15, 2.30, 2.05]
        nm = net_margin(p, c_post, r_post)
        kappa = [c / r for c, r in zip(c_post, r_post)]
        nm2 = [(pp - kk) / pp for pp, kk in zip(p, kappa)]
        tau = [1.0 - r for r in r_post]
        gross = [(pp - cc) / pp for pp, cc in zip(p, c_post)]
        a = [(rr * pp - cc) / pp for pp, cc, rr in zip(p, c_post, r_post)]
        err = max(_max_diff(nm, nm2), _max_diff(nm, net_from_gross_margin(gross, tau)),
                  _max_diff(nm, net_from_profit_per_gross(a, r_post)))
        row(rows, "margin_basis_identities", case="net_revenue", max_error=err,
            tolerance=3.0e-14,
            detail="net=(rp-c)/(rp)=(p-kappa)/p; gross and profit/gross are pointwise only")
    except Exception as exc:
        failures.append({"attempt": "margin_basis_identities", "status": "failed", "error": repr(exc)})
        row(rows, "margin_basis_identities", case="net_revenue", passed=False, detail=repr(exc))

    # One/multiple leader and follower fixtures.  B/C identify the underlying
    # choice variable; F/S identify the role partition used by the check.
    # Keep B in this role fixture so BF/CF can exercise a two-firm core.  The
    # separate post-policy fixture above intentionally merges B into A.
    owner_response = ["A", "A", "B", "C"]
    role_sets = {
        "B": ("B", [], ["A", "B", "C"]),
        "C": ("C", [], ["A", "B", "C"]),
        "BF": ("B", ["A", "B"], ["C"]),       # two core leaders; fringe C
        "CF": ("C", ["A", "B"], ["C"]),
        "BS": ("B", ["A"], ["B", "C"]),      # one leader; followers B,C
        "CS": ("C", ["A"], ["B", "C"]),
    }
    for demand_name, demand in demands.items():
        for mode_name, (choice_mode, leaders, followers) in role_sets.items():
            if choice_mode == "B":
                x0 = list(base_prices)
            else:
                x0 = demand.quantity(base_prices)
            for firm in leaders + followers:
                try:
                    x_original = best_response(x0, demand, choice_mode, owner_response, firm,
                                               c_post, r_post, fixed, "original")
                    x_effective = best_response(x0, demand, choice_mode, owner_response, firm,
                                                c_post, r_post, fixed, "effective")
                    action_error = _max_diff(x_original, x_effective)
                    idx = _firm_products(owner_response, firm)
                    rf = r_post[idx[0]]
                    grad_o = finite_gradient(
                        lambda z, f=firm: firm_profit(z, demand, choice_mode, owner_response, f,
                                                       c_post, r_post, fixed, "original"),
                        x_original,
                    )
                    grad_e = finite_gradient(
                        lambda z, f=firm: firm_profit(z, demand, choice_mode, owner_response, f,
                                                       c_post, r_post, fixed, "effective"),
                        x_effective,
                    )
                    foc_scale_error = _max_abs(grad_o[j] - rf * grad_e[j] for j in idx)
                    foc_level = _max_abs(grad_o[j] for j in idx)
                    row(rows, "follower_response_equivalence", demand=demand_name,
                        mode=mode_name,
                        case=(f"leader:{firm}" if firm in leaders else f"follower:{firm}"),
                        max_error=action_error,
                        secondary_error=foc_scale_error, tolerance=3.0e-5,
                        passed=action_error <= 3.0e-5 and foc_scale_error <= 5.0e-5,
                        detail="direct physical-cost payoff versus effective-cost payoff; fixed costs retained")
                    row(rows, "response_foc", demand=demand_name, mode=mode_name,
                        case=(f"leader:{firm}" if firm in leaders else f"follower:{firm}"),
                        max_error=foc_level, secondary_error=foc_scale_error,
                        tolerance=5.0e-5,
                        passed=foc_scale_error <= 5.0e-5,
                        detail="finite-difference own-choice FOC scaling")
                except Exception as exc:
                    failures.append({"attempt": f"response:{demand_name}:{mode_name}:{firm}",
                                     "status": "failed", "error": repr(exc)})
                    row(rows, "follower_response_equivalence", demand=demand_name,
                        mode=mode_name,
                        case=(f"leader:{firm}" if firm in leaders else f"follower:{firm}"),
                        passed=False, detail=repr(exc))

    # Heterogeneous-within-firm counterexample: the direct physical objective
    # carries product-specific r_j weights, while a naive unweighted kappa game
    # drops them.  The policy is intentionally outside the restricted route.
    try:
        owner_h = ["H", "H", "R"]
        r_h = [0.58, 0.94, 0.86]
        k_h = [1.10, 1.10, 1.40]
        c_h = [k * r for k, r in zip(k_h, r_h)]
        d_h = Demand("logit", (1.90, 1.65, 1.20), 1.25, 8.0)
        p_h = [2.20, 2.10, 2.25]
        zero_fixed = {"H": 0.0, "R": 0.0}
        direct = best_response(p_h, d_h, "B", owner_h, "H", c_h, r_h,
                               zero_fixed, "original")
        naive = best_response(p_h, d_h, "B", owner_h, "H", c_h, r_h,
                              zero_fixed, "effective")
        separation = _max_diff(direct, naive)
        row(rows, "heterogeneous_counterexample", demand="logit", mode="B", case="H_two_products",
            max_error=separation, tolerance=1.0e-4, passed=separation > 1.0e-4,
            detail="r_H differs by product; unweighted kappa objective changes multiproduct optimum")
    except Exception as exc:
        failures.append({"attempt": "heterogeneous_counterexample", "status": "failed", "error": repr(exc)})
        row(rows, "heterogeneous_counterexample", demand="logit", mode="B", case="H_two_products",
            passed=False, detail=repr(exc))

    if not failures:
        failures.append({"attempt": "none", "status": "none", "error": "No numerical failures; expected scipy import was not attempted."})

    fields = ["test", "demand", "mode", "case", "max_error", "secondary_error", "tolerance", "pass", "detail"]
    with (out_dir / "oracle_validation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with (out_dir / "failed_attempts.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["attempt", "status", "error"])
        writer.writeheader()
        writer.writerows(failures)
    return rows, failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=pathlib.Path,
                        default=pathlib.Path(__file__).resolve().parent)
    args = parser.parse_args()
    rows, failures = run(args.out_dir)
    failed_rows = [r for r in rows if not bool(int(r["pass"]))]
    print(f"wrote {len(rows)} validation rows to {args.out_dir}")
    print(f"validation_failures={len(failed_rows)} retained_attempts={len(failures)}")
    return 1 if failed_rows else 0


if __name__ == "__main__":
    raise SystemExit(main())
