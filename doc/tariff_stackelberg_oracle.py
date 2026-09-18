#!/usr/bin/env python3
"""Bounded independent one-leader Stackelberg checks.

This reference-only script imports the local standard-library primitives in
tariff_uniform_oracle.py.  For each fixed leader action it solves the two
followers' Nash response using repeated *actual* physical-payoff best
responses.  It then optimizes the multiproduct leader's reduced payoff through
that solved response and repeats the calculation in effective-cost units.

The four deliberately small cases are Logit/CES x Bertrand/Cournot.  They use
one leader (A, which owns products 1 and 2) and two followers (B and C).  No
production package, R, scipy, HMC, or network operation is used.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
from typing import Callable, Dict, List, Mapping, Sequence, Tuple

import tariff_uniform_oracle as ref


OWNER = ["A", "A", "B", "C"]
LEADER = "A"
FOLLOWERS = ["B", "C"]
RETENTION = [0.83, 0.83, 0.83, 0.88]
KAPPA = [1.20, 1.50, 1.05, 1.35]
R_PRE = [0.78, 0.78, 0.91, 0.86]
SHOCK = [0.10, -0.05, 0.08, 0.00]
FIXED = {"A": 0.35, "B": 0.25, "C": 0.20}


def _costs() -> Tuple[List[float], List[float]]:
    c_pre, c_post, _ = ref.post_effective_costs(KAPPA, R_PRE, RETENTION, SHOCK)
    return c_post, RETENTION[:]


def _max_abs(values: Sequence[float]) -> float:
    return max((abs(float(v)) for v in values), default=0.0)


def _max_diff(x: Sequence[float], y: Sequence[float]) -> float:
    return _max_abs([a - b for a, b in zip(x, y)])


def _scale(x: Sequence[float]) -> float:
    return max(1.0, _max_abs(x))


def _replace(x: Sequence[float], j: int, value: float) -> List[float]:
    out = list(x)
    out[j] = value
    return out


def finite_gradient_scaled(
    fn: Callable[[Sequence[float]], float],
    x: Sequence[float],
    indices: Sequence[int],
    relative_step: float = 3.0e-3,
) -> List[float]:
    """Central finite differences with a scale-aware step for nested solves."""
    out: List[float] = []
    for j in indices:
        h = max(3.0e-4, abs(float(x[j])) * relative_step)
        def central(step: float) -> float:
            plus = _replace(x, j, x[j] + step)
            minus = _replace(x, j, max(1.0e-7, x[j] - step))
            return (fn(plus) - fn(minus)) / (plus[j] - minus[j])
        coarse = central(h)
        fine = central(0.5 * h)
        # Remove the leading O(h^2) truncation term from the nested reduced
        # profit derivative without changing the independent solve.
        out.append((4.0 * fine - coarse) / 3.0)
    return out


def follower_conditional_foc(
    x: Sequence[float],
    demand: ref.Demand,
    mode: str,
    j: int,
    c: Sequence[float],
    retention: Sequence[float],
    basis: str,
) -> float:
    """Analytic own FOC for a single-product follower.

    Demand and inverse-demand derivatives are written here independently of
    the package and of the finite-difference optimizer.  The caller only uses
    this for B/C followers with one product each.
    """
    prices, quantities = ref._prices_quantities(x, demand, mode)
    p, q = prices[j], quantities[j]
    rr, cc = retention[j], c[j]
    if mode == "B":
        if demand.kind == "logit":
            dqdx = -demand.slope * q * (1.0 - q / demand.scale)
        elif demand.kind == "ces":
            revenue_share = p * q / demand.scale
            dqdx = q / p * (-demand.slope + (demand.slope - 1.0) * revenue_share)
        else:
            raise ValueError("unsupported demand for follower FOC")
        if basis == "original":
            return rr * q + (rr * p - cc) * dqdx
        kappa = cc / rr
        return q + (p - kappa) * dqdx

    if mode == "C":
        if demand.kind == "logit":
            outside = demand.scale - sum(quantities)
            dpdx = -(1.0 / q + 1.0 / outside) / demand.slope
        elif demand.kind == "ces":
            gamma = demand.slope
            quality = demand.delta_or_quality[j]
            z = p * (q / quality) ** (1.0 / gamma)
            outside_revenue_share = (
                demand.outside_price ** (1.0 - gamma) * z ** gamma / demand.scale
            )
            dpdx = (
                -p / (gamma * q)
                - (gamma - 1.0) * p * p
                / (gamma * demand.scale * (1.0 + (gamma - 1.0) * outside_revenue_share))
            )
        else:
            raise ValueError("unsupported demand for follower FOC")
        if basis == "original":
            return rr * p - cc + rr * q * dpdx
        kappa = cc / rr
        return p - kappa + q * dpdx
    raise ValueError("mode must be B or C")


def precise_follower_response(
    x_start: Sequence[float],
    demand: ref.Demand,
    mode: str,
    firm: str,
    c: Sequence[float],
    retention: Sequence[float],
    basis: str,
    rounds: int = 12,
    bisection_iterations: int = 65,
) -> Tuple[List[float], float, bool]:
    """Bisection best response for a one-product follower."""
    products = ref._firm_products(OWNER, firm)
    if len(products) != 1:
        raise ValueError("precise follower FOC requires one product per follower")
    j = products[0]
    x = list(x_start)
    bracketed = True
    for _ in range(rounds):
        lo, hi = ref._action_bounds(x, j, demand, mode)

        def foc(value: float) -> float:
            return follower_conditional_foc(
                _replace(x, j, value), demand, mode, j, c, retention, basis,
            )

        flo, fhi = foc(lo), foc(hi)
        if not (math.isfinite(flo) and math.isfinite(fhi)):
            raise ValueError("non-finite follower FOC bracket")
        if flo <= 0.0:
            x[j] = lo
            bracketed = False
        elif fhi >= 0.0:
            x[j] = hi
            bracketed = False
        else:
            for _ in range(bisection_iterations):
                mid = 0.5 * (lo + hi)
                fmid = foc(mid)
                if fmid > 0.0:
                    lo = mid
                else:
                    hi = mid
            x[j] = 0.5 * (lo + hi)
    final_foc = abs(follower_conditional_foc(x, demand, mode, j, c, retention, basis))
    lo, hi = ref._action_bounds(x, j, demand, mode)
    curvature_step = max(1.0e-5, abs(x[j]) * 1.0e-5)
    left = max(lo + 1.0e-9, x[j] - curvature_step)
    right = min(hi - 1.0e-9, x[j] + curvature_step)
    if not bracketed or not (lo < x[j] < hi) or not right > left:
        raise ValueError("follower FOC solution is not an interior bracketed response")
    curvature = (
        follower_conditional_foc(_replace(x, j, right), demand, mode, j, c, retention, basis)
        - follower_conditional_foc(_replace(x, j, left), demand, mode, j, c, retention, basis)
    ) / (right - left)
    if not math.isfinite(curvature) or curvature >= 0.0:
        raise ValueError("follower own FOC lacks negative local curvature")
    return x, final_foc, bracketed


def follower_nash(
    start: Sequence[float],
    demand: ref.Demand,
    mode: str,
    c: Sequence[float],
    retention: Sequence[float],
    basis: str,
    sweeps: int = 5,
    response_rounds: int = 5,
    golden_iterations: int = 20,
    residual_check: bool = True,
) -> Tuple[List[float], int, float]:
    """Solve followers' conditional Nash response by actual best responses."""
    x = list(start)
    used = 0
    for used in range(1, sweeps + 1):
        before = list(x)
        for firm in FOLLOWERS:
            response, _, _ = precise_follower_response(
                x, demand, mode, firm, c, retention, basis,
                rounds=response_rounds, bisection_iterations=golden_iterations + 30,
            )
            for j in ref._firm_products(OWNER, firm):
                x[j] = response[j]
        if _max_diff(before, x) <= 2.0e-9 * _scale(x):
            break

    if not residual_check:
        return x, used, float("nan")

    # A final independent residual check solves each follower once more and
    # records the largest action change.  This is reported, not hidden.
    residual = 0.0
    for firm in FOLLOWERS:
        response, _, _ = precise_follower_response(
            x, demand, mode, firm, c, retention, basis,
            rounds=max(response_rounds, 8), bisection_iterations=max(golden_iterations + 35, 55),
        )
        j = ref._firm_products(OWNER, firm)[0]
        # Report the FOC at the returned Nash profile x.  The response solve
        # above is only the independent best-response residual check; using
        # its own returned FOC here would falsely certify any starting x.
        residual = max(
            residual,
            abs(follower_conditional_foc(x, demand, mode, j, c, retention, basis)),
        )
    return x, used, residual


def reduced_value(
    candidate: Sequence[float],
    demand: ref.Demand,
    mode: str,
    c: Sequence[float],
    retention: Sequence[float],
    basis: str,
    start_followers: Sequence[float],
    strict: bool = False,
) -> Tuple[float, List[float], float]:
    # Start followers from the supplied baseline while preserving the trial
    # leader action.  Overwriting these leader entries with the baseline would
    # make the reduced objective spuriously constant in leader choices.
    follower_start = list(start_followers)
    for j in ref._firm_products(OWNER, LEADER):
        follower_start[j] = candidate[j]
    if strict:
        solved, _, residual = follower_nash(
            follower_start, demand, mode, c, retention, basis,
            sweeps=8, response_rounds=8, golden_iterations=30,
            residual_check=True,
        )
    else:
        solved, _, residual = follower_nash(
            follower_start, demand, mode, c, retention, basis,
            residual_check=False,
        )
    value = ref.firm_profit(
        solved, demand, mode, OWNER, LEADER, c, retention, FIXED, basis,
    )
    return value, solved, residual


def leader_reduced_solve(
    start: Sequence[float],
    demand: ref.Demand,
    mode: str,
    c: Sequence[float],
    retention: Sequence[float],
    basis: str,
    leader_sweeps: int = 5,
    golden_iterations: int = 28,
) -> Tuple[List[float], List[float], float, int]:
    """Optimize A's own actions while fully resolving B/C at every trial."""
    x = list(start)
    leader_indices = ref._firm_products(OWNER, LEADER)
    start_followers = list(start)
    for used in range(1, leader_sweeps + 1):
        before = list(x)
        for j in leader_indices:
            lo, hi = ref._action_bounds(x, j, demand, mode)

            def objective(value: float) -> float:
                trial = _replace(x, j, value)
                try:
                    result, _, _ = reduced_value(
                        trial, demand, mode, c, retention, basis, start_followers,
                    )
                    return result
                except (ArithmeticError, ValueError, OverflowError):
                    return -1.0e300

            x[j] = ref.golden_max(objective, lo, hi, iterations=golden_iterations)
        if _max_diff(before, x) <= 2.0e-8 * _scale(x):
            break

    final, followers, residual = reduced_value(
        x, demand, mode, c, retention, basis, start_followers, strict=True,
    )
    # reduced_value returns the fully solved profile in `followers`; retaining
    # the separate leader vector makes the reduced-game construction explicit.
    return followers, followers, residual, used


def reaction_derivatives(
    leader_profile: Sequence[float],
    demand: ref.Demand,
    mode: str,
    c: Sequence[float],
    retention: Sequence[float],
    basis: str,
) -> List[List[float]]:
    """Finite-difference d follower-Nash-actions / d leader-actions."""
    leader_indices = ref._firm_products(OWNER, LEADER)
    follower_indices = [j for j, firm in enumerate(OWNER) if firm in FOLLOWERS]
    out = [[0.0 for _ in leader_indices] for _ in follower_indices]
    for col, j in enumerate(leader_indices):
        # The larger step keeps a 1e-8 response-solve residual from being
        # magnified into a false derivative discrepancy.
        h = max(1.0e-4, abs(leader_profile[j]) * 1.0e-4)
        plus = _replace(leader_profile, j, leader_profile[j] + h)
        minus = _replace(leader_profile, j, max(1.0e-7, leader_profile[j] - h))
        _, pplus, _ = reduced_value(
            plus, demand, mode, c, retention, basis, leader_profile, strict=True,
        )
        _, pminus, _ = reduced_value(
            minus, demand, mode, c, retention, basis, leader_profile, strict=True,
        )
        denominator = plus[j] - minus[j]
        for row, product in enumerate(follower_indices):
            out[row][col] = (pplus[product] - pminus[product]) / denominator
    return out


def flatten(matrix: Sequence[Sequence[float]]) -> List[float]:
    return [value for row in matrix for value in row]


def add_row(
    rows: List[Dict[str, object]],
    case: str,
    demand: str,
    mode: str,
    check: str,
    error: float,
    scale: float,
    tolerance: float,
    passed: bool,
    detail: str,
) -> None:
    rows.append({
        "case": case,
        "demand": demand,
        "mode": mode,
        "check": check,
        "max_error": f"{error:.12g}",
        "scale": f"{scale:.12g}",
        "tolerance": f"{tolerance:.12g}",
        "pass": int(bool(passed)),
        "detail": detail,
    })


def run(out_dir: pathlib.Path) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    c_post, retention = _costs()
    demands = {
        "logit": ref.Demand("logit", (1.55, 1.35, 1.10, 0.95), 1.15, 10.0),
        "ces": ref.Demand("ces", (1.20, 1.05, 0.90, 0.80), 3.0, 18.0),
    }
    prices = [2.35, 2.25, 2.20, 2.10]
    rows: List[Dict[str, object]] = []
    failures: List[Dict[str, object]] = []

    for demand_name, demand in demands.items():
        for mode in ("B", "C"):
            case = "one_leader_two_followers_multiproduct_leader"
            try:
                start = prices if mode == "B" else demand.quantity(prices)
                leader_search = {
                    "leader_sweeps": 6 if mode == "B" else 5,
                    "golden_iterations": 42 if mode == "B" else 28,
                }
                physical, _, physical_residual, physical_sweeps = leader_reduced_solve(
                    start, demand, mode, c_post, retention, "original", **leader_search,
                )
                effective, _, effective_residual, effective_sweeps = leader_reduced_solve(
                    start, demand, mode, c_post, retention, "effective", **leader_search,
                )
                leader_idx = ref._firm_products(OWNER, LEADER)
                follower_idx = [j for j, firm in enumerate(OWNER) if firm in FOLLOWERS]
                action_scale = _scale(physical)

                # Guard against the reduced-value construction silently
                # discarding a trial leader action.  The leader coordinates
                # must pass through exactly, and at least one follower or the
                # reduced profit must react to a finite leader perturbation.
                trial_step = max(2.0e-4, abs(physical[leader_idx[0]]) * 2.0e-4)
                trial = _replace(physical, leader_idx[0], physical[leader_idx[0]] + trial_step)
                trial_value, trial_profile, _ = reduced_value(
                    trial, demand, mode, c_post, retention, "original", start,
                )
                leader_trial_error = abs(trial_profile[leader_idx[0]] - trial[leader_idx[0]])
                follower_trial_change = _max_diff(
                    [physical[j] for j in follower_idx],
                    [trial_profile[j] for j in follower_idx],
                )
                base_value = ref.firm_profit(
                    physical, demand, mode, OWNER, LEADER, c_post, retention, FIXED, "original",
                )
                profit_trial_change = abs(trial_value - base_value)
                if leader_trial_error > 1.0e-12 * action_scale:
                    raise AssertionError(
                        "trial leader action was overwritten before follower solving"
                    )
                if max(follower_trial_change, profit_trial_change) <= 1.0e-8:
                    raise AssertionError(
                        "trial leader action produced no follower/reduced-profit response"
                    )
                add_row(
                    rows, case, demand_name, mode, "leader_trial_passthrough_and_sensitivity",
                    leader_trial_error, action_scale, 1.0e-12 * action_scale,
                    True,
                    f"follower_change={follower_trial_change:.6g}; profit_change={profit_trial_change:.6g}",
                )
                leader_action_error = _max_diff(
                    [physical[j] for j in leader_idx],
                    [effective[j] for j in leader_idx],
                )
                follower_action_error = _max_diff(
                    [physical[j] for j in follower_idx],
                    [effective[j] for j in follower_idx],
                )
                action_tol = 1.0e-5 * action_scale
                add_row(
                    rows, case, demand_name, mode, "leader_reduced_action_equivalence",
                    leader_action_error, action_scale, action_tol,
                    leader_action_error <= action_tol,
                    f"A owns products 1,2; physical sweeps={physical_sweeps}; effective sweeps={effective_sweeps}",
                )
                add_row(
                    rows, case, demand_name, mode, "follower_nash_action_equivalence",
                    follower_action_error, action_scale, action_tol,
                    follower_action_error <= action_tol,
                    "B,C solved by repeated actual best responses at every leader trial",
                )

                residual_scale = 1.0
                residual_tol = 1.0e-10
                add_row(
                    rows, case, demand_name, mode, "physical_follower_nash_residual",
                    physical_residual, residual_scale, residual_tol,
                    physical_residual <= residual_tol,
                    "final independent physical-payoff best-response residual",
                )
                add_row(
                    rows, case, demand_name, mode, "effective_follower_nash_residual",
                    effective_residual, residual_scale, residual_tol,
                    effective_residual <= residual_tol,
                    "final independent effective-cost best-response residual",
                )

                # Reduced-profit FOC scaling is assessed at each independently
                # optimized profile.  Both gradients are recomputed through a
                # fresh fully solved follower response.
                def physical_reduced(z: Sequence[float]) -> float:
                    return reduced_value(
                        z, demand, mode, c_post, retention, "original", start,
                        strict=True,
                    )[0]

                def effective_reduced(z: Sequence[float]) -> float:
                    return reduced_value(
                        z, demand, mode, c_post, retention, "effective", start,
                        strict=True,
                    )[0]

                # Compare the two payoff gradients at the same solved leader
                # profile.  Comparing gp(physical) to ge(effective) would mix
                # the positive-scaling identity with equilibrium-search error.
                gp = finite_gradient_scaled(physical_reduced, physical, leader_idx)
                ge = finite_gradient_scaled(effective_reduced, physical, leader_idx)
                rf = retention[leader_idx[0]]
                foc_error = _max_abs(a - rf * b for a, b in zip(gp, ge))
                foc_scale = max(1.0, _max_abs(gp))
                foc_tol = 1.0e-5 * foc_scale
                foc_level_error = max(_max_abs(gp), _max_abs(ge))
                foc_level_scale = max(1.0, foc_level_error)
                foc_level_tol = 1.0e-5 * foc_level_scale
                add_row(
                    rows, case, demand_name, mode, "leader_reduced_foc_residual",
                    foc_level_error, foc_level_scale, foc_level_tol,
                    foc_level_error <= foc_level_tol,
                    "absolute finite-difference own-choice gradient at reduced optimum",
                )
                add_row(
                    rows, case, demand_name, mode, "leader_reduced_foc_scaling",
                    foc_error, foc_scale, foc_tol, foc_error <= foc_tol,
                    "d reduced physical payoff = r_A d reduced effective payoff",
                )

                # Likewise, evaluate both follower reaction maps at one common
                # leader profile; equilibrium action differences are reported
                # separately above.
                dp = reaction_derivatives(
                    physical, demand, mode, c_post, retention, "original",
                )
                de = reaction_derivatives(
                    physical, demand, mode, c_post, retention, "effective",
                )
                reaction_error = _max_diff(flatten(dp), flatten(de))
                reaction_scale = max(1.0, _max_abs(flatten(dp)))
                reaction_tol = 1.0e-5 * reaction_scale
                add_row(
                    rows, case, demand_name, mode,
                    "follower_reaction_derivative_equivalence", reaction_error,
                    reaction_scale, reaction_tol, reaction_error <= reaction_tol,
                    "finite difference of fully re-solved follower Nash actions",
                )
            except Exception as exc:
                failures.append({
                    "case": case, "demand": demand_name, "mode": mode,
                    "status": "failed", "error": repr(exc),
                })
                add_row(
                    rows, case, demand_name, mode, "case_execution", math.inf, 1.0,
                    0.0, False, repr(exc),
                )

    failed_checks = [r for r in rows if not bool(int(r["pass"]))]
    for failed in failed_checks:
        failures.append({
            "case": failed["case"],
            "demand": failed["demand"],
            "mode": failed["mode"],
            "status": "validation_failed",
            "error": (
                f"{failed['check']}: max_error={failed['max_error']} "
                f"> tolerance={failed['tolerance']}"
            ),
        })
    if not failures:
        failures.append({
            "case": "none", "demand": "-", "mode": "-", "status": "none",
            "error": "No numerical failures; all four mandatory cases executed.",
        })

    fields = ["case", "demand", "mode", "check", "max_error", "scale", "tolerance", "pass", "detail"]
    with (out_dir / "full_stackelberg_validation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with (out_dir / "full_stackelberg_failures.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case", "demand", "mode", "status", "error"])
        writer.writeheader()
        writer.writerows(failures)
    return rows, failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=pathlib.Path,
                        default=pathlib.Path(__file__).resolve().parent)
    args = parser.parse_args()
    rows, _ = run(args.out_dir)
    failed = [r for r in rows if not bool(int(r["pass"]))]
    print(f"wrote {len(rows)} full Stackelberg validation rows to {args.out_dir}")
    print(f"validation_failures={len(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
