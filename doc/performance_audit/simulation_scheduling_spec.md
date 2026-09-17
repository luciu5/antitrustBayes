# Posterior simulation scheduling: equivalence gate

This is an implementation specification, not an accepted optimization or a
performance result. The existing scheduler remains the reference until the
measurements below pass.

## Task identity and randomness

Select posterior indices and realize the quality/ordered-k draw plan once,
before dispatch. A task is `(unit_id, draw_position)`, where draw position
indexes that fixed plan. Scheduling must never redraw a quality or a role,
change posterior weights, or alter the order of returned rows. Preserve the
existing draw-parallel scheduler as the default during the experiment.

Compare three schedulers: draws within each unit; units with sequential draws;
and a flattened unit-by-draw task list. Construct immutable unit inputs once.
Return the same result/error record for every task, including failures. Do not
drop failed draws or silently normalize policy summaries over a different set.

## Warm starts

The unit scheduler may use the previous successful post-policy equilibrium as
the next task's `priceStart`, passed through the public `antitrust::simulate()`
argument. Use it only when labels, product count, active subset, and finite
price requirements match. A previous draw supplies a numerical initial value,
never costs, demand, ownership, or any other economic state.

An explicit user `priceStart` takes precedence. Do not reorder draws by
structural proximity in the first experiment. If later tried, permute task
execution only and restore the original draw positions before assembly.

A solver error, non-finite equilibrium, or residual exceeding the existing
accepted criterion triggers a canonical-start retry with the same parameters,
scenario and unchanged solver controls. Retain the initial failure and retry
outcome in diagnostics. Warnings indicating nonconvergence require the same
treatment even when a solver returns a vector. If multiple equilibria produce
materially different policy outputs, hold the optimization rather than choose
the faster result.

Native hard-role simulations use `coordination` and need a separately verified
public warm-start capability. Do not modify coordination or pass unsupported
arguments under this task. Scheduling can still change independently of warm
starts; report the limitation explicitly if its public API cannot consume one.

## Measurement and acceptance

Use 200 fixed posterior positions for both a small single-market fixture and
a multi-market panel. Include sampled/collapsed quality and exact BF/CF/BS/CS
where supported. Reuse the same draw plan for observed and latent-quality
anchoring, each compared against its own unchanged reference.

Separate extraction, unit preparation, structural construction, policy solve,
and result assembly. Count solver calls, iterations and function evaluations;
report canonical retries, failed tasks, and process-tree peak PSS when practical.
Compare product prices/shares, consumer and producer surplus, and total welfare
under merger, +10% costs, and combined policies at existing numerical tolerances.
Small fixtures may favor serial dispatch because worker overhead dominates.
Retain timing regressions; only integrate a scheduling change when its measured
benefit justifies the additional code.
