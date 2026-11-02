These compact artifacts record the executed baseline profile. The scripts are
execution snapshots using temporary output paths and locally ignored fit archives.
Use a matching installed baseline package and adjust paths when reproducing them.
Heavy RDS objects and compiler artifacts are deliberately excluded.

The original `helmert_timing_v2.csv` is an independent R reference calculation.
Subsequent `helmert_matched_*` artifacts compare compiled native kernels with
matched compiler commands. Those are fixed-point target/gradient measurements,
not HMC ESS/second evidence. Standalone GQ timing includes context setup.

The panel fixture runner records full normalized targets with and without
Jacobians, full gradients, finite differences and complete fixed-draw output.
Its first initialization attempt dropped dimensions on length-one vectors;
`helmert_panel_initial_runner_failures.csv` retains those failed runs. The repaired
runner assigns through `[]` to preserve the Stan initialization shapes.
