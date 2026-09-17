# Finite-difference cancellation diagnostic

The original seven-case gate remains in `extra_helmert_comparison.csv`:
six pass and CES auction exceeds the unchanged absolute tolerance 1e-10.
At the common point its normalized/raw targets agree exactly, the full
analytic gradients differ by 8.88e-16, and constrained/GQ outputs differ by
2.78e-17. A central finite difference with h=1e-5 differs by 7.11e-10.

The independent diagnostic repeats both kernels in separate R processes
on a frozen step grid from 1e-2 to 1e-5, retaining every value. For the
offending quality coordinate, subtraction of target values around -70
amplifies a roughly 1.42e-14 rounding difference by 1/(2h). This explains
the observed 7.11e-10 difference without asserting that a failed gate passed.
At h=1e-5, both kernels' finite differences also miss their own analytic
gradients by more than 1e-10: the original measurement is insufficiently
accurate for that tolerance.

The secondary test uses fourth-order Richardson differentiation,
`(4 D(h/2)-D(h))/3`, at h=1e-3. It removes the central rule's leading h^2
truncation error while reducing cancellation compared with h=1e-5. The
usual balance is O(h^4) truncation plus O(machine_epsilon * abs(target)/h)
rounding; the step grid demonstrates the resulting accuracy here.
All three tested coordinates agree between kernels exactly, and each
kernel agrees with its own analytic derivatives within 8.87e-12.

`resolve_fd_gate.R` preserves the original failure, requires every other
target/gradient/schema/centering gate, and applies the same 1e-10 threshold
to both the refined between-kernel difference and the two analytic errors.
It resolves this identified measurement failure; it does not provide a
generic exception for future failed cases. No prior, HMC, economic equation,
solver tolerance or default changes. These are numerical derivative checks,
not sampler or model-comparison evidence.
