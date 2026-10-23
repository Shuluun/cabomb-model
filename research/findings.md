# Research findings

The follow-up work connects early release forecasting with the delivery model's
intracellular response. It uses a shared cumulative-response function, with a
different clock in each setting: time since the in-vitro experiment began, or
age since a carrier entered a cell.

## Early measurements and late release

The retrospective development panel contains 14 curves from four studies.
Study-equal mean late RMSE was 7.830354 percentage points for the original
four-family early-SSE selector, 13.758546 for rolling selection within those
four families, and 8.339084 for the expanded rolling selector. Respectively,
2, 4 and 3 of the 14 curves exceeded a late RMSE of 10 percentage points.

The useful result is that better early fitting and a larger function library
did not consistently improve late prediction. Some individual curves benefited
from a flexible plateau or a simple slope baseline. Selecting among them from
a short prefix remains the harder problem. This panel had already been seen
during development, so a new independent panel is the next useful test.

## Connection to CABOMB

For an internalization flux J(s), cumulative response is
`Q(t) = integral_0^t J(s) F(t-s) ds`. Each cohort starts its own waiting clock
when it enters the cell. An exponential F reproduces the fixed-rate first-order
degradation response; changing its shape explores other waiting profiles.

With the original 200 nm, medium-EPR/MPS V1.3 input held fixed, three Weibull
shapes with the same four-hour mean changed the early response timing. Across
the sampled 4, 12, 24, 48 and 72-hour times, the ranking was YSA > OCN–YSA >
unmodified > OCN for this response metric. By 72 hours the shape-dependent
span was below 1% relative to the exponential result for all four designs.
This is a fixed-scenario structural comparison; it motivates checking early
measurements and robustness across other scenarios.

## Numerical record

The historical bridge passed its mass ledger, non-negativity, first-order
identity and quadrature checks. The largest exponential-reference difference
was about 7.43e-8 normalized units, and the largest mid/fine quadrature difference
about 1.15e-7. New package runs are checked separately in
`docs/validation_report.md`.

One strict historical forecast-replay check failed: the maximum deviation was
5.620903849e-6 percentage points against a frozen 1e-6 tolerance. The recorded
diagnosis traced most of this numerical difference to retaining versus omitting
a structural-zero residual row in optimization. The original threshold was
retained. This replay check and the current method comparison answer different
questions; the failed check is not relabeled as a pass.

## Next research step

The most useful extension is a prospectively specified test on new time series,
followed by uncertainty propagation from the release response to delivery
metrics. Matched CABOMB dissolution, internalization and cargo-release assays
would then test whether the same response function can connect those processes.
The convolution itself is established mathematics; a publishable contribution
would need new evidence about prediction or decision robustness.

The historical literature series are not shipped with this preview. The
synthetic command demonstrates the complete portable workflow, and the bridge
can be rerun using the bundled model parameters.
