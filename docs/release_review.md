# Release review

An independent AI agent reviewed the code, documentation, research outputs and
local release archives on 2026-10-23. The following issues were corrected in
the release copy:

| Area | Change |
| --- | --- |
| Parameter selection | Added `--registry`; runs save the selected CSV and its SHA256. |
| Parameter documentation | Synchronized four pulmonary-retention upper bounds with the executable registry. |
| Analysis scope | Removed inactive routes from default sequential-dual sensitivity and sampling; removed the repeated inactive grid dimension. |
| Generated interpretation | Feasibility dimensions and win fractions now come from actual output rather than fixed prose. |
| Reproduction | Added source-cell mappings, aggregation, time conversion, splits and checksums for all 14 historical literature curves. |
| Numerical portability | Separated strict historical trajectory replay, same-backend output checks and independent-solver convergence; retained the original failed Ubuntu replay record. |

The reviewer independently compared the findings with the saved curve-level
metrics and bridge outputs, inspected archive structure and links, and ran
the new registry and analysis-scope tests. No unresolved major or medium
release issue remained after those fixes. The final full-suite execution and
installed-package checks are recorded in `validation_report.md`.
The subsequent numerical review compared tighter LSODA and Radau solutions
and verified the unchanged core-model and parameter hashes recorded in
`../research/numerical_reference.json`.

This review is an AI-assisted software and research-record check. It is separate
from external human peer review or biological experiments.
