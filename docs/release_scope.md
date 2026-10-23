# Release scope

This standalone personal project is prepared from the existing CABOMB model
and the subsequent release-method study. The source working files are retained
separately from this distribution.

The public candidate includes the amount-based ODE model, parameter registries,
tests, optional Streamlit interface, empirical forecasting methods, an original
synthetic example, and an age-since-internalization response analysis. Model
equations are unchanged; resource loading and run manifests are made portable.
The explicit `--registry` option selects edited parameter files and records
their snapshots. Default sequential sensitivity/sampling excludes inactive
routes, and the dual-vs-YSA-only grid removes a repeated inactive dimension.

The literature-dependent legacy `calcium_no_experiment.py` workflow is outside
this release. Its source datasets and historical workstation orchestration are
not included. The independent stoichiometric calculator remains available.

The 14-curve comparison is reported as retrospective development, not a new
held-out validation. Published raw data, extracted literature curves and figures
based on those curves are excluded from this preview. Users can import legally
obtained series using the documented CSV schema. This public package alone does
not reconstruct all historical literature analyses.

Historical protocol files document the prior experiment, including its original
workflow restrictions. They are not instructions to install or run this package.
New packaging checks are recorded separately from historical scientific results.
