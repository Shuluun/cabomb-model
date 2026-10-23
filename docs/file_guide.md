# Code file guide

| File | Role |
| --- | --- |
| `src/cabomb_model/__init__.py` | CABOMB first-version PBPK-inspired ODE model. |
| `src/cabomb_model/__main__.py` | Module entry point for `python -m cabomb_model`. |
| `src/cabomb_model/analysis.py` | Analysis workflows for CABOMB V1. |
| `src/cabomb_model/calcium_equivalent.py` | Conditional calcium supply equivalent, independent of the delivery ODE. |
| `src/cabomb_model/calcium_supply_plotting.py` | Figures for conditional budgets and normalized time; no sample rate is assumed. |
| `src/cabomb_model/calcium_supply_time.py` | Conditional single-batch supply; not cytosolic free or retained total calcium. |
| `src/cabomb_model/cli.py` | Command line interface for the CABOMB V1 model. |
| `src/cabomb_model/designs.py` | Carrier design switches for the frozen CABOMB V1 model. |
| `src/cabomb_model/fluxes.py` | Shared instantaneous flux definitions for CABOMB V1.1. |
| `src/cabomb_model/io.py` | Small CSV and filesystem helpers. |
| `src/cabomb_model/metrics.py` | Mass ledgers and output metrics for CABOMB simulations. |
| `src/cabomb_model/model.py` | Frozen CABOMB V1 ODE right-hand side and rate functions. |
| `src/cabomb_model/model_modes.py` | Explicit model-version switches for CABOMB simulations. |
| `src/cabomb_model/parameters.py` | Parameter registry loading and executable parameter validation. |
| `src/cabomb_model/plotting.py` | SVG plotting helpers for CABOMB generated figures. |
| `src/cabomb_model/registry_artifacts.py` | Generated artifacts derived from the parameter registry. |
| `src/cabomb_model/reporting.py` | Report generation helpers. |
| `src/cabomb_model/resource_paths.py` | Locate resources without searching unrelated parent repositories. |
| `src/cabomb_model/scenarios.py` | EPR, MPS, IFP, and pH scenarios for exploratory model runs. |
| `src/cabomb_model/solver.py` | Numerical integration for the CABOMB V1 ODE model. |
| `src/cabomb_model/states.py` | State definitions for the frozen CABOMB V1 ODE system. |
| `src/cabomb_model/v13_analysis.py` | Analysis workflows for the optional V1.3 explicit-circulation mode. |
| `src/cabomb_model/v13_circulation.py` | V1.3 explicit-circulation ODE, solver, and mass-ledger metrics. |
| `src/cabomb_model/v13_plotting.py` | SVG plotting helpers for the V1.3 explicit-circulation mode. |
| `src/cabomb_model/v13_reporting.py` | Markdown report generation for the optional V1.3 explicit-circulation mode. |
| `src/cabomb_model/v13_states.py` | State definitions for the optional V1.3 explicit-circulation mode. |
| `src/cabomb_release/__init__.py` | Early-curve forecasting and assumed intracellular response kernels. |
| `src/cabomb_release/__main__.py` | Package entry point or exports. |
| `src/cabomb_release/bridge.py` | Fixed V1.3 input with age-since-internalization response scenarios. |
| `src/cabomb_release/cli.py` | Portable forecasting commands with separate fitting and late evaluation. |
| `src/cabomb_release/methods.py` | Empirical bounded percentage forecasts; no mechanistic parameter inference. |
| `src/cabomb_release/plotting.py` | Figures for the model-generated waiting-time comparison. |
| `app/streamlit_app.py` | Streamlit interface for CABOMB V1.1, V1.2, and V1.3 exploratory modes. |
| `tests/test_analysis.py` | Regression and behavior tests. |
| `tests/test_calcium_equivalent.py` | Unit, domain and CLI checks for the independent stoichiometric calculator. |
| `tests/test_calcium_supply_time.py` | Physical accounting and explicit time-unit contracts, without sample assumptions. |
| `tests/test_cli.py` | Regression and behavior tests. |
| `tests/test_fluxes.py` | Regression and behavior tests. |
| `tests/test_model_core.py` | Regression and behavior tests. |
| `tests/test_model_modes.py` | Regression and behavior tests. |
| `tests/test_mps_v12.py` | Regression and behavior tests. |
| `tests/test_parameter_registry.py` | Regression and behavior tests. |
| `tests/test_registry_selection.py` | User-selected registries must affect runs and remain traceable. |
| `tests/test_release_mathematics.py` | Mathematical tests retained from the retrospective methods study. |
| `tests/test_release_workflow.py` | Information boundaries and standalone resource behavior. |
| `tests/test_reporting.py` | Regression and behavior tests. |
| `tests/test_sequential_targeting_v12.py` | Regression and behavior tests. |
| `tests/test_v13_ablation.py` | Integration checks for the scoped comparison and its saved evidence. |
| `tests/test_v13_analysis_scope.py` | Guard against inactive sweep dimensions and data-independent reporting. |
| `tests/test_v13_explicit_circulation.py` | Tests for the optional V1.3 explicit-circulation mode. |

Numeric resources under `src/*/resources/` are included in the wheel. The
master CSV resource matches the readable registry in `data/processed/`.
The literature-dependent legacy no-experiment workflow is outside this release;
its original files remain in the working project.
