# Running the model and release analysis

Install and activate the environment using the README commands. Core operation
needs NumPy, SciPy, pandas, Matplotlib and openpyxl. Streamlit is optional.

## Delivery model

```bash
python -m cabomb_model run-baseline --model-version v1.3 --output-dir outputs/baseline
python -m cabomb_model run-analysis --model-version v1.3 --output-dir outputs/analysis
python -m cabomb_model build-report --model-version v1.3 --output-dir outputs/report_run --report-path reports/model_report.md --html-path reports/model_report.html
```

The full analysis is longer than the baseline. Its report command reruns the
analysis. V1.1 and V1.2 remain available through the same version flag.
The model commands can overwrite files in their selected output directory;
choose a new directory when retaining a previous run.

To change model parameters, copy `data/processed/parameters_master.csv`, edit
the values and evidence notes, and explicitly select the edited file:

```bash
python -m cabomb_model validate-registry --registry data/local/my_parameters.csv
python -m cabomb_model run-baseline --model-version v1.3 --registry data/local/my_parameters.csv --output-dir outputs/custom_baseline
```

`--registry` is available for baseline, analysis and report commands. Each saves
`registry_used.csv` and `registry_source.json` with the selected file's SHA256.
Without the flag, the installed package uses its bundled default; editing the
readable CSV in `data/processed/` alone does not change an installed run. Python
users can pass a path to `v13_default_parameters(path)` or use the scoped
`cabomb_model.resource_paths.parameter_registry(path)` context manager.

## Your own release series

Prepare three UTF-8 CSVs with these columns:

| File | Columns | Meaning |
| --- | --- | --- |
| `early.csv` | `study_id,curve_id,time_h,release_percent` | Known early observations |
| `future_times.csv` | `study_id,curve_id,time_h` | Times to predict, without values |
| `late.csv` | `study_id,curve_id,time_h,release_percent` | Held-out values used only for scoring |

Each early curve must begin with a justified structural `(0,0)` and have at
least four strictly positive observation times. Do not invent a zero for an
assay with pre-existing release. Percentages must be finite and between 0 and
100. Duplicate study/curve/time keys are rejected. Future times must be strictly
after the early window; late times must exactly match them. Source aggregation
and normalization should be fixed before inspecting prediction errors.

Place personal files under `data/local/`, which is ignored by Git, and run:

```bash
python -m cabomb_release forecast --early data/local/early.csv --future data/local/future_times.csv --output-dir outputs/my_forecast
python -m cabomb_release evaluate --run-dir outputs/my_forecast --late data/local/late.csv
```

Forecasting saves every optimizer start, rolling decisions and predictions
before the evaluation command reads late values. Integrity hashes detect
changes to forecast artifacts before evaluation. They record provenance and
do not by themselves prove that a study was blinded.

The methods are eight bounded empirical families, persistence and an anchored
last-three-point slope, with three selection rules. The rolling rules use the
same early target points across methods; insufficient folds trigger an explicit
fallback. Full details are in `research/protocol/specification.md`.

## Intracellular waiting-time analysis

```bash
python -m cabomb_release bridge --output-dir outputs/bridge
```

This fixes the 200 nm, medium-EPR/MPS V1.3 scenario and convolves its
internalization flux with Weibull cumulative waiting responses. Shape values
0.5, 1 and 2 share a four-hour mean, an assumed scenario tied to the original
first-order rate of 0.25 per hour. The ledger is normalized carrier-equivalent
response, with `gross + pending = cumulative internalized input`.
This fixed reference command uses bundled model parameters. Custom delivery
experiments use the model CLI's `--registry` option above.

## Conditional calcium equivalent

```bash
python -m cabomb_model.calcium_equivalent --mass-pg 0.1 --dissolved-fraction 0.5 --volume-pl 2 --output-dir outputs/calcium_equivalent
```

These values are example conditions. This is a stoichiometric supply-equivalent
calculation, not a measurement of free cytosolic calcium.

## Input and output locations

Parameter registries are under `data/processed/`; package resources contain
matching CSV copies. Synthetic inputs are written under the selected run's
`inputs/` subdirectory. Forecast output includes predictions, decisions,
optimization diagnostics and a lock file; late metrics are under `evaluation/`.
Bridge output includes its parameter snapshot, fluxes, response time courses
and numerical checks. README files contain instructions rather than datasets.
