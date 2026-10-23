# Reproducing the release work

The public commands implement the forecasting methods and the fixed-input
intracellular response analysis. The earlier 14-curve comparison used published
source workbooks; obtain them from the articles listed in `REFERENCES.md`.
`data_provenance.json` records every curve's worksheet, cell mapping, aggregation,
time conversion, split and original source-file checksum. It contains no observed
response values.

## Historical literature panel

1. Obtain each study's Source Data workbook and retain its original file.
   The source checksum identifies the snapshot used in the earlier run.
2. For C02, C09 and C16, use the inclusive row bounds, time column and value
   columns in `data_provenance.json`. Read numerical rows, group any repeated
   numerical time, and average the listed same-time numerical value columns.
   These columns are not assumed to represent longitudinal replicate identities.
   Convert time using `time_to_h`; retain the recorded percentage denominator.
3. Use the recorded structural zero and the first `n_train` positive unique
   times as the early window. Remaining positive times form the late window.
   The source row and selection notes identify blocks intentionally excluded
   as duplicated formulations.
4. For the two CaCO₃ series, use Figure1 C189:G200: column C is hours, D/E are
   the pH 7.4 mean/SD, and F/G the pH 5.0 mean/SD. Use the explicitly listed
   early and late times. Retain author means without synthesizing replicates.
5. Save `early.csv`, `future_times.csv` and `late.csv` using the schema in
   `../docs/run_guide.md`, then run the separate `forecast` and `evaluate`
   commands. Keep source metadata and checksums with the local files.
6. For the study-equal summary, average each method's curve RMSE values within
   each study, then average those four study means. Count curves whose late
   RMSE exceeds 10 percentage points separately. An unresolved curve keeps
   the aggregate incomplete rather than being silently dropped.

The saved historical findings remain a retrospective development result. The
package's synthetic example is the quickest full execution check, while these
source mappings support reproducing the literature comparison with obtained
workbooks. The old strict historical replay check and its tolerance are retained
in `findings.md`; they are distinct from rerunning the current implementation.

## Delivery-model link

```bash
python -m cabomb_release bridge --output-dir outputs/bridge
```

The command records its model parameters, internalization input and waiting
response settings. Mass-ledger, first-order identity and three-grid quadrature
checks are saved with the resulting CSVs and figure. It uses the bundled model
scenario described in the historical protocol.

## Current delivery analyses

The default local sensitivity and Monte Carlo analyses concern the sequential
OCN–YSA design. They exclude `alpha_ysa_direct_access` and `k_ysa_on`, whose
routes are inactive for that combination. The dual-versus-YSA-only feasibility
grid retains `k_ysa_on` because it affects YSA-only: five parameters at three
levels give 243 sampled settings. Every grid row records the model modes.
The optional parallel route remains part of the separate mechanism comparison.

The software release removes an inactive grid dimension from the earlier
729-row analysis. The unique default grid results and their win fractions are
preserved; the uncertainty samples are regenerated after filtering the inactive
parameters. New outputs are kept separately from the historical local check.
