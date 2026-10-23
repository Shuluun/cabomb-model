# Data sources

The model uses the numeric parameter registries in `data/processed/`. Each row
records its source, units, bounds and evidence status. Assumed values and
project-supplied reference values remain labeled in the registry. Identical CSVs
are bundled inside the installed Python package so runs do not depend on the
working directory. The XLSX file is a readable companion to the master CSV.

The prior release study used two CaCO₃ series and twelve series from three
other studies. The publications below own the underlying measurements.

| Study | Primary publication | Recorded access / license information |
| --- | --- | --- |
| C02 | [10.1038/s41467-026-72589-0](https://doi.org/10.1038/s41467-026-72589-0) | article CC BY-NC-ND 4.0; author Source Data publicly attached; no separate workbook license located; local research extraction only, no redistribution permission inferred |
| C09 | [10.1038/s41467-023-36330-5](https://doi.org/10.1038/s41467-023-36330-5) | article CC BY 4.0; author Source Data attached; local snapshots retained with attribution |
| C16 | [10.1038/s41467-022-32043-3](https://doi.org/10.1038/s41467-022-32043-3) | article CC BY 4.0; author Source Data attached; local attribution retained |
| CaCO3_2024 | [10.1038/s41467-024-54402-y](https://doi.org/10.1038/s41467-024-54402-y) | Article CC BY-NC-ND 4.0. Separate source-data permission not established. Not bundled. |

These source attachments and extracted curves are not bundled in this preview.
The NC-ND entries have no independently established source-workbook permission
in the project record. For the other entries, access notes are preserved from
the earlier source audit. Obtain files from the primary article's Source Data
section and check their current terms before redistributing them.

To analyze an obtained series, record its DOI, workbook/sheet/cell range,
denominator, time conversion, aggregation and file SHA256. Convert it into the
CSV schema in `docs/run_guide.md`, using hours and the reported release percent.
Do not treat percent values from different assays as one biological quantity.
For CaCO₃ Figure 1g, the original DCNP/DCNP-Ab label and denominator ambiguity
remain recorded in the research notes.

The shipped executable demonstration is mathematically generated. Its formula
and inputs are saved separately from forecasts when the example runs. No
synthetic points are presented as laboratory observations.
