# Synthetic forecasting example

```bash
python -m cabomb_release synthetic --output-dir outputs/synthetic
```

The command creates a deterministic, noise-free Weibull series, separates its
early inputs from late observations, fits and saves forecasts, then evaluates
them and draws a figure. The formula, units and synthetic status are recorded
in `inputs/generation.json`. This example checks the workflow and illustrates
the task; experimental performance comes from a separate evaluation.
