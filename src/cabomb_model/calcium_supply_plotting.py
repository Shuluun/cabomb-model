"""Figures for conditional budgets and normalized time; no sample rate is assumed."""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


def _save(fig: plt.Figure, path: Path, caption: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0.16, 1, 0.94))
    fig.text(0.07, 0.025, caption, fontsize=8)
    fig.savefig(path.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def plot_supply_results(static: list[dict], time: list[dict], points: list[dict], interpolated: list[dict], output: Path) -> None:
    """Render deterministic mathematical curves and separate literature interpolation."""
    with mpl.rc_context({"svg.fonttype": "none", "svg.hashsalt": "calcium-no-experiment-v0.1",
                         "font.size": 10, "axes.spines.top": False, "axes.spines.right": False}):
        fig, ax = plt.subplots(figsize=(7.8, 5.5))
        for fraction, color, style in ((1.0, "#4477AA", "-"), (0.5, "#AA3377", "--"), (0.1, "#777777", ":")):
            rows = [row for row in static if row["dissolved_fraction"] == fraction]
            ax.loglog([r["q_pg_per_pl"] for r in rows], [r["supply_equivalent_um"] for r in rows],
                      color=color, linestyle=style, label=f"f = {fraction:g}", linewidth=2)
        ax.set_xlabel("q = initial internalized pure CaCO3 mass / cytosolic reference volume (pg/pL)", fontsize=8)
        ax.set_ylabel("Carrier calcium supply equivalent (µM)")
        ax.set_title("Conditional stoichiometric budget")
        ax.legend(title="Dissolved fraction")
        _save(fig, output / "static_budget", "q axis is an adjustable display window, not a biological range or recommended dose.\n"
              "Cumulative carrier supply normalized by reference volume; not free or retained cytosolic calcium.")

        fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
        x = [r["x_half_times"] for r in time]
        axes[0].plot(x, [r["dissolved_fraction"] for r in time], color="#4477AA", linewidth=2)
        axes[1].plot(x, [r["normalized_rate_per_half_time"] for r in time], color="#AA3377", linewidth=2)
        axes[0].set_title("a  Cumulative normalized supply")
        axes[1].set_title("b  Normalized supply rate")
        axes[0].set_ylabel("C_eq(t) / C_eq_complete = F(x)")
        axes[1].set_ylabel("dF/dx (per normalized half time)")
        axes[0].set_ylim(0, 1)
        axes[1].set_ylim(bottom=0)
        for ax in axes:
            ax.set_xlabel("x = t / t_half (dimensionless)")
            ax.set_xlim(0, max(x))
        _save(fig, output / "normalized_time", "Assumed first-order dissolution; one batch already internalized at t = 0; fixed volume; complete endpoint.\n"
              "No hour-scale half time is assigned. Actual rate = (C_eq_complete / t_half) × dF/dx.")

        fig, ax = plt.subplots(figsize=(7.8, 5.2))
        ax.plot([r["time_h"] for r in interpolated], [r["reported_dissolution_percent"] for r in interpolated],
                color="#777777", label="Piecewise-linear interpolation")
        ax.scatter([float(r["time_h"]) for r in points], [float(r["reported_dissolution_percent"]) for r in points],
                   color="#4477AA", zorder=3, label="Published Table 2 values")
        ax.set_title("Literature system only: Lee et al. (2015), N-Cal")
        ax.set_xlabel("Time in the original assay (h)")
        ax.set_ylabel("Reported dissolution (%)")
        ax.set_xlim(float(points[0]["time_h"]), float(points[-1]["time_h"]))
        ax.set_ylim(0, 1)
        ax.legend()
        _save(fig, output / "literature_interpolation", "Simulated gastric fluid, initial pH 1.5, 37 °C, 5 mg/mL. Three published points; no invented error bars.\n"
              "DOI: 10.2147/IJN.S79403. No extrapolation, endpoint rescaling or sample half-time fitting.")
