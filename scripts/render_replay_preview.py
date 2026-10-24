"""Render a GitHub GIF from the earlier viewer's saved data; no ODE rerun."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[1]
VIEWER = ROOT / "visualizations/organ_replay"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qa-dir", type=Path, help="Optional directory for four still frames.")
    args = parser.parse_args()
    data = json.loads((VIEWER / "data.json").read_text(encoding="utf-8"))
    points = np.asarray(data["points"], dtype=float)
    times = points[:, 0]
    assert points.shape == (400, 4)
    assert np.isfinite(points).all() and (points[:, 1:] >= 0).all()
    assert np.all(np.diff(times) > 0) and times[0] == 0 and times[-1] == 48
    bg, panel, text, muted, cyan, gold = (
        "#101b2b", "#19283b", "#e9f0f4", "#9fb0bf", "#72d6df", "#f4d28a"
    )
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.2), dpi=95, facecolor=bg)
    fig.subplots_adjust(left=.07, right=.975, bottom=.23, top=.64, wspace=.32)
    fig.text(.07, .91, "CABOMB | Carrier and cargo time replay", color=text,
             fontsize=19, weight="bold")
    fig.text(.07, .82, "Earlier V1.2 baseline  /  100 nm  /  OCN–YSA", color=muted, fontsize=11)
    clock = fig.text(.975, .82, "", ha="right", color=gold, fontsize=13, weight="bold")
    titles = ["Liver-associated carrier", "Tumor-associated carrier", "Active intracellular cargo"]
    lines, dots, markers, values = [], [], [], []
    for i, ax in enumerate(axes):
        ax.set_facecolor(panel)
        ax.set_title(titles[i], color=text, fontsize=11, pad=37, loc="left")
        ax.set(xlim=(0, 48), ylim=(0, data["displayCeilings"][i]), xticks=[0, 16, 32, 48])
        ax.set_xlabel("Time (h)", color=muted)
        ax.set_ylabel("% initial carrier" if i < 2 else "% initial cargo", color=muted, fontsize=9)
        ax.tick_params(colors=muted, labelsize=9)
        for spine in ax.spines.values():
            spine.set_color("#35465a")
        ax.grid(axis="y", color="#35465a", alpha=.6, linewidth=.6)
        ax.plot(times, points[:, i+1], color=muted, alpha=.23, lw=1.3)
        lines.append(ax.plot([], [], color=cyan, lw=2.4)[0])
        dots.append(ax.plot([], [], "o", color=gold, ms=5)[0])
        markers.append(ax.axvline(0, color=gold, alpha=.6, lw=.9, ls="--"))
        values.append(ax.text(0, 1.05, "", transform=ax.transAxes, color=gold,
                              fontsize=18, weight="bold"))
    fig.text(.07, .09, "Model-generated regional amounts. Tumor carrier includes internalized carrier.",
             color=muted, fontsize=9)
    fig.text(.07, .045, "Preview exported from the original 2026-10-17 interactive viewer data · CABOMB",
             color=muted, fontsize=8)

    def update(t: float):
        clock.set_text(f"{t:04.1f} h")
        selected = times < t
        for i in range(3):
            value = float(np.interp(t, times, points[:, i+1]))
            lines[i].set_data(np.r_[times[selected], t], np.r_[points[selected, i+1], value])
            dots[i].set_data([t], [value])
            markers[i].set_xdata([t, t])
            values[i].set_text(f"{value:.3f}%")
        return [clock, *lines, *dots, *markers, *values]

    frames = np.r_[np.linspace(0, 48, 61), np.repeat(48., 10)]
    animation = FuncAnimation(fig, update, frames=frames, interval=100, blit=False)
    animation.save(VIEWER / "time_replay.gif", writer=PillowWriter(fps=10))
    if args.qa_dir:
        args.qa_dir.mkdir(parents=True, exist_ok=True)
        for t in [0, 8, 24, 48]:
            update(t)
            fig.savefig(args.qa_dir / f"time_{t:02d}h.png", facecolor=bg)
    plt.close(fig)
    print(f"Saved {VIEWER / 'time_replay.gif'}")


if __name__ == "__main__":
    main()
