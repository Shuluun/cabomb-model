"""Generated artifacts derived from the parameter registry."""

from __future__ import annotations

from pathlib import Path

from cabomb_model.io import ensure_dir
from cabomb_model.parameters import REQUIRED_REGISTRY_COLUMNS, load_registry

PREVIEW_FILENAME = "parameters_master_preview.png"
PREVIEW_ROWS = 12
MAX_CELL_CHARS = 54


def _clip_cell(value: object, *, max_chars: int = MAX_CELL_CHARS) -> str:
    text = str(value).encode("ascii", "replace").decode("ascii")
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 3] + "..."


def write_parameter_registry_preview(
    output_dir: str | Path = "outputs",
    *,
    registry_path: str | Path | None = None,
) -> Path:
    """Render a deterministic PNG preview of the parameter registry."""

    records = load_registry(registry_path)
    if not records:
        raise ValueError("Cannot render parameter registry preview from an empty registry.")

    columns = [column for column in REQUIRED_REGISTRY_COLUMNS if column in records[0]]
    preview_records = records[:PREVIEW_ROWS]
    cell_text = [
        [_clip_cell(record.get(column, "")) for column in columns]
        for record in preview_records
    ]

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    row_count = len(preview_records) + 1
    fig_width = 24
    fig_height = max(4.0, 0.34 * row_count)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis("off")

    column_widths = [
        0.055,
        0.075,
        0.16,
        0.055,
        0.055,
        0.055,
        0.065,
        0.075,
        0.1,
        0.08,
        0.075,
        0.085,
        0.085,
        0.09,
        0.075,
        0.055,
        0.13,
        0.12,
    ]
    total_width = sum(column_widths[: len(columns)])
    normalized_widths = [width / total_width for width in column_widths[: len(columns)]]

    table = ax.table(
        cellText=cell_text,
        colLabels=columns,
        cellLoc="left",
        colLoc="left",
        colWidths=normalized_widths,
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(4.6)
    table.scale(1.0, 1.45)

    for (row, _col), cell in table.get_celld().items():
        cell.set_edgecolor("#D9E2F3")
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#1F4E79")
            cell.get_text().set_color("#FFFFFF")
            cell.get_text().set_weight("bold")
        elif row % 2 == 0:
            cell.set_facecolor("#EAF2F8")
        else:
            cell.set_facecolor("#FFFFFF")

    preview_path = ensure_dir(Path(output_dir) / "figures") / PREVIEW_FILENAME
    fig.savefig(
        preview_path,
        dpi=160,
        bbox_inches="tight",
        pad_inches=0.08,
        metadata={"Software": "cabomb_model"},
    )
    plt.close(fig)
    return preview_path
