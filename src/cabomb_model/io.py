"""Small CSV and filesystem helpers."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Mapping


def ensure_dir(path: str | Path) -> Path:
    """Create a directory and return it as a Path."""

    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_csv(path: str | Path, rows: Iterable[Mapping[str, object]]) -> Path:
    """Write dictionaries to CSV using their first-row field order."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows_list = list(rows)
    if not rows_list:
        output_path.write_text("", encoding="utf-8")
        return output_path
    fieldnames = list(rows_list[0].keys())
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_list)
    return output_path

