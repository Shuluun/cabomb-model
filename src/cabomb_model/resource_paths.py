"""Locate resources without searching unrelated parent repositories."""
from pathlib import Path
from contextlib import contextmanager
from contextvars import ContextVar
from typing import Iterator

_registry: ContextVar[Path | None] = ContextVar("cabomb_parameter_registry", default=None)

def resource_path(name: str) -> Path:
    return Path(__file__).resolve().parent / "resources" / name

def selected_registry_path() -> Path:
    """Return the explicitly selected registry or the bundled default."""
    return _registry.get() or resource_path("parameters_master.csv")

@contextmanager
def parameter_registry(path: str | Path | None) -> Iterator[None]:
    """Select one registry for a command and restore the caller's selection."""
    selected = Path(path).expanduser().resolve() if path is not None else None
    if selected is not None and not selected.is_file():
        raise FileNotFoundError(f"Parameter registry not found: {selected}")
    token = _registry.set(selected)
    try:
        yield
    finally:
        _registry.reset(token)

def source_root() -> Path | None:
    candidate = Path(__file__).resolve().parents[2]
    if (candidate / "pyproject.toml").is_file() and (candidate / "src/cabomb_model/resource_paths.py").resolve() == Path(__file__).resolve():
        return candidate
    return None
