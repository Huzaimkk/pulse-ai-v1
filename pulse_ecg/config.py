"""Load config.yaml — the single home of every tunable constant."""

from functools import lru_cache
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "config.yaml"


@lru_cache(maxsize=None)
def load(path: str | Path = DEFAULT) -> dict:
    with open(path, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError(f"{path}: expected a mapping at top level")
    return cfg


def path(key: str, cfg_path: str | Path = DEFAULT) -> Path:
    """A `paths:` entry resolved against the directory holding config.yaml."""
    cfg_path = Path(cfg_path)
    return (cfg_path.parent / load(cfg_path)["paths"][key]).resolve()
