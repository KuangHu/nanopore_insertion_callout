"""Site-specific locations for the Cycle code, read from environment variables.

The shell entry points load these variables from ``config/site.env`` (or ``config/igi.env``) via
``config/load_site.sh``. Python code calls the helpers below instead of hardcoding paths. When a
variable is unset, the helpers fall back to the IGI locations used before the 2026-11 migration, so
existing IGI workflows keep working unchanged.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

IGI_DEFAULTS = {
    "CYCLE_DATA_ROOT": "/groups/rubin/projects/kuang/out/IS_cycle",
    "CYCLE_HMM_DIR": "/home/kuangh/scripts/IS110/hmm",
    "CYCLE_ECOR_DIR": "/groups/rubin/projects/kuang/out/IS110/ECOR_batch",
}


def _env_path(name: str) -> Path:
    value = os.environ.get(name) or IGI_DEFAULTS[name]
    return Path(value).expanduser()


def data_root() -> Path:
    """Root of the IS_cycle output tree (``CYCLE_DATA_ROOT``)."""
    return _env_path("CYCLE_DATA_ROOT")


def hmm_dir() -> Path:
    """Directory holding DEDD.hmm and Tnp20.hmm (``CYCLE_HMM_DIR``; else ``<repo>/data/hmm`` if present)."""
    if os.environ.get("CYCLE_HMM_DIR"):
        return Path(os.environ["CYCLE_HMM_DIR"]).expanduser()
    local = REPO_ROOT / "data" / "hmm"
    if (local / "DEDD.hmm").exists():
        return local
    return Path(IGI_DEFAULTS["CYCLE_HMM_DIR"])


def ecor_dir() -> Path:
    """Root of the ECOR side-project inputs and outputs (``CYCLE_ECOR_DIR``)."""
    return _env_path("CYCLE_ECOR_DIR")


def batch_metadata_dir() -> Path:
    """Directory with batch_NNN.tsv sample lists (``CYCLE_BATCH_METADATA_DIR``; default ``<repo>/data/batches``)."""
    value = os.environ.get("CYCLE_BATCH_METADATA_DIR")
    return Path(value).expanduser() if value else REPO_ROOT / "data" / "batches"
