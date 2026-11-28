"""What made an output directory: code commit, environment, tool versions and parameters.

Every pipeline script calls ``write_provenance`` for the directories it writes, which
gets ``provenance_<step>.json``. The environment lock lists the conda packages
(name-version-build) of the running interpreter's prefix.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TOOLS = {
    "minimap2": ["minimap2", "--version"],
    "samtools": ["samtools", "--version"],
    "sniffles": ["sniffles", "--version"],
    "miniasm": ["miniasm", "-V"],
    "racon": ["racon", "--version"],
    "prodigal": ["prodigal", "-v"],
    "hmmsearch": ["hmmsearch", "-h"],
    "blastn": ["blastn", "-version"],
}


def _run(cmd: list[str]) -> str:
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    lines = (result.stdout or result.stderr).strip().splitlines()
    if cmd[0] == "hmmsearch":  # the version is in the banner
        lines = [line for line in lines if line.startswith("# HMMER")]
    return lines[0].strip() if lines else ""


def write_provenance(out_dirs, step: str, parameters: dict, extra: dict | None = None) -> None:
    """Write ``provenance_<step>.json`` into each of ``out_dirs``."""
    meta = Path(sys.prefix) / "conda-meta"
    info = {
        "step": step,
        "time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code_commit": _run(["git", "-C", str(REPO), "rev-parse", "HEAD"]) or "unknown",
        "code_changed_files": subprocess.run(["git", "-C", str(REPO), "status", "--short"], capture_output=True,
                                             text=True).stdout.splitlines(),
        "command": sys.argv,
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in parameters.items()},
        "python": sys.version.split()[0],
        "tools": {name: _run(cmd) for name, cmd in TOOLS.items()},
        "environment_prefix": sys.prefix,
        "environment_lock": sorted(p.stem for p in meta.glob("*.json")) if meta.is_dir() else [],
        **(extra or {}),
    }
    for directory in out_dirs:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        (directory / f"provenance_{step}.json").write_text(json.dumps(info, indent=2, default=str))
