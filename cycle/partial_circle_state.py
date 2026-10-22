"""Completion records for resumable partial-circle jobs (no tool dependencies)."""

import hashlib
import json
from pathlib import Path


def file_stamp(path):
    path = Path(path)
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def input_signature(records, fastq, parameters):
    return {
        "version": 1,
        "records_sha256": hashlib.sha256(Path(records).read_bytes()).hexdigest(),
        "fastq": file_stamp(fastq),
        "parameters": parameters,
    }


def output_paths(directory, sample_id):
    directory = Path(directory)
    return [directory / f"{sample_id}_partial_circle_{suffix}" for suffix in (
        "reads.tsv", "summary.tsv", "summary.json",
    )]


def is_complete(directory, sample_id, signature):
    try:
        state = json.loads((Path(directory) / "complete.json").read_text())
        return (state["input"] == signature
                and state["outputs"] == [file_stamp(p) for p in output_paths(directory, sample_id)]
                and all(p.stat().st_size > 0 for p in output_paths(directory, sample_id)))
    except (OSError, ValueError, KeyError, TypeError):
        return False


def mark_complete(directory, sample_id, signature):
    paths = output_paths(directory, sample_id)
    if not all(p.is_file() and p.stat().st_size > 0 for p in paths):
        raise ValueError(f"Incomplete partial-circle outputs for {sample_id}")
    if not isinstance(json.loads(paths[-1].read_text()), list):
        raise ValueError(f"Invalid partial-circle summary for {sample_id}")
    state = {"input": signature, "outputs": [file_stamp(p) for p in paths]}
    pending = Path(directory) / "complete.json.tmp"
    pending.write_text(json.dumps(state, indent=2), encoding="utf-8")
    pending.replace(Path(directory) / "complete.json")
