"""Upload cleaned batch outputs to HuggingFace dataset repo.

Includes per batch:
- is_formatter_output/   IS records JSONs (annotated, guide, raw)
- partial_circle_output/ Partial circle calls
- sniffles_output/       Sniffles VCFs and tables
- metadata_for_sniffles.tsv (if exists)
- chunks/*_all.tsv       Sample metadata

Excludes:
- reference_genomes/     (recoverable from NCBI)
- chunks/chunk_*/sniffles_output, chunks/chunk_*/circle_output (redundant or large)
- BAMs, FASTQs (already cleaned)
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repository root
from cycle.site_config import data_root  # noqa: E402

from huggingface_hub import HfApi


def upload_batch(repo_id: str, batch_dir: Path, batch_name: str):
    """Upload a single batch directory to HF."""
    api = HfApi()

    # Build allow patterns relative to repo path (we'll upload to /<batch_name>/)
    allow_patterns = [
        f"{batch_name}/is_formatter_output/**",
        f"{batch_name}/partial_circle_output/**",
        f"{batch_name}/sniffles_output/**",
        f"{batch_name}/metadata_for_sniffles.tsv",
        f"{batch_name}/chunks/*_all.tsv",
    ]
    ignore_patterns = [
        "**/sra_downloads/**",
        "**/alignments/**",
        "**/*.bam", "**/*.bam.bai",
        "**/assembly/**",
        "**/reference_genomes/**",
        "**/chunks/chunk_*/circle_output/**",  # redundant per-chunk circle output
    ]

    # Build a temp staging tree by symlinking (avoids copying)
    staging = Path(os.environ.get("TMPDIR", "/tmp")) / "hf_staging" / batch_name
    if staging.exists():
        import shutil
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    # The Hub client does not follow directory symlinks, so every file is linked on its own.
    # Sniffles calls may sit at the batch level or under chunks/chunk_N/.
    sources = [batch_dir / sub for sub in ("is_formatter_output", "partial_circle_output", "sniffles_output")]
    sources += sorted((batch_dir / "chunks").glob("chunk_*/sniffles_output"))
    for src in sources:
        for path in src.rglob("*"):
            if path.is_file() and "assembly" not in path.relative_to(src).parts and not path.name.endswith((".bam", ".bai")):
                dst = staging / path.relative_to(batch_dir)
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.symlink_to(path)

    meta = batch_dir / "metadata_for_sniffles.tsv"
    if meta.exists():
        (staging / "metadata_for_sniffles.tsv").symlink_to(meta)

    # Copy chunks/_all.tsv files
    chunks_src = batch_dir / "chunks"
    if chunks_src.exists():
        chunks_dst = staging / "chunks"
        chunks_dst.mkdir(exist_ok=True)
        for tsv in chunks_src.glob("*_all.tsv"):
            (chunks_dst / tsv.name).symlink_to(tsv)

    print(f"Uploading {batch_name} from {staging}...")
    api.upload_folder(
        folder_path=str(staging),
        repo_id=repo_id,
        repo_type="dataset",
        path_in_repo=batch_name,
        ignore_patterns=ignore_patterns,
    )
    print(f"  Done.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-id", default=os.environ.get("HF_DATASET_REPO"),
                        help="HuggingFace dataset repo (or set HF_DATASET_REPO)")
    parser.add_argument("--base-dir", default=str(data_root()))
    parser.add_argument("--batches", nargs="+", required=True,
                        help="Batch numbers (e.g. 6 7 8) or names (e.g. batch_006)")
    args = parser.parse_args()
    if not args.repo_id:
        parser.error("--repo-id or HF_DATASET_REPO is required")

    api = HfApi()
    # Create repo if needed
    try:
        api.create_repo(args.repo_id, repo_type="dataset", exist_ok=True, private=True)
        print(f"Repo {args.repo_id} ready")
    except Exception as e:
        print(f"Repo creation: {e}")

    base = Path(args.base_dir)
    for b in args.batches:
        if b.startswith("batch_"):
            name = b
        else:
            name = f"batch_{int(b):03d}"
        bdir = base / name
        if not bdir.exists():
            print(f"Skipping {name}: not found")
            continue
        upload_batch(args.repo_id, bdir, name)


if __name__ == "__main__":
    main()
