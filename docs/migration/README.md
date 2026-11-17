# Migration: IGI cluster -> Lawrencium (2026-11)

This directory documents how the code was made portable and how the data were packaged.
Data packages, checksums, verification reports and chat exports are **not** in git; they live in
the migration bundle on IGI (`/groups/rubin/projects/kuang/migration_20261117/`, see its
`MIGRATION_INDEX.md`).

## What changed in the code (branch `migration/lawrencium`)

- All cluster-specific values moved to a site file: `config/igi.env` (the old hardcoded values, used
  when nothing else is configured) and `config/lawrencium.env.example` (template for
  `config/site.env`). `config/load_site.sh` loads them; `cycle/site_config.py` is the Python side.
- `submit_one_batch.sh` and `submit_all.sh` read code, data, log, environment and scheduler settings
  from the site file. With the IGI defaults, `submit_one_batch.sh` emits byte-identical `sbatch`
  arguments to the previous version (checked with a stub `sbatch`, 119/119 arguments).
- Standalone SLURM scripts (`scripts/*.sh`, `sideproject/**/*.sh`) dropped partition, qos, node
  exclusions and log paths from their `#SBATCH` headers; submit them with `bin/cycle-sbatch`.
- Python scripts locate the repository from their own path and read data/asset locations from
  `cycle.site_config`; the IS110 HMM directory is `CYCLE_HMM_DIR`.
- `upload_batch_to_hf.py` no longer defaults to a private dataset name (`--repo-id` or
  `HF_DATASET_REPO` is required).
- Environment inventory: `envs/` (exact pinned lists, portable specs, tool versions).
- Every original cluster-specific reference and what was done with it: `path_review.tsv`
  (263 references in 50 files, classified by purpose).

See `docs/lawrence_setup.md` for configuring and running on Lawrencium.

## Packaged components (in the migration bundle)

| Component | Content |
|---|---|
| IS_cycle results | `IS_cycle` after the 2026-11-17 cleanup as per-batch / per-folder `tar.zst`, one deduplicated reference-genome archive, member lists, `restore_refs.py`, `verify_restore.py` |
| IS_cycle external link targets | dereferenced copies of 220 symlinks in `is110_top20/` that pointed into a shared NCBI genome database |
| cycle_rna_analysis | the sibling analysis project, per top-level folder |
| IS110_downstream | the sibling project plus `out/is110_analysis` (target of two of its symlinks) |
| assets | IS110 HMMs and `Pfam_transposase.hmm`, with original paths and checksums; checksum of the shared Pfam-A release |
| code snapshot | byte-exact tarball of the IGI working tree before any migration change (including untracked files) |
| environment inventory | `envs/` exports |
| chat export | Claude Code transcripts and memory for this project (private) |

## Verification approach

Each component was checked independently of the code that built it: SHA-256 and `zstd -t` of every
archive, extraction into a new directory, rebuilding per-batch reference files, re-hashing the
restored tree with a separate inventory tool, and comparing every path, type, size, SHA-256, mode,
symlink target and file mtime against an inventory of the source. A before/after source
inventory confirms the sources did not change while being packaged.

## Milestones

1. **Ready to transfer**: all components built and verified on IGI (integrity + independent restore).
2. **Restored and usable on Lawrencium**: checksums verified after transfer, archives extracted,
   `verify_restore.py` and the inventory comparison pass on the target, environment rebuilt,
   `config/site.env` filled in, test suite and a one-sample pipeline run succeed.
3. **Source eligible for cleanup**: only after milestone 2 is confirmed; the IGI sources are kept
   untouched until then.

`READS_FQ_REPRODUCIBILITY.md` documents the legacy `reads.fq` files removed before packaging and
how they can be recreated.
