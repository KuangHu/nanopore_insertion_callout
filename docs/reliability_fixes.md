# Reliability fixes from the 2026-10-22 review

This change addresses the five reproduced orchestration/input bugs and the
reproduced F3/F4 negative cases from the review of `b5101b4`.

## Changes

- FASTQ lookup accepts string or Path directories, requires a complete sample
  identifier, excludes unfinished downloads/directories and rejects ambiguous
  fallback matches. Standard filename precedence is retained.
- `submit_one_batch.sh` gives each chunk its own `reference_genomes/`, including
  `resolve_status.tsv`, downloads and indices. This trades some duplicate
  reference storage for independent jobs without shared-file races.
- All three minimap2/samtools mapping pipes check both processes, drain stderr
  without a pipe deadlock, and publish BAM files only after success. Indexed
  mappings also require successful indexing. Alignment jobs exit nonzero when
  a sample fails or lacks its required reference/FASTQ.
- Both submission scripts use `run_partial_circle_batch.py`. Failed children
  are reported and make the batch fail after remaining samples are attempted.
  The single-sample command owns the sample subdirectory, avoiding the former
  `sample/sample/` nesting. Phase 4 also waits for its own Phase 1 download job.
- Partial-circle completion requires all three result files and an atomic
  `complete.json` marker tied to the input records, FASTQ file metadata and
  parameters. Empty valid results are successful; missing/invalid inputs and
  mapping errors are failures. Interrupted or changed-input jobs rebuild their
  reference/BAM intermediates. FASTQ identity uses path, size and mtime, not a
  full digest of potentially very large reads files.
- Formatter localization rejects clipped, incomplete, secondary, supplementary,
  reference-skipped or unknown/zero-confidence alignments. A complete primary
  query alignment is required before treating its interval as the full element.
- Tail-head support counts actual query-aligned bases on each side and requires
  a match block across the junction itself. Deletions/reference skips cannot
  substitute for a bridge; indels away from the junction remain permitted.
  Hard clipping no longer counts as observed genomic flank sequence.

## Existing results and scope

Existing formatter outputs, BAM caches and old biological conclusions are not
automatically revalidated. Use a fresh output directory for a clean comparison
of the updated pipeline. Old partial-circle directories without a current
completion marker are rerun; nested legacy outputs are not promoted to verified
results. Keep chunk reference directories and metadata when cleaning large
FASTQ/BAM intermediates.

The F3 change fixes the reproduced truncation/flag handling, not the entire
quality assessment: it does not establish genomic locus identity or add a
sequence-identity quality model. Full-query validation is intentionally
conservative and can reject clipped true candidates. F1, F2, F5 and the other
unaddressed items in `audit_F1_F10.md` remain open. F4's CIGAR handling is covered
by synthetic tests; biological sensitivity/specificity has not been measured.

## Tests

With Python 3.11+ and `pandas`, `biopython`, `pysam` installed:

```bash
python -m unittest discover -s tests -v
bash -n submit_all.sh
bash -n submit_one_batch.sh
```

The regression suite uses real subprocesses for pipe exits, stderr pressure and
timeouts, and captures real shell-generated submission commands for interleaved
reference-resolution tests. It mocks SLURM, external aligners and BAM records.
Local Windows validation used real pandas but an explicit pysam import stub;
no real BAM I/O, sequence alignment, cluster job or production-data benchmark
was performed.
