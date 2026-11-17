# Deleted legacy `reads.fq` files: provenance and reproducibility

On 2026-11-17, during the pre-migration cleanup, 1,787,475 files (2.94 TB) named
`batch_00[0-4]/is_formatter_output/<run>/assembly/<uuid12>/reads.fq` were deleted from the IGI
`IS_cycle` tree (3,091 runs). The deletion manifest (path and byte size of every file) and per-file
deletion log are kept in the migration bundle (`provenance/cleanup_20261117/`).

This note documents what those files were, where their inputs still are, and how they can be
recreated. It also records a pilot that recreated a sample of them.

## What the files were

`reads.fq` was an intermediate written by the IS formatter (`cycle/is_formatter/formatter.py`).
For one IS element (Sniffles UUID), it held every primary, mapped read from the run's
circle-detection BAM whose reference was one of that UUID's bait contigs (`{uuid}__...`),
de-duplicated by read name, with original qualities. It was written only when an element had at
least 3 such reads (`DEFAULT_MIN_READS_FOR_ASSEMBLY`) and was the input to the miniasm + minipolish
local assembly. The formatter's results (`<run>_is_records.json` / `.tsv`, plus the later
`_annotated.json` and `_guide.json`) were kept and are in the results archive.

Path mapping: `<uuid12>` is the first 12 characters of the Sniffles UUID, which is `is_id` in the
kept records and the `UUID` column of the kept per-run Sniffles table.

## Provenance (from the IGI SLURM logs)

| Runs | Bait FASTA and circle BAM built by | Date |
|---:|---|---|
| 3,082 | batch "step3" jobs 1029690–1029694 (one per batch 000–004) | 2026-02-20 |
| 8 | batch_000 jobs 1025741 / 1027052 | 2026-02-17 / 18 |
| 1 | no bait-creation record found (batch_002) | unknown |

The 2026-02-20 jobs built per-run baits with both tail-head concatemers and insertion-context
sequences (`[upstream 500 bp][IS][downstream 500 bp]` from the run's reference genome), mapped the
run's FASTQ with `minimap2 -a --secondary=no -x map-ont`, and ran the formatter on every element.
Later formatter passes (2026-02-21 and 2026-03-03) reused the same baits and BAMs and only
rewrote a subset of the files (elements with tail-head reads), from the same reads.

Those jobs ran working-tree code that was never committed verbatim. The closest commit is
`03072ec` (2026-02-21): same per-run circle detection and bait design; it adds the
`circle_evidence.example_th_read` field and a tail-head filter that `--no-require-th` turns off.

## Inputs and where they remain

| Input | Location | Status |
|---|---|---|
| Raw reads | public SRA/ENA runs; run accession = folder name; sample lists in `data/batches/batch_00[0-4].tsv` | 100 of 100 randomly sampled runs present in NCBI SRA on 2026-11-17 (97 with ENA FASTQ, 3 only via the SRA route) |
| Per-run Sniffles table (UUIDs, consensus) | `IS_cycle/batch_00X/sniffles_output/<run>/<run>.table.txt` | kept; in the IS_cycle results archive |
| Reference genome per run | `IS_cycle/batch_00X/metadata_for_sniffles.tsv` (run -> FASTA) and `reference_genomes/<acc>/*.fna` | kept; FASTA in the deduplicated reference archive |
| Code | commit `03072ec` of this repository | in git |
| Software | `opfi` env (minimap2 2.30, pysam, miniasm 0.3, minipolish 0.2.0, kingfisher 0.4.1) + samtools 1.17 | `envs/` |
| Parameters | `cycle/circle_detect/config.py` and `cycle/is_formatter/config.py` at `03072ec`: min junction overlap 100 bp, min consensus 200 bp, entropy >= 1.7, insertion-context flank 500 bp, formatter min reads 3, flank 80 bp, assembly timeout 120 s; formatter run with `--no-require-th` | thread count does not change which reads map |
| Download route | kingfisher `-m ena-ftp aws-http prefetch -f fastq.gz` (pipeline default) | reads.fq stores read names and qualities, so use the same route |

## Recreating one run

```bash
git -C <repo> archive 03072ec | tar -x -C code_03072ec
R=work/<run>; mkdir -p $R/sniffles_output/<run> $R/alignments
cp IS_cycle/<batch>/sniffles_output/<run>/<run>.table.txt $R/sniffles_output/<run>/
ref=$(awk -F'\t' '$1=="<run>"{print $3}' IS_cycle/<batch>/metadata_for_sniffles.tsv)   # remap to the new IS_cycle root
printf "sample_id\treference\tstatus\n<run>\t$ref\tok\n" > $R/alignments/alignment_status.tsv
(head -1 data/batches/<batch>.tsv; grep -P "^<run>\t" data/batches/<batch>.tsv) > $R/meta.tsv
cd code_03072ec && python scripts/run_pipeline.py --metadata $R/meta.tsv --outdir $R \
   --fastq-dir $R/sra_downloads --align-dir $R/alignments \
   --steps download circle format --no-require-th --parallel 1 --threads 8
# then delete $R/sra_downloads and the circle BAMs (stream-and-delete)
```

The pilot script that does exactly this is in the migration bundle (`tools/repro_pilot_readsfq_v2.sh`).

## Pilot result (2026-11-17)

| Run | Batch | Species | Deleted `reads.fq` | Recreated, same path and byte size |
|---|---|---|---:|---:|
| SRR13575581 | 003 | *S. pyogenes* | 13 | 13 |
| SRR25466773 | 000 | *E. coli* | 7 | 7 |
| SRR31647025 | 000 | *S. aureus* | 7 | 7 |
| SRR31235479 | 004 | *E. coli* | 5 | 5 |

All 32 files were recreated at the same `assembly/<uuid12>/` paths with identical byte sizes. The
recreated `is_records.json` matched the kept file in every field except the added
`example_th_read`. A first pilot with the 2026-03-03 code (`97039a6`) did **not** reproduce them
(it built tail-head-only baits and filtered elements), which is why the code version above matters.

Limits: 4 of 3,091 runs (32 of 1,787,475 files), all small; batches 001 and 002 were not sampled;
the original file contents are gone, so byte size (not checksum) is the comparison; 9 runs had baits
from other jobs. Raw-read availability can change at SRA/ENA.

## Results archive versus reproducible rerun package

| | Results archive (transferred) | Reproducible rerun package |
|---|---|---|
| Purpose | continue analysis on Lawrencium without recomputing | regenerate removed intermediates or rerun steps |
| Contains | IS_cycle outputs as they stand after cleanup: per-run Sniffles tables and VCFs, circle summaries and reads tables, IS records JSON, partial-circle calls, novelty and clustering results, catalogs, delivery zips, one copy of each reference FASTA | this repository at a stated commit, `envs/` (exact and portable), `config/`, sample lists (`data/batches`), assets (HMMs, Pfam), and the run/commit/parameter mapping above |
| Does not contain | raw reads, BAMs, minimap2 `.mmi` indices, legacy `reads.fq` | any data that must be downloaded (SRA reads) or recomputed |
| Recreating | n/a | `.mmi`: `minimap2 -d` from the FASTA; `reads.fq`: procedure above; BAMs and raw reads: download + align |
