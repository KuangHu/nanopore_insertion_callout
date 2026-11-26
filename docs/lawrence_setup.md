# Running Cycle on Lawrencium

The pipeline was developed on the IGI cluster (`igi.biotite`). All cluster-specific locations and
scheduler settings now live in one site file, so moving to Lawrencium needs no code edits.
See `docs/migration/README.md` for the full migration procedure and the transfer packages.

## 1. Code and environment

```bash
git clone https://github.com/KuangHu/nanopore_insertion_callout.git ~/tools/Cycle
cd ~/tools/Cycle && git checkout migration/lawrencium   # until merged into main

# Recreate the conda environment (details and the exact pinned inventory: envs/README.md)
mamba env create -n opfi -f envs/opfi.environment.yml
# samtools is NOT in opfi on IGI (it came from /shared/software, v1.17); add it explicitly:
mamba install -n opfi -c bioconda samtools=1.17
```

## 2. Site settings

```bash
cp config/lawrencium.env.example config/site.env    # git-ignored
$EDITOR config/site.env                              # replace every CHANGE_ME
```

| Variable | Purpose | IGI value (config/igi.env) |
|---|---|---|
| `CYCLE_DATA_ROOT` | IS_cycle output tree | `/groups/rubin/projects/kuang/out/IS_cycle` |
| `CYCLE_LOG_DIR` | SLURM stdout/stderr | `$HOME/logs` |
| `CYCLE_CONDA_ENV`, `CYCLE_ENV_PREFIX` | conda env name and prefix. Every job finds the environment by the prefix: the batch submitters prepend `$CYCLE_ENV_PREFIX/bin` to PATH and the standalone scripts use `conda run -p "$CYCLE_ENV_PREFIX"`; the name is informational | `opfi`, `/home/kuangh/miniconda3/envs/opfi` |
| `CYCLE_HMM_DIR` | `DEDD.hmm`, `Tnp20.hmm` for the IS110 filter | `/home/kuangh/scripts/IS110/hmm` |
| `CYCLE_ISFINDER_FASTA` | ISfinder reference for novelty annotation | `$CYCLE_DATA_ROOT/batch_000/is_reference/ISfinder_raw.fna` |
| `CYCLE_BATCH_METADATA_DIR` | `batch_NNN.tsv` sample lists | `<repo>/data/batches` |
| `CYCLE_GTDB_GENOMES_DIR`, `CYCLE_GTDB_BLASTDB`, `CYCLE_ECOR_DIR` | side-project inputs | shared GTDB r220, (missing) BLAST db, `out/IS110/ECOR_batch` |
| `CYCLE_SLURM_ACCOUNT` | `--account` (required on Lawrencium) | empty |
| `CYCLE_SLURM_PARTITION`, `CYCLE_SLURM_QOS` | `--partition`, `--qos` | `standard`, `standard` |
| `CYCLE_SLURM_CPUS`, `CYCLE_SLURM_MEM` | per-node job size; thread counts scale from CPUs | `48`, `192G` |
| `CYCLE_SLURM_AVX2_EXCLUDE` | nodes to avoid for mmseqs2 (Phase 3 and `--needs-avx2` jobs) | 11 `node-48-256g-*` nodes |
| `CYCLE_SLURM_EXTRA` | any other sbatch flags, e.g. `--constraint=...` | empty |

Precedence: exported environment variables, then `config/site.env` (or the file named by
`CYCLE_SITE_ENV`), then `config/igi.env`. Python code reads the same variables through
`cycle/site_config.py`.

## 3. Submitting work

```bash
bash submit_one_batch.sh 95 --nodes 3                        # batch pipeline (all 4 phases)
bin/cycle-sbatch scripts/is110_run_all.sh                    # standalone SLURM scripts
bin/cycle-sbatch --needs-avx2 sideproject/active_is110/run_collect.sh
CYCLE_DRY_RUN=1 bin/cycle-sbatch scripts/is110_build_json.sh # print the sbatch command only
```

Standalone scripts no longer carry partition, qos, account or log paths in their `#SBATCH` headers;
`bin/cycle-sbatch` adds them from the site file. Extra sbatch options after the script name override
the header (for example `--cpus-per-task=32`), and thread counts inside the scripts follow
`SLURM_CPUS_PER_TASK`.

## 4. Data directory structure

```
IS_cycle/
  batch_000/
    sra_downloads/       # FASTQs (deleted after processing)
    reference_genomes/   # legacy submit_all.sh layout
    is_formatter_output/ # final IS records + JSONs (kept)
    chunks/
      chunk_N/
        reference_genomes/ # chunk-local references, indices, status (keep)
        sra_downloads/     # removable after all phases are verified
        alignments/        # retain status metadata; BAMs are removable
  system_clustering_batch_000/
  novelty_batch_000/
```

`submit_one_batch.sh` isolates references per chunk so concurrent jobs cannot
overwrite each other's state or partially downloaded files. Do not delete the
whole `chunks/` tree when retaining reproducible reference mappings. See
[reliability fixes](reliability_fixes.md) for completion markers and rerun behavior.

The batches transferred from IGI (000–094) predate per-chunk references and keep
`batch_NNN/reference_genomes/`; their minimap2 `.mmi` indices were removed before
transfer and are rebuilt by the `index` step (or `minimap2 -d`) when a batch is re-aligned.

## 5. Quick test

```bash
source config/load_site.sh
head -4 data/batches/batch_006.tsv > /tmp/batch_test.tsv
python scripts/run_pipeline.py --metadata /tmp/batch_test.tsv \
  --outdir "$CYCLE_DATA_ROOT/batch_test" \
  --steps download resolve index align sniffles circle format --parallel 1 --threads 8
```
