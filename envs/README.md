# Software environments

Two conda environments are used. Files here were exported on 2026-11-17 from the IGI cluster.

| File | What it is | Use it to |
|---|---|---|
| `opfi.explicit.txt`, `cresil.explicit.txt` | **Exact inventory**: every package with build string, channel URL and MD5 (`conda list --explicit --md5`), linux-64 only | rebuild the identical environment: `conda create -n opfi --file envs/opfi.explicit.txt` |
| `opfi.environment.yml`, `cresil.environment.yml` | All packages with versions, no build strings, incl. pip section (`conda env export --no-builds`) | portable re-solve on another linux-64 cluster: `mamba env create -n opfi -f envs/opfi.environment.yml` |
| `opfi.from-history.yml`, `cresil.from-history.yml` | Only the packages that were explicitly requested (`--from-history`) | minimal re-solve if the full spec no longer resolves |
| `opfi.pip-freeze.txt`, `cresil.pip-freeze.txt` | pip view of each environment (`pip list --format=freeze`) | check pip-only packages |
| `tool_versions.tsv` | Versions of the external programs the code calls, and where each came from | verify a rebuilt environment |

- `opfi` (Python 3.12.3): the batch pipeline, IS110 analyses and side projects.
- `cresil` (Python 3.6.7): the CReSIL/NanoCircle broad-cycle pipeline in `cycle_rna_analysis`
  (needs the CReSIL `rsplit("_", 2)` patch described in that project's notes).

## Not in the conda environments

- **samtools**: on IGI the pipeline used the system install `/shared/software/bin/samtools`
  (1.17, htslib 1.17); `opfi` has no samtools package. On a new cluster install it into the env:
  `mamba install -n opfi -c bioconda samtools=1.17`.
- **Reference assets** (not software, but required): IS110 HMMs `DEDD.hmm` / `Tnp20.hmm`
  (`CYCLE_HMM_DIR`), `Pfam_transposase.hmm`, and the shared Pfam-A 2020 release used by
  eggNOG-mapper; see `docs/migration/README.md` (assets) for checksums and sources.

The root-level `environment.yml` and `requirements.txt` are the February 2026 exports and are kept
for history; the files in this directory supersede them.
