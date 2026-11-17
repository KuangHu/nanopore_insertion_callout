# Manifest

## Repository layout

| Path | Contents |
|---|---|
| `cycle/` | Python package: download, reference resolution/alignment (`preprocess/`), Sniffles2 and tldr runners, circle detection, IS formatter, partial circles, ORF annotator, guide finder, system clustering, novelty annotator, IS110 filter, nesting detector, visualizer; `site_config.py` (site paths) |
| `scripts/` | Entry points: `run_pipeline.py` (stages download -> format), per-module runners `run_*.py`, IS110 collection scripts `is110_*.py/.sh`, status/QA tools, HuggingFace upload |
| `submit_one_batch.sh` | Current batch submitter (4 SLURM phases per batch, N chunks) |
| `submit_all.sh` | Previous multi-batch submitter (kept for reference) |
| `bin/cycle-sbatch` | Submits standalone SLURM scripts with the site's scheduler settings |
| `config/` | Site settings: `igi.env` (IGI defaults), `lawrencium.env.example`, `load_site.sh` |
| `sideproject/` | One-off analyses: active IS110 collection, IS110 split finder (GTDB), IS621 search, ECOR samples |
| `data/` | Inputs: `sra_metadata*.tsv` (SRA search results), `batches/batch_NNN.tsv` (sample lists per batch) |
| `envs/` | Environment inventory (exact and portable specs, tool versions) |
| `tests/` | Unit and integration tests (`pytest tests`) |
| `notebooks/` | Exploratory notebook (junction reads) |
| `docs/` | Pipeline parts, audit F1-F10, reliability fixes, Lawrencium setup, `migration/` |
| `PIPELINE_*.md`, `DATA_STRUCTURE.md` | Methods, notes, summaries, output layout (IGI paths are historical) |

## Outputs under `IS_cycle/` and what produces them

| Output | Produced by |
|---|---|
| `batch_NNN/chunks/chunk_N/sniffles_output/`, `.../circle_output/`, `batch_NNN/is_formatter_output/*_is_records.json` | `submit_one_batch.sh` Phase 1 -> `scripts/run_pipeline.py --steps download resolve index align sniffles circle format` |
| `batch_NNN/reference_genomes/` | `run_pipeline.py` resolve + index steps (NCBI datasets, `minimap2 -d`) |
| `batch_NNN/is_formatter_output/*_annotated.json`, `*_guide.json` | Phase 2 -> `scripts/run_orf_annotator.py`, `scripts/run_guide_finder.py` |
| `system_clustering_batch_NNN/`, `novelty_batch_NNN/` | Phase 3 -> `scripts/run_system_clustering.py`, `scripts/run_novelty_annotator.py` |
| `batch_NNN/partial_circle_output/` | Phase 4 -> `scripts/run_partial_circle_batch.py` |
| `batch_000/is_reference/` | `run_pipeline.py` is_ref step (ISfinder) |
| `batch_000…004/` (legacy layout, batch-level `sniffles_output/`) | `submit_all.sh` and earlier per-step runs (Feb-Mar 2026) |
| `is110_all_006_026/`, `_006_052/`, `_006_070/` | `scripts/is110_run_all.sh`, `is110_collect_all.py`, `is110_step345.py`, `is110_build_json.py` |
| `active_is110/`, `active_is110_000_094/` | `sideproject/active_is110/collect_active_is110.py`, `collect_active_is110_000_094.py` |
| `is110_split_gtdb/` | `sideproject/run_split_finder.sh`, `run_is110_blastn.sh` |
| `is621_search/` | `sideproject/is621_search/*.sh` |
| `combined_006_010/`, `visualizations_guide12_30*/` | `scripts/combine_and_visualize.py`, `scripts/visualize_guide12_30.py` |
| `is110_all_000_084/`, `sniffles_is110_006_078/`, `combined_is110_006_078/`, `all_insertions_*`, `assembly_catalog_v3/`, `insertions_with_th_junction_000_094/`, `insertions_th_expanded_000_094/`, `cycle_calls_006_078/`, delivery zips | sibling project `cycle_rna_analysis` (`scripts/` there; packaged as its own component) |
| `is110_gtdb_output/`, `is110_ecoli_blast/`, `is110_ncbi/`, `is110_top20/`, `is110_combined_006_020*/`, `is110_006_026/`, `is110_all_batch_000*/`, `is110_circular_batch_000/` | early side analyses; producing script not found in this repository or the sibling projects |

Raw reads, BAMs, minimap2 indices and the legacy assembly `reads.fq` files are not kept; see
`docs/migration/READS_FQ_REPRODUCIBILITY.md` for how removed intermediates can be regenerated.
