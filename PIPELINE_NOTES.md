# IS110 Circular Intermediate Detection Pipeline — Methods & Versions

End-to-end pipeline for detecting IS110-family insertion sequences and their circular intermediates from long-read sequencing data (Oxford Nanopore / PacBio).

---

## Overview

```
SRA FASTQ → minimap2 → Sniffles2 → CircleFinder (concatemer-bait)
            → Local assembly (miniasm+minipolish) → Prodigal ORFs
            → HMMER (DEDD + Tnp20) → IS110 transposase
            → Guide finder (flanking ↔ noncoding)
            → CReSIL + NanoCircle (partial circles)
```

All steps run under conda environment `opfi` unless otherwise noted. Cycle calling (CReSIL/NanoCircle) runs under `cresil`.

---

## 1. Download — Public sequencing reads

**Tool:** kingfisher (v0.4.1)
**Source:** ENA/SRA via `ena-ftp`, `aws-http`, `prefetch` (in order)
**Output:** `.fastq.gz` per accession

Metadata retrieved via NCBI Entrez (BioProject, organism, instrument, total bases) and stored in per-chunk TSVs.

---

## 2. Reference Genome Resolution

**Tool:** NCBI Datasets CLI (v18.13.0)
**Source:** RefSeq / GenBank assemblies matched to sample organism
**Output:** `.fna` (assembly) + `.fna.fai` index

For each sample, the best-matching complete reference genome is chosen (organism + strain match preferred).

---

## 3. Alignment — Reads to reference

**Tool:** minimap2 (v2.30-r1287) + samtools (v1.17)
**Preset:** `map-ont` for Nanopore, `map-pb` for PacBio
**Options:** `-a --MD -t N` piped through `samtools sort` → sorted indexed BAM
**Timeout:** 2 hours per sample to prevent runaway alignments

---

## 4. Insertion Calling

**Tool:** Sniffles2 (v2.7.2)
**Options:** `--allow-overwrite --minsupport 3 --minsvlen 100`
**Output:** VCF + `.table.txt` with per-insertion UUID, chromosome, position, consensus sequence, supporting reads

The `.table.txt` is the primary input for downstream steps. Each insertion gets a unique UUID and a consensus sequence built from supporting reads.

---

## 5. Circle Detection — Concatemer Bait (Full Circle)

**Tool:** Custom module `cycle.circle_detect.CircleFinder` (using pysam v0.23.3, minimap2)
**Method:**
1. Build `[IS][IS]` concatemer bait from Sniffles consensus
2. Re-align reads to bait with minimap2 (`-x map-ont`)
3. Classify junction-spanning reads:
   - **tail-head** (`n_tail_head_reads`): both halves of read align to inside of bait, spanning the IS-IS junction → **direct circular intermediate evidence**
   - **genome-head** / **tail-genome**: half in bait, half in flanking → insertion site reads
4. Minimum overlap: 100 bp on each side of junction (configurable)

**Output:** `<sample>_circle_summary.tsv` + `<sample>_circle_reads.tsv`

---

## 6. Local Assembly & Flanking Extraction

**Tools:**
- minimap2 (v2.30, overlap mode for read self-alignment)
- miniasm (v0.3-r179)
- minipolish (v0.2.0)
- mafft / exonerate (consensus polishing)

**Module:** `cycle.is_formatter.ISFormatter`
**Output:** `<sample>_is_records.json`

For each Sniffles IS with TH reads ≥ 1:
- Extract reads spanning the IS region
- miniasm + minipolish to produce a cleaner IS consensus + flanking
- Falls back to longest read if assembly fails
- 80 bp flanking on each side stored separately

---

## 7. ORF Annotation

**Tool:** Prodigal (v2.6.3, Feb 2016)
**Mode:** `-p meta` (metagenomic, no training required)
**Module:** `cycle.orf_annotator.ORFAnnotator`
**Output:** `<sample>_is_records_annotated.json` (records + ORFs with protein sequences)

---

## 8. IS110 Transposase Identification

**Tool:** HMMER (v3.4, hmmsearch)
**HMM databases:**
- `DEDD_Tnp_IS110` (DEDD recombinase domain, HMMER3/f v3.3, Nov 2019)
- `Transposase_20` (Tnp20, HMMER3/f v3.3, Nov 2019)
- Source: ISfinder-derived HMMs at `/home/kuangh/scripts/IS110/hmm/`
- E-value threshold: **1e-5**

**Module:** `cycle.is110_filter.IS110Filter`
**Logic:**
- **Case 1**: single ORF carries both DEDD and Tnp20 hits → IS110 transposase
- **Case 2**: two adjacent ORFs (gap ≤ 300 bp), one with DEDD, one with Tnp20 → split IS110
- Records with only one domain are rejected

**Output:** `is110_ids.json` — UUIDs of IS elements carrying IS110-family transposase

---

## 9. Guide Finder — Flanking ↔ Internal Noncoding Alignment

**Tool:** Custom `cycle.guide_finder.ShortAlignmentFinder` (pure-Python seed-and-extend)
**Tiered rule:**
- 9 bp exact match (0 mismatches), OR
- 11 bp with ≤ 1 mismatch (12 bp minimum for mismatch allowance)
- Forward + reverse-complement search

**Output:** `<sample>_is_records_guide.json` with per-record `guide_hits` and `guide_summary` (best length, total hits).

---

## 10. Partial Circle Detection — Split-Read Back-Jump

**Tool:** Custom `cycle.partial_circle.PartialCircleDetector` (pysam, minimap2)
**Method:**
1. Build single-copy IS reference from guide JSONs
2. Map reads with minimap2 (`-x map-ont`), name-sorted BAM
3. For each split read with 2+ alignment segments on same IS:
   - Detect **back-jump** (segment B's reference start < segment A's reference end)
   - The partial circle = `[seg_b.ref_start, seg_a.ref_end]`
4. Single-linkage clustering of breakpoints (≤ 5 bp tolerance), require ≥ 2 supporting reads

**Output:** `<sample>_partial_circle_summary.json` + `.tsv`

---

## 11. CReSIL + NanoCircle — Genome-Wide Circular DNA

For broader (non-IS) circular DNA detection:

**Tools** (conda env `cresil`):
- **CReSIL** (Identification of extrachromosomal circular DNA from long-read sequencing; from `cresil_nanocircle/CReSIL`)
- **NanoCircle** (Simple + Chimeric cycle calling; from `cycle_rna_analysis/NanoCircle`)
- bedtools (v2.31.1) for coverage BED

**Pipeline:** `pilot_one_sample.sh` per sample
- Re-download FASTQ → align to reference → CReSIL + NanoCircle in parallel → per-sample `<sample>_cycle_junctions.tsv`

---

## 12. System Clustering (Optional)

**Tool:** MMseqs2 + Louvain community detection
**Module:** `cycle.system_clustering.SystemClusterer`
**Purpose:** Group transposons sharing protein families across samples

---

## Outputs Summary

| Stage | File | Per |
|-------|------|-----|
| Sniffles table | `<sample>.table.txt` | sample |
| Circle summary | `<sample>_circle_summary.tsv` | sample |
| Raw IS records | `<sample>_is_records.json` | sample |
| Annotated (with ORFs) | `<sample>_is_records_annotated.json` | sample |
| Guide hits | `<sample>_is_records_guide.json` | sample |
| Partial circles | `<sample>_partial_circle_summary.json` | sample |
| IS110 IDs (cross-batch) | `is110_ids.json` | study |
| Cycle junctions (CReSIL/NC) | `<sample>_cycle_junctions.tsv` | sample |

---

## Software Versions (snapshot)

| Tool | Version | Purpose |
|------|---------|---------|
| kingfisher | 0.4.1 | SRA download |
| NCBI Datasets CLI | 18.13.0 | Reference download |
| minimap2 | 2.30-r1287 | Long-read alignment |
| samtools | 1.17 | BAM manipulation |
| Sniffles2 | 2.7.2 | SV / insertion calling |
| pysam | 0.23.3 | BAM parsing |
| miniasm | 0.3-r179 | Assembly |
| minipolish | 0.2.0 | Consensus polishing |
| Prodigal | 2.6.3 | ORF prediction |
| HMMER | 3.4 | HMM search |
| pyhmmer | 0.11.2 | Python HMMER bindings |
| BLAST+ (blastp) | 2.17.0+ | Protein similarity |
| bedtools | 2.31.1 | Coverage / interval ops |
| MMseqs2 | (latest in `opfi` env) | Protein clustering |
| CReSIL | (from `cresil_nanocircle/CReSIL`) | Long-read eccDNA |
| NanoCircle | (from `cycle_rna_analysis/NanoCircle`) | Long-read circles |
| Prokka | not used — replaced by Prodigal+HMMER | — |

## Database Versions

| Database | Version | Purpose |
|----------|---------|---------|
| `DEDD_Tnp_IS110.hmm` | HMMER3/f v3.3 (Nov 2019) | IS110 DEDD recombinase |
| `Transposase_20.hmm` | HMMER3/f v3.3 (Nov 2019) | IS110 Tnp20 transposase |
| RefSeq / GenBank | retrieved on-demand via Datasets CLI | reference genomes |
| ISfinder | (consulted for IS family context, not directly queried in cycle pipeline) | — |

---

## Filtering Choices

- Minimum IS consensus length: 500 bp
- Minimum IS consensus entropy: 1.0 bit (low-complexity filter)
- Sniffles minimum support: 3 reads
- Circle junction minimum overlap: 100 bp / side
- Guide hit rule: 9/9 exact OR 11/12 (≤1 mm)
- Partial circle: ≥2 supporting reads, ≥50 bp, ≤80% of IS length
- IS110 transposase HMM E-value: 1e-5
- ORF–ORF adjacency for split DEDD/Tnp20: ≤300 bp
