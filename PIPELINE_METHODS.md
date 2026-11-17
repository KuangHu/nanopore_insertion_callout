# Methods

## Sequencing data acquisition and reference selection

Public long-read sequencing data (Oxford Nanopore and PacBio) were retrieved from the European Nucleotide Archive (ENA) and the Sequence Read Archive (SRA) using `kingfisher` (v0.4.1), attempting `ena-ftp`, `aws-http`, and NCBI `prefetch` in order. For each run, sample metadata (BioProject, organism, instrument, total bases) was obtained from the SRA metadata service. Reference genomes for each run were downloaded with the NCBI Datasets command-line tool (v18.13.0), preferring complete RefSeq assemblies matched to the reported organism and strain.

## Read alignment

Long reads were aligned to the matched reference genome with `minimap2` (v2.30-r1287) using the `map-ont` preset for Oxford Nanopore reads and the `map-pb` preset for PacBio reads, with the `-a --MD` flags so that supplementary alignment and per-base match information were retained. Alignments were piped directly into `samtools sort` (v1.17) to produce coordinate-sorted, indexed BAM files. A two-hour per-sample wall-clock limit was enforced on the alignment step to prevent runaway processes on pathological inputs.

## Structural variant calling and insertion consensus

Insertion-type structural variants were called from the sorted BAMs with Sniffles2 (v2.7.2) using `--minsupport 3 --minsvlen 100` and `--allow-overwrite`. The output `.table.txt` file was used downstream: each insertion received a unique UUID together with its chromosome, position, the consensus insertion sequence built from supporting reads, and the supporting-read count.

## Detection of full circular intermediates via concatemer-bait junction reads

To detect circular intermediates produced by IS-element excision, a `[IS][IS]` concatemer "bait" sequence was constructed for each Sniffles insertion by duplicating the consensus end to end. All reads were re-aligned to this bait with `minimap2` (`-x map-ont`), and reads spanning the head-to-tail junction were identified using a custom Python module (`cycle.circle_detect`) built on `pysam` (v0.23.3). A read was classified as a tail-head junction read when both halves of the alignment lay within the IS body and together spanned the central IS-IS junction with at least 100 bp on each side, providing direct evidence of a circular intermediate. Reads with one segment inside the bait and one segment in the genomic flanking were classified as genome-head or tail-genome reads, marking the integration site. Per-insertion summaries (`*_circle_summary.tsv`) were used to filter insertions with at least one tail-head read for downstream assembly.

## Local assembly of IS body and flanking sequences

For each insertion supported by at least one tail-head read, the supporting reads were extracted and assembled de novo with `miniasm` (v0.3-r179) using overlap detection by `minimap2`, followed by polishing with `minipolish` (v0.2.0). This yielded a clean per-insertion consensus together with up to 80 bp of upstream and downstream genomic flanking on each side. When local assembly failed, the longest supporting read was used as the consensus. The output was stored as one structured JSON record per insertion (`*_is_records.json`) containing the IS body, flanking sequences, and circular-intermediate read counts.

## Open reading frame annotation

Open reading frames (ORFs) within each IS consensus were predicted with Prodigal (v2.6.3) run in metagenomic mode (`-p meta`), which performs gene prediction without requiring a training set. Translated protein sequences and ORF coordinates were appended to each IS record (`*_is_records_annotated.json`).

## Identification of IS110-family transposases

IS110-family transposases were identified by searching all predicted proteins with HMMER (v3.4, `hmmsearch`) against two profile HMMs: a DEDD-recombinase profile (`DEDD_Tnp_IS110`) and a transposase profile (`Transposase_20`), both built and distributed in HMMER3/f format (v3.3, November 2019), at an E-value threshold of 1e-5. An insertion was classified as carrying an IS110-family transposase when either (i) a single ORF was hit by both HMMs or (ii) two ORFs separated by no more than 300 bp were hit, one by each HMM (split DEDD-Tnp20 architectures). Insertions matching only one of the two domains were rejected.

## Guide-region detection

To detect candidate RNA-guide sequences encoded within each IS, a tiered short-alignment search was performed between the genomic flanking regions and the noncoding internal regions of the IS body using a custom Python implementation (`cycle.guide_finder`). Both forward and reverse-complement orientations were searched. A hit was retained when it was either a perfect 9 bp match (zero mismatches) or an at-least-11 bp match with no more than one mismatch (mismatches were only allowed in alignments of 12 bp or longer). Per-insertion guide hits and a summary (best length, total hits) were appended to each record (`*_is_records_guide.json`).

## Detection of partial circular intermediates

Partial-circle intermediates — sub-regions of an IS that excise to form smaller circles — were detected from split-read back-jump signals (`cycle.partial_circle`). For each sample, a single-copy reference of all called IS bodies was built and reads were re-aligned with `minimap2` (`-x map-ont`); name-sorted BAMs were parsed with `pysam`. For each read with two or more alignment segments on the same IS, segments were ordered by query position and consecutive pairs were inspected for a back-jump (segment B beginning before segment A ends on the IS reference). Each back-jump defines a candidate partial circle spanning the interval `[start_B, end_A]` on the IS. Candidates were filtered to circles of at least 50 bp and at most 80 % of the IS length, with at least 10 bp of overlap on each side of the junction. Breakpoints from independent reads were clustered by single linkage with a 5 bp tolerance, and clusters supported by at least two reads were reported (`*_partial_circle_summary.json`).

## Genome-wide circular DNA detection

For comparison with non-IS circular DNA, each sample was additionally processed with CReSIL (Identification of extrachromosomal circular DNA from long-read sequencing) and NanoCircle (simple + chimeric circle calling), as packaged in `/home/kuangh/tools/cycle_rna_analysis/cresil_nanocircle`. Coverage tracks for downstream normalisation were generated with `bedtools` (v2.31.1) using `genomecov -bg` and merged within 1 kb. Per-cycle junction sequences were extracted from each call.

## Sequencing-depth normalisation

To enable cross-sample comparison of circular-intermediate abundance, a per-sample on-target sequencing depth was computed as the ratio of the SRA-reported total sequenced bases to the total length of regions covered by at least one read in the coverage BED. Per-sample circle rates were then defined as the number of called cycles (or supporting junction reads) divided by this on-target depth.

## Software and database versions

The pipeline was run under a conda environment named `opfi` for the main steps and `cresil` for the genome-wide circular-DNA tools. Software versions were: `kingfisher` 0.4.1; NCBI Datasets CLI 18.13.0; `minimap2` 2.30-r1287; `samtools` 1.17; `pysam` 0.23.3; Sniffles2 2.7.2; `miniasm` 0.3-r179; `minipolish` 0.2.0; Prodigal 2.6.3; HMMER 3.4; `pyhmmer` 0.11.2; NCBI BLAST+ `blastp` 2.17.0+; `bedtools` 2.31.1. Profile HMM databases used were the IS110 DEDD-recombinase profile (`DEDD_Tnp_IS110`) and the IS110 transposase profile (`Transposase_20`), both in HMMER3/f format version 3.3 (November 2019), originally derived from ISfinder family models.
