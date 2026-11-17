# IS110 Detection and Circular Intermediate Pipeline (no guide step)

## Pipeline

Long-read sequencing data (Oxford Nanopore or PacBio) were retrieved from the
European Nucleotide Archive (ENA) and the NCBI Sequence Read Archive (SRA)
using `kingfisher` (v0.4.1), which attempted `ena-ftp`, `aws-http`, and NCBI
`prefetch` in order until a run downloaded successfully; per-run metadata
(BioProject, organism, instrument, total bases) was obtained from the SRA
metadata service. For each run, a matched closed reference genome was
downloaded from NCBI RefSeq / GenBank using the NCBI Datasets command-line
tool (`datasets`, v18.13.0), preferring complete assemblies matched to the
reported organism and strain.

Long reads were aligned to the matched reference with `minimap2` (v2.30-r1287,
`map-ont` or `map-pb` preset) and the alignments piped through `samtools sort`
(v1.17) to produce a coordinate-sorted, indexed BAM. Non-reference insertions
were called from the BAM with Sniffles2 (v2.7.2, `--minsupport 3 --minsvlen 100`);
each insertion received a UUID together with its genomic coordinates, supporting
reads, and a consensus insertion sequence built from those reads.

Each Sniffles consensus was then locally assembled with `miniasm` (v0.3-r179)
and polished with `minipolish` (v0.2.0), and up to 80 bp of upstream and
downstream genomic flanking was recovered around the insertion site. Open
reading frames inside the assembled insertion were predicted with Prodigal
(v2.6.3, `-p meta`), yielding translated protein sequences and coordinates.

Predicted proteins were searched with HMMER `hmmsearch` (v3.4, E-value 1e-5)
against two profile HMMs in HMMER3/f format (v3.3, November 2019): the
IS110-family DEDD recombinase profile (`DEDD_Tnp_IS110.hmm`) and the
IS110-family transposase profile (`Transposase_20.hmm`). Both HMMs were
originally derived from ISfinder family models and are archived locally at
`/home/kuangh/scripts/IS110/hmm/`. An insertion was classified as an
IS110-family element when either (i) a single ORF was hit by both HMMs, or
(ii) two ORFs separated by no more than 300 bp were hit by the two HMMs
respectively, corresponding to the split DEDD–Tnp20 architecture typical of
IS110-family transposases. Insertions matching only one of the two domains
were rejected.

To detect full circular intermediates, a `[IS][IS]` concatemer "bait" was
constructed for each IS by duplicating the Sniffles consensus end to end. All
reads were re-aligned to the bait with `minimap2` (`-x map-ont`), and reads
spanning the head-to-tail junction with at least 100 bp of alignment on each
side were counted as tail-head junction reads. These reads provide direct
evidence of the full circular intermediate produced during IS110 excision.
Reads that aligned half in the bait and half in the genomic flanking were
recorded separately as genome-head or tail-genome reads, marking the
integration boundary. Junction extraction was implemented in a custom Python
module (`cycle.circle_detect`) built on `pysam` (v0.23.3).

Partial circular intermediates — sub-regions of the IS that excise to form
smaller circles — were detected from split-read back-jump signals
(`cycle.partial_circle`). A single-copy IS reference was built from all called
insertions and reads were re-mapped with `minimap2` (`-x map-ont`); the
name-sorted BAM was parsed with `pysam` to group primary and supplementary
alignments by read. For each read with two or more segments aligning to the
same IS, consecutive segments were ordered by query position and inspected for
a back-jump on the IS reference — that is, segment B beginning on the IS
reference before segment A ends. Each back-jump defined a candidate partial
circle spanning `[start_B, end_A]`. Candidates were retained when the circle
was at least 50 bp and no more than 80 % of the IS length, with at least 10 bp
of overlap on each side of the junction. Independent read breakpoints were
clustered by single linkage with a 5 bp tolerance, and clusters supported by
at least two reads were reported as partial circle calls.

The final output per insertion was a structured JSON record containing the
Sniffles consensus (and, when local assembly succeeded, the polished assembled
sequence with genomic flanking), Prodigal-predicted ORFs, HMM-derived IS110
domain annotations (DEDD, Tnp20, and whether they share an ORF), full
circular-intermediate evidence (tail-head, genome-head, tail-genome read
counts), and clustered partial-circle calls with per-cluster breakpoint
coordinates and supporting-read counts.

## Data sources

- **Long-read FASTQs** — ENA / SRA, retrieved via `kingfisher` (v0.4.1)
- **Reference genomes** — NCBI RefSeq / GenBank, retrieved via NCBI Datasets CLI (v18.13.0)
- **IS110 profile HMMs** — `DEDD_Tnp_IS110.hmm` and `Transposase_20.hmm`,
  HMMER3/f v3.3 (Nov 2019), derived from ISfinder family models
