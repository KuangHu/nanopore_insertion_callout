# IS110 Circular Intermediate Detection — Method Summary

## Pipeline

Long-read sequencing data (Oxford Nanopore or PacBio) in FASTQ format are aligned
to a closed reference genome with `minimap2` (v2.30-r1287, `map-ont` or `map-pb`
preset) followed by `samtools sort` (v1.17) to produce a coordinate-sorted BAM.
Non-reference insertions are called from the BAM with Sniffles2 (v2.7.2,
`--minsupport 3 --minsvlen 100`), yielding for each insertion a UUID, genomic
coordinates, supporting reads, and a consensus insertion sequence. Each
consensus is then assembled and polished locally with `miniasm` (v0.3-r179) and
`minipolish` (v0.2.0), and up to 80 bp of upstream and downstream genomic
flanking is recovered.

Open reading frames are predicted within each assembled insertion with Prodigal
(v2.6.3, `-p meta`). Translated proteins are searched with HMMER `hmmsearch`
(v3.4, E-value 1e-5) against two profile HMMs in HMMER3/f format (v3.3,
November 2019): the IS110-family DEDD recombinase profile
(`DEDD_Tnp_IS110.hmm`) and the IS110-family transposase profile
(`Transposase_20.hmm`). An insertion is annotated as an IS110-family element
when either (i) a single ORF is hit by both HMMs, or (ii) two ORFs separated
by no more than 300 bp are hit by the two HMMs respectively (split DEDD–Tnp20
architectures).

To detect circular intermediates, a `[IS][IS]` concatemer "bait" is constructed
by duplicating each Sniffles consensus end to end. All reads are re-aligned to
the bait with `minimap2` (`-x map-ont`), and reads spanning the head-to-tail
junction with at least 100 bp on each side are counted as tail-head junction
reads, providing direct evidence of a full circular intermediate. Reads with one
half in the bait and one half in the genomic flanking are recorded as
genome-head or tail-genome reads, marking the integration boundary. Junction
counting is implemented in a custom Python module (`cycle.circle_detect`,
`pysam` v0.23.3).

Partial circular intermediates — sub-regions of the IS that excise to form
smaller circles — are detected from split-read back-jump signals. A single-copy
IS reference is built from the called insertions and reads are re-mapped with
`minimap2` (`-x map-ont`). For each read with two or more alignment segments on
the same IS, consecutive segments are inspected for a back-jump in which segment
B begins on the IS reference before segment A ends; the partial circle
corresponds to the interval `[start_B, end_A]`. Candidate circles between 50 bp
and 80 % of the IS length, with at least 10 bp of overlap on each side of the
junction, are clustered across reads by single linkage (5 bp tolerance) and
retained when supported by two or more independent reads
(`cycle.partial_circle`, `pysam` v0.23.3).

The output for each insertion is a structured JSON record containing the
assembled IS sequence, flanking regions, ORFs with their IS110 domain
annotations, tail-head junction read counts (full circle evidence), and
clustered partial-circle calls (partial circle evidence).
