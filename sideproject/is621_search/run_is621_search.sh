#!/bin/bash
#SBATCH --job-name=is621_search
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch sideproject/is621_search/run_is621_search.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

OUT=$CYCLE_DATA_ROOT/is621_search
QUERY=$CYCLE_ROOT/sideproject/is621_search/is621_query.faa
mkdir -p $OUT

# NOTE: is110_all_006_070/all_proteins_hmm.faa was deleted in the approved 2026-11-17 cleanup (superseded);
# to re-run, point DB_FASTA at is110_all_000_084/all_proteins_hmm.faa.
DB_FASTA=$CYCLE_DATA_ROOT/is110_all_006_070/all_proteins_hmm.faa
echo "Database: $DB_FASTA  ($(grep -c '^>' $DB_FASTA) proteins)"

echo "=== Building BLAST DB (no parse_seqids) ==="
makeblastdb -in $DB_FASTA -dbtype prot -out $OUT/db

echo "=== BLASTp ==="
blastp \
    -query $QUERY \
    -db $OUT/db \
    -out $OUT/is621_hits.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-5 \
    -num_threads ${SLURM_CPUS_PER_TASK:-16} \
    -max_target_seqs 100000

echo "=== Results ==="
echo "Total hits: $(wc -l < $OUT/is621_hits.tsv)"
echo "Unique IS elements: $(cut -f2 $OUT/is621_hits.tsv | sed 's/__.*$//' | sort -u | wc -l)"

echo ""
echo "Top 20 hits by identity:"
sort -k3,3 -rn $OUT/is621_hits.tsv | head -20 | awk -F'\t' '{printf "  %.1f%%  cov=%d/%d  %s\n", $3, $4, $13, $2}'
