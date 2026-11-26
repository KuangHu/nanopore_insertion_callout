#!/bin/bash
#SBATCH --job-name=is110_blastn
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch sideproject/run_is110_blastn.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

DB=$CYCLE_GTDB_BLASTDB
QUERY=$CYCLE_DATA_ROOT/is110_all_006_026/is110_consensus.fna
OUT_DIR=$CYCLE_DATA_ROOT/is110_split_gtdb
mkdir -p $OUT_DIR

echo "=== BLASTn: 2276 IS110 vs GTDB ==="
blastn \
    -query $QUERY \
    -db $DB \
    -out $OUT_DIR/is110_vs_gtdb.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-10 \
    -num_threads ${SLURM_CPUS_PER_TASK:-48} \
    -max_target_seqs 500 \
    -max_hsps 5

echo "=== Results ==="
echo "Hits: $(wc -l < $OUT_DIR/is110_vs_gtdb.tsv)"
echo "Unique IS110 with hits: $(cut -f1 $OUT_DIR/is110_vs_gtdb.tsv | sort -u | wc -l)"
echo "Unique genomes hit: $(cut -f2 $OUT_DIR/is110_vs_gtdb.tsv | sort -u | wc -l)"

echo "=== Finding split hits ==="
conda run -p "$CYCLE_ENV_PREFIX" python $CYCLE_ROOT/sideproject/find_splits_from_blast.py \
    --blast-results $OUT_DIR/is110_vs_gtdb.tsv \
    --output-dir $OUT_DIR

echo "=== Done ==="
