#!/bin/bash
#SBATCH --job-name=gtdb_blastdb
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch --needs-avx2 sideproject/build_gtdb_blastdb.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

DB_DIR=$(dirname "$CYCLE_GTDB_BLASTDB")

echo "=== Building BLAST database (combined FASTA already exists) ==="
ls -lh $DB_DIR/gtdb_reps_r220.fna

# Clean up chunks from previous run
rm -f $DB_DIR/chunk_* $DB_DIR/genome_list.txt

echo "=== Running makeblastdb ==="
makeblastdb \
    -in $DB_DIR/gtdb_reps_r220.fna \
    -dbtype nucl \
    -out $DB_DIR/gtdb_reps_r220 \
    -title "GTDB r220 representative genomes" \
    -parse_seqids

echo "=== Done ==="
ls -lh $DB_DIR/gtdb_reps_r220.*
