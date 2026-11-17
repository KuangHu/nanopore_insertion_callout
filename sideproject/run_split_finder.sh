#!/bin/bash
#SBATCH --job-name=is110_split
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch sideproject/run_split_finder.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

cd "$CYCLE_ROOT"
conda run -n "$CYCLE_CONDA_ENV" python sideproject/is110_split_finder.py \
    --is110-fasta $CYCLE_DATA_ROOT/is110_all_006_026/is110_consensus.fna \
    --genome-dir $CYCLE_GTDB_GENOMES_DIR \
    --output-dir $CYCLE_DATA_ROOT/is110_split_gtdb \
    --threads ${SLURM_CPUS_PER_TASK:-48} \
    --batch-size 2000
