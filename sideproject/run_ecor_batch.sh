#!/bin/bash
#SBATCH --job-name=ecor_batch
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch sideproject/run_ecor_batch.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

cd "$CYCLE_ROOT"
conda run -n "$CYCLE_CONDA_ENV" python sideproject/run_ecor_batch.py \
    --fastq-dir $CYCLE_ECOR_DIR/fastq \
    --ref-dir $CYCLE_ECOR_DIR/reference_genomes \
    --output-dir $CYCLE_ECOR_DIR/output \
    --threads ${SLURM_CPUS_PER_TASK:-48}
