#!/bin/bash
#SBATCH --job-name=ecor_pc
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch sideproject/run_ecor_partial_circle.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

conda run -p "$CYCLE_ENV_PREFIX" python $CYCLE_ROOT/sideproject/run_ecor_partial_circle.py
