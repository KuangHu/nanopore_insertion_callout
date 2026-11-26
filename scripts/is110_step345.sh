#!/bin/bash
#SBATCH --job-name=is110_step345
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch scripts/is110_step345.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

conda run -p "$CYCLE_ENV_PREFIX" python $CYCLE_ROOT/scripts/is110_step345.py
