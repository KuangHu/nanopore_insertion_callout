#!/bin/bash
#SBATCH --job-name=is110_collect
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch scripts/is110_collect_simple.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

conda run -p "$CYCLE_ENV_PREFIX" python $CYCLE_ROOT/scripts/is110_collect_simple.py
