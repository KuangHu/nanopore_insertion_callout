#!/bin/bash
#SBATCH --job-name=is621_all_long
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch sideproject/is621_search/run_visualize_all.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

conda run -p "$CYCLE_ENV_PREFIX" python $CYCLE_ROOT/sideproject/is621_search/visualize_all_long_is621.py
