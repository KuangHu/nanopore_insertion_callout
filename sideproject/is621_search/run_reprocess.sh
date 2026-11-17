#!/bin/bash
#SBATCH --job-name=is621_reprocess
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch sideproject/is621_search/run_reprocess.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

conda run -n "$CYCLE_CONDA_ENV" python $CYCLE_ROOT/sideproject/is621_search/reprocess_long_is621.py
