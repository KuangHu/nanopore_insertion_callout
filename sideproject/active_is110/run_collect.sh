#!/bin/bash
#SBATCH --job-name=active_is110
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch --needs-avx2 sideproject/active_is110/run_collect.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

conda run -n "$CYCLE_CONDA_ENV" python $CYCLE_ROOT/sideproject/active_is110/collect_active_is110.py
