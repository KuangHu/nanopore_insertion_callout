#!/bin/bash
#SBATCH --job-name=plot_splits
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch sideproject/run_plot_splits.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

conda run -n "$CYCLE_CONDA_ENV" python $CYCLE_ROOT/sideproject/plot_split_example.py
