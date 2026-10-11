#!/bin/bash
#SBATCH --job-name=active_is110
#SBATCH --output=/home/kuangh/logs/active_is110_%j.out
#SBATCH --error=/home/kuangh/logs/active_is110_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --qos=standard
#SBATCH --partition=standard
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14,node-48-256g-17,node-48-256g-18,node-48-256g-19,node-48-256g-20

conda run -n opfi python /home/kuangh/tools/Cycle/sideproject/active_is110/collect_active_is110.py
