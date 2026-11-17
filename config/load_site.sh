# shellcheck shell=bash
# Load site settings for the Cycle pipeline. Source this file; do not execute it.
#
# Precedence (highest first):
#   1. variables already exported in the environment
#   2. the file named by $CYCLE_SITE_ENV, else $CYCLE_ROOT/config/site.env (git-ignored, per cluster)
#   3. $CYCLE_ROOT/config/igi.env (the values used on IGI before the 2026-11 migration)
# Every value in the env files is written as VAR="${VAR:-default}", so step 1 always wins.
#
# After sourcing, CYCLE_SBATCH_SITE is a bash array with the site's sbatch options
# (--partition, --qos, --account, extra flags), ready to splice into sbatch calls.

if [ -z "${CYCLE_ROOT:-}" ]; then
    CYCLE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
export CYCLE_ROOT

_cycle_site_file="${CYCLE_SITE_ENV:-$CYCLE_ROOT/config/site.env}"
if [ ! -f "$_cycle_site_file" ]; then
    if [ -n "${CYCLE_SITE_ENV:-}" ]; then
        echo "[cycle] CYCLE_SITE_ENV=$CYCLE_SITE_ENV does not exist" >&2
        return 1 2>/dev/null || exit 1
    fi
    _cycle_site_file="$CYCLE_ROOT/config/igi.env"
    [ -n "${CYCLE_QUIET:-}" ] || echo "[cycle] no config/site.env; using IGI defaults from config/igi.env" >&2
fi
# shellcheck source=/dev/null
. "$_cycle_site_file"
unset _cycle_site_file

: "${CYCLE_BATCH_METADATA_DIR:=$CYCLE_ROOT/data/batches}"
: "${CYCLE_ISFINDER_FASTA:=$CYCLE_DATA_ROOT/batch_000/is_reference/ISfinder_raw.fna}"
export CYCLE_DATA_ROOT CYCLE_LOG_DIR CYCLE_CONDA_ENV CYCLE_ENV_PREFIX CYCLE_HMM_DIR \
       CYCLE_BATCH_METADATA_DIR CYCLE_ISFINDER_FASTA CYCLE_ECOR_DIR CYCLE_GTDB_GENOMES_DIR CYCLE_GTDB_BLASTDB \
       CYCLE_SLURM_ACCOUNT CYCLE_SLURM_PARTITION CYCLE_SLURM_QOS CYCLE_SLURM_CPUS CYCLE_SLURM_MEM \
       CYCLE_SLURM_AVX2_EXCLUDE CYCLE_SLURM_EXTRA

CYCLE_SBATCH_SITE=()
[ -n "${CYCLE_SLURM_PARTITION:-}" ] && CYCLE_SBATCH_SITE+=("--partition=$CYCLE_SLURM_PARTITION")
[ -n "${CYCLE_SLURM_QOS:-}" ] && CYCLE_SBATCH_SITE+=("--qos=$CYCLE_SLURM_QOS")
[ -n "${CYCLE_SLURM_ACCOUNT:-}" ] && CYCLE_SBATCH_SITE+=("--account=$CYCLE_SLURM_ACCOUNT")
if [ -n "${CYCLE_SLURM_EXTRA:-}" ]; then
    # shellcheck disable=SC2206
    CYCLE_SBATCH_SITE+=($CYCLE_SLURM_EXTRA)
fi
