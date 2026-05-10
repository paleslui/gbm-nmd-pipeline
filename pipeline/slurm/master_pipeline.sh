#!/usr/bin/env bash
#SBATCH --job-name=gbm-master
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --mail-type=ALL
#SBATCH --output=__BASE__/logs/Slurm-%j.out
#SBATCH --error=__BASE__/logs/Slurm-%j.err

# =============================================================================
# GBM NMD-Neoantigen Pipeline — Master Script
# =============================================================================
# Submits three SLURM jobs in sequence using afterok dependencies.
# Usage: sbatch master_pipeline.sh
#
# Output structure (all in $RESULTS/run_${RUN_TS}/):
#   1_gbm_analysis/   Stage 1: mutation landscape (gbm_analysis.py)
#   2_pvacseq/        Stage 2: VEP + pVACseq filtered Nextflow run
#   3_nmd_analysis/   Stage 3: per-sample NMD scoring + cohort summary
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source __BASE__/config.sh

RUN_TS=$(date '+%Y%m%d_%H%M%S')
RUN_DIR=$RESULTS/run_${RUN_TS}

echo "============================================================"
echo "  GBM NMD-Neoantigen Pipeline"
echo "  Run timestamp: ${RUN_TS}"
echo "  Run dir:       ${RUN_DIR}"
echo "  Base:          ${BASE}"
echo "============================================================"

mkdir -p $LOGS $RUN_DIR

JOB1=$(sbatch \
    --partition=$SLURM_PARTITION \
    --constraint=$SLURM_CONSTRAINT \
    --mail-user=$SLURM_MAIL \
    --export=ALL,RUN_TS=${RUN_TS},RUN_DIR=${RUN_DIR} \
    --parsable \
    ${SLURM_SCRIPTS}/slurm_python.sh)
echo "[SUBMITTED] Job 1 (gbm-analysis):  SLURM job ${JOB1}"

JOB2=$(sbatch \
    --partition=$SLURM_PARTITION \
    --constraint=$SLURM_CONSTRAINT \
    --mail-user=$SLURM_MAIL \
    --dependency=afterok:${JOB1} \
    --export=ALL,RUN_TS=${RUN_TS},RUN_DIR=${RUN_DIR} \
    --parsable \
    ${SLURM_SCRIPTS}/slurm_pvacseq_filtered.sh)
echo "[SUBMITTED] Job 2 (gbm-nextflow):  SLURM job ${JOB2} (depends on ${JOB1})"

JOB3=$(sbatch \
    --partition=$SLURM_PARTITION \
    --constraint=$SLURM_CONSTRAINT \
    --mail-user=$SLURM_MAIL \
    --dependency=afterok:${JOB2} \
    --export=ALL,RUN_TS=${RUN_TS},RUN_DIR=${RUN_DIR} \
    --parsable \
    ${SLURM_SCRIPTS}/slurm_nmd.sh)
echo "[SUBMITTED] Job 3 (gbm-nmd):       SLURM job ${JOB3} (depends on ${JOB2})"

echo ""
echo "============================================================"
echo "  Pipeline queued. Monitor: watch -n 30 'squeue -u $USER'"
echo "  Outputs all under: ${RUN_DIR}/"
echo "    Stage 1: ${RUN_DIR}/1_gbm_analysis/"
echo "    Stage 2: ${RUN_DIR}/2_pvacseq/"
echo "    Stage 3: ${RUN_DIR}/3_nmd_analysis/"
echo "  Logs:    $LOGS/"
echo "============================================================"
