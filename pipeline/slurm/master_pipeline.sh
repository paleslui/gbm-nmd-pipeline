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
#   Job 1: gbm-analysis  — Stages 1-4: mutation landscape (gbm_analysis.py)
#   Job 2: gbm-nextflow  — Stages 5-6: VEP annotation + pVACseq
#   Job 3: gbm-nmd       — Stage 7:    NMD sensitivity scoring
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source __BASE__/config.sh

RUN_TS=$(date '+%Y%m%d_%H%M%S')

echo "============================================================"
echo "  GBM NMD-Neoantigen Pipeline"
echo "  Run timestamp: ${RUN_TS}"
echo "  Base directory: ${BASE}"
echo "============================================================"

mkdir -p $LOGS

JOB1=$(sbatch \
    --partition=$SLURM_PARTITION \
    --constraint=$SLURM_CONSTRAINT \
    --mail-user=$SLURM_MAIL \
    --export=ALL,RUN_TS=${RUN_TS} \
    --parsable \
    ${SLURM_SCRIPTS}/slurm_python.sh)
echo "[SUBMITTED] Job 1 (gbm-analysis):  SLURM job ${JOB1}"

JOB2=$(sbatch \
    --partition=$SLURM_PARTITION \
    --constraint=$SLURM_CONSTRAINT \
    --mail-user=$SLURM_MAIL \
    --dependency=afterok:${JOB1} \
    --export=ALL,RUN_TS=${RUN_TS} \
    --parsable \
    ${SLURM_SCRIPTS}/slurm_pvacseq.sh)
echo "[SUBMITTED] Job 2 (gbm-nextflow):  SLURM job ${JOB2} (depends on ${JOB1})"

JOB3=$(sbatch \
    --partition=$SLURM_PARTITION \
    --constraint=$SLURM_CONSTRAINT \
    --mail-user=$SLURM_MAIL \
    --dependency=afterok:${JOB2} \
    --export=ALL,RUN_TS=${RUN_TS} \
    --parsable \
    ${SLURM_SCRIPTS}/slurm_nmd.sh)
echo "[SUBMITTED] Job 3 (gbm-nmd):       SLURM job ${JOB3} (depends on ${JOB2})"

echo ""
echo "============================================================"
echo "  Pipeline queued. Monitor: watch -n 30 'squeue -u $USER'"
echo "  Outputs:"
echo "    Stages 1-4: $RESULTS/gbm_analysis_${RUN_TS}/"
echo "    Stages 5-6: $RESULTS/nextflow_${RUN_TS}/"
echo "    Stage 7:    $RESULTS/nmd_scoring_*_${RUN_TS}/"
echo "    Logs:       $LOGS/"
echo "============================================================"