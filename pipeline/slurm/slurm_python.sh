#!/usr/bin/env bash
#SBATCH --job-name=gbm-analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=__BASE__/logs/Slurm-%j.out
#SBATCH --error=__BASE__/logs/Slurm-%j.err

# =============================================================================
# Stage 1: gbm_analysis.py — mutation landscape, TMZ signature, paired analysis
# Output: $RUN_DIR/1_gbm_analysis/
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source __BASE__/config.sh

module purge
module load DefaultModules
module load gcc/9.4.0-pe5.34
module load lsfm-init-miniconda/1.0.0
source __BASE__/miniforge3/etc/profile.d/conda.sh
conda activate __BASE__/conda_envs/nf_pvacseq

# Defaults if run standalone (not via master_pipeline.sh)
RUN_TS=${RUN_TS:-$(date '+%Y%m%d_%H%M%S')}
RUN_DIR=${RUN_DIR:-$RESULTS/run_${RUN_TS}}
OUT_DIR=$RUN_DIR/1_gbm_analysis

echo "[START] Stage 1 (gbm_analysis): $(date)"
echo "[INFO] Run timestamp: ${RUN_TS}"
echo "[INFO] Output: ${OUT_DIR}"

mkdir -p $OUT_DIR

python $PIPELINE/gbm_analysis.py \
    --vcf_dir $VCF_DIR \
    --out_dir $OUT_DIR \
    --hla_dir $HLA_DIR

echo "[DONE] Stage 1 (gbm_analysis): $(date)"
conda deactivate
