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

# ── level-3 reproducibility fix: never trust bare `python` on the compute node.
#    `module load lsfm-init-miniconda/1.0.0` puts the cluster miniconda3-4.12.0
#    (Python 3.9, no scipy) on PATH; if `conda activate` above ever silently
#    fails, bare `python` falls through to it and Stage 1 dies with
#    "No module named 'scipy'". We (1) hard-fail if activation didn't set
#    CONDA_PREFIX, (2) smoke-test the critical imports, and (3) always invoke
#    the interpreter by absolute path ($CONDA_PREFIX/bin/python). Do NOT revert.
if [ -z "${CONDA_PREFIX:-}" ] || [ ! -x "$CONDA_PREFIX/bin/python" ]; then
    echo "[ERROR] CONDA_PREFIX not set or python missing — env activation failed" >&2
    exit 1
fi
echo "[INFO] Using python: $CONDA_PREFIX/bin/python"
$CONDA_PREFIX/bin/python -c "import scipy, pandas, numpy, matplotlib, seaborn, cyvcf2" 2>&1 || {
    echo "[ERROR] Critical Python imports failed — nf_pvacseq env is broken" >&2
    exit 1
}

# Defaults if run standalone (not via master_pipeline.sh)
RUN_TS=${RUN_TS:-$(date '+%Y%m%d_%H%M%S')}
RUN_DIR=${RUN_DIR:-$RESULTS/run_${RUN_TS}}
OUT_DIR=$RUN_DIR/1_gbm_analysis

echo "[START] Stage 1 (gbm_analysis): $(date)"
echo "[INFO] Run timestamp: ${RUN_TS}"
echo "[INFO] Output: ${OUT_DIR}"

mkdir -p $OUT_DIR

$CONDA_PREFIX/bin/python $PIPELINE/gbm_analysis.py \
    --vcf_dir $VCF_DIR \
    --out_dir $OUT_DIR \
    --hla_dir $HLA_DIR || {
    # level-3 fix: propagate failure so afterok dependency stops the pipeline
    # instead of cascading Stage 2/3 onto missing Stage 1 outputs.
    echo "[ERROR] Stage 1 gbm_analysis.py FAILED — aborting" >&2
    exit 1
}

echo "[DONE] Stage 1 (gbm_analysis): $(date)"
conda deactivate
