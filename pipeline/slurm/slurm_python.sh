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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source __BASE__/config.sh

module purge
module load DefaultModules
module load gcc/9.4.0-pe5.34
module load lsfm-init-miniconda/1.0.0
conda activate base

pip install --quiet --user cyvcf2 pandas matplotlib seaborn

RUN_TS=${RUN_TS:-$(date '+%Y%m%d_%H%M%S')}

echo "[START] Stage 1-4 (gbm_analysis): $(date)"
echo "[INFO] Run timestamp: ${RUN_TS}"

python $PIPELINE/gbm_analysis.py \
    --vcf_dir $VCF_DIR \
    --out_dir $RESULTS/gbm_analysis_${RUN_TS} \
    --hla_dir $HLA_DIR

echo "[DONE] Stage 1-4 (gbm_analysis): $(date)"
conda deactivate