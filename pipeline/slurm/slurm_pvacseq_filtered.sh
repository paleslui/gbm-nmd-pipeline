#!/usr/bin/env bash
#SBATCH --job-name=gbm-nf-filt
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=96:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=__BASE__/logs/Slurm-%j.out
#SBATCH --error=__BASE__/logs/Slurm-%j.err

# =============================================================================
# Run Nextflow pVACseq pipeline on FILTERED VCFs (NMD-relevant variants only).
# Run the filter first (pipeline/scripts/filter_nmd_relevant.sh) to populate
# data/vcf_filtered_decompressed/ before launching this script.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source __BASE__/config.sh

# Override: use FILTERED VCFs (validated against gold standard 11_T)
VCF_DIR=__BASE__/data/vcf_filtered_decompressed

module purge
module load DefaultModules
module load gcc/9.4.0-pe5.34
module load lsfm-init-miniconda/1.0.0
source __BASE__/miniforge3/etc/profile.d/conda.sh
conda activate __BASE__/conda_envs/nf_pvacseq

RUN_TS=${RUN_TS:-$(date '+%Y%m%d_%H%M%S')}
NF_OUTDIR=$RESULTS/nextflow_filtered_${RUN_TS}

echo "[START] Filtered pVACseq run: $(date)"
echo "[INFO] Input VCFs: $VCF_DIR"
echo "[INFO] Output: ${NF_OUTDIR}"

mkdir -p $TMP
cd $TMP

if [ -d "$VEP_CACHE/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38" ]; then
    VEP_CACHE_ARG="--vep_cache $VEP_CACHE --vep_cache_version $VEP_CACHE_VERSION"
else
    echo "[INFO] VEP cache not found — Nextflow will download v${VEP_CACHE_VERSION} to $VEP_CACHE"
    VEP_CACHE_ARG="--vep_cache_version $VEP_CACHE_VERSION"
fi

if [ -d "$VEP_PLUGINS" ] && [ "$(ls -A $VEP_PLUGINS)" ]; then
    VEP_PLUGINS_ARG="--vep_plugins $VEP_PLUGINS"
else
    echo "[INFO] VEP plugins not found — Nextflow will download them"
    VEP_PLUGINS_ARG=""
fi

__BASE__/conda_envs/nf_pvacseq/bin/nextflow run $NF_PIPELINE/main.nf \
    -profile conda \
    -c $NF_PIPELINE/slurm.config \
    -work-dir $NEXTFLOW_WORK \
    --input $VCF_DIR \
    --hla_csv $HLA_CSV \
    --fasta $FASTA \
    --outdir ${NF_OUTDIR} \
    $VEP_CACHE_ARG \
    $VEP_PLUGINS_ARG \
    --pvacseq_iedb $IEDB \
    --pvacseq_algorithm "$PVACSEQ_ALGORITHMS" \
    --validate_params false

echo "[DONE] Filtered pVACseq run: $(date)"
conda deactivate
