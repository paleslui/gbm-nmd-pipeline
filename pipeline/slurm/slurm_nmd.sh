#!/usr/bin/env bash
#SBATCH --job-name=gbm-nmd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
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

pip install --quiet --user pandas

RUN_TS=${RUN_TS:-$(date '+%Y%m%d_%H%M%S')}
NF_OUTDIR=$RESULTS/nextflow_${RUN_TS}

echo "[START] Stage 7 (NMD scoring): $(date)"
echo "[INFO] Reading pVACseq results from: ${NF_OUTDIR}"

FILTERED_TSVS=$(find ${NF_OUTDIR}/pvactools -name "*.filtered.tsv" 2>/dev/null)

if [ -z "$FILTERED_TSVS" ]; then
    echo "[ERROR] No filtered TSV files found in ${NF_OUTDIR}/pvactools"
    exit 1
fi

echo "[INFO] Found $(echo "$FILTERED_TSVS" | wc -l) filtered TSV(s)"

for TSV in $FILTERED_TSVS; do
    SAMPLE=$(basename $(dirname $(dirname $TSV)))
    echo "[INFO] Scoring sample: ${SAMPLE}"

    VCF_SAMPLE=$(echo $SAMPLE | tr -d '_')
    VEP_VCF=${NF_OUTDIR}/ensemblvep/${VCF_SAMPLE}-ensemble-annotated_vep.vcf.gz

    if [ ! -f "$VEP_VCF" ]; then
        VEP_VCF=""
    fi

    if [ -z "$VEP_VCF" ]; then
        echo "[WARN] No VEP VCF found for ${SAMPLE} — rule-based scoring only"
    else
        echo "[INFO] VEP VCF: ${VEP_VCF}"
    fi

    OUT_DIR=$RESULTS/nmd_scoring_${SAMPLE}_${RUN_TS}

    if [ -n "$VEP_VCF" ]; then
        python $PIPELINE/nmd_scoring.py \
            --pvacseq_tsv ${TSV} \
            --vep_vcf ${VEP_VCF} \
            --out_dir ${OUT_DIR}
    else
        python $PIPELINE/nmd_scoring.py \
            --pvacseq_tsv ${TSV} \
            --out_dir ${OUT_DIR}
    fi

    echo "[INFO] Results: ${OUT_DIR}"
done

echo "[DONE] Stage 7 (NMD scoring): $(date)"
conda deactivate