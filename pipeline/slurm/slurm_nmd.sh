#!/usr/bin/env bash
#SBATCH --job-name=gbm-nmd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=__BASE__/logs/Slurm-%j.out
#SBATCH --error=__BASE__/logs/Slurm-%j.err

# =============================================================================
# Stage 3: nmd_scoring.py per sample + nmd_cohort_summary.py
# Output: $RUN_DIR/3_nmd_analysis/{per_sample,cohort}/
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
#    (Python 3.9, no scipy) on PATH; if `conda activate` silently fails, bare
#    `python` falls through to it. Hard-fail on bad activation, smoke-test the
#    critical imports, and always call the interpreter by absolute path. Do NOT
#    revert. See slurm_python.sh for the full rationale.
if [ -z "${CONDA_PREFIX:-}" ] || [ ! -x "$CONDA_PREFIX/bin/python" ]; then
    echo "[ERROR] CONDA_PREFIX not set or python missing — env activation failed" >&2
    exit 1
fi
echo "[INFO] Using python: $CONDA_PREFIX/bin/python"
$CONDA_PREFIX/bin/python -c "import scipy, pandas, numpy, matplotlib, seaborn, cyvcf2" 2>&1 || {
    echo "[ERROR] Critical Python imports failed — nf_pvacseq env is broken" >&2
    exit 1
}

# Defaults if run standalone
RUN_TS=${RUN_TS:-$(date '+%Y%m%d_%H%M%S')}
RUN_DIR=${RUN_DIR:-$RESULTS/run_${RUN_TS}}
NF_OUTDIR=$RUN_DIR/2_pvacseq
NMD_DIR=$RUN_DIR/3_nmd_analysis
PER_SAMPLE_DIR=$NMD_DIR/per_sample
COHORT_DIR=$NMD_DIR/cohort

echo "[START] Stage 3 (NMD scoring): $(date)"
echo "[INFO] Reading pVACseq results from: $NF_OUTDIR"
echo "[INFO] Per-sample output:           $PER_SAMPLE_DIR"
echo "[INFO] Cohort output:               $COHORT_DIR"

# Backward compat: tolerate the old layout (results/nextflow_filtered_${RUN_TS}/)
if [ ! -d "$NF_OUTDIR/pvactools" ] && [ -d "$RESULTS/nextflow_filtered_${RUN_TS}/pvactools" ]; then
    NF_OUTDIR=$RESULTS/nextflow_filtered_${RUN_TS}
    echo "[WARN] Falling back to old layout: $NF_OUTDIR"
fi

FILTERED_TSVS=$(find ${NF_OUTDIR}/pvactools -name "*.filtered.tsv" 2>/dev/null)
if [ -z "$FILTERED_TSVS" ]; then
    echo "[ERROR] No filtered TSV files found in ${NF_OUTDIR}/pvactools"
    exit 1
fi
echo "[INFO] Found $(echo "$FILTERED_TSVS" | wc -l) filtered TSV(s)"

mkdir -p $PER_SAMPLE_DIR

# ── Per-sample NMD scoring ──────────────────────────────────────────────────
for TSV in $FILTERED_TSVS; do
    SAMPLE=$(basename $(dirname $(dirname $TSV)))
    echo "[INFO] Scoring sample: ${SAMPLE}"

    VCF_SAMPLE=$(echo $SAMPLE | tr -d '_')
    VEP_VCF=${NF_OUTDIR}/ensemblvep/${VCF_SAMPLE}-ensemble-annotated_vep.vcf.gz
    [ ! -f "$VEP_VCF" ] && VEP_VCF=""

    if [ -z "$VEP_VCF" ]; then
        echo "[WARN] No VEP VCF for ${SAMPLE} — rule-based scoring only"
    fi

    OUT_DIR=$PER_SAMPLE_DIR/${SAMPLE}
    mkdir -p $OUT_DIR

    if [ -n "$VEP_VCF" ]; then
        $CONDA_PREFIX/bin/python $PIPELINE/nmd_scoring.py \
            --pvacseq_tsv $TSV \
            --vep_vcf $VEP_VCF \
            --out_dir $OUT_DIR || { echo "[ERROR] nmd_scoring.py failed for ${SAMPLE}" >&2; exit 1; }
    else
        $CONDA_PREFIX/bin/python $PIPELINE/nmd_scoring.py \
            --pvacseq_tsv $TSV \
            --out_dir $OUT_DIR || { echo "[ERROR] nmd_scoring.py failed for ${SAMPLE}" >&2; exit 1; }
    fi
done

# ── Cohort summary ──────────────────────────────────────────────────────────
echo "[INFO] Generating cohort-level summary report..."
mkdir -p $COHORT_DIR
$CONDA_PREFIX/bin/python $PIPELINE/nmd_cohort_summary.py \
    --input_dir $PER_SAMPLE_DIR \
    --out_dir $COHORT_DIR \
    --stage1_summary $RUN_DIR/1_gbm_analysis/sample_mutation_summary.tsv \
    --paired_overlap $RUN_DIR/1_gbm_analysis/paired_variant_overlap.tsv || {
    echo "[ERROR] Stage 3 nmd_cohort_summary.py FAILED — aborting" >&2
    exit 1
}

echo "[INFO] Cohort report: $COHORT_DIR/cohort_report.html"
echo "[DONE] Stage 3 (NMD scoring): $(date)"
conda deactivate
