#!/usr/bin/env bash
# =============================================================================
# GBM NMD-Neoantigen Pipeline — One-Time Setup Script
# =============================================================================
# Run this ONCE before your first pipeline run:
#   tmux new -s gbm-setup
#   bash setup.sh 2>&1 | tee logs/setup.log
#
# This script will:
#   1. Validate config.sh paths
#   2. Install miniforge into $BASE/miniforge3
#   3. Create nf_pvacseq conda env (Nextflow + Java via mamba)
#   4. Build all Nextflow pipeline conda environments using mamba
#   5. Download VEP cache v113 into $RESOURCES/vep_cache
#   6. Patch vep/main.nf to unset PERL5LIB (Perl conflict fix)
#   7. Create required conda symlinks
#   8. Update pipeline scripts with correct paths
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

mkdir -p $LOGS
exec > >(tee -a $LOGS/setup.log) 2>&1

echo "============================================================"
echo "  GBM NMD-Neoantigen Pipeline — Setup"
echo "  $(date)"
echo "  Base: ${BASE}"
echo "============================================================"

# ── Step 1: Validate required paths ──────────────────────────────────────────
echo ""
echo "[STEP 1] Validating paths in config.sh..."

ERRORS=0
[ ! -d "$DATA/vcf" ]    && echo "[ERROR] VCF dir not found: $DATA/vcf"        && ERRORS=$((ERRORS+1))
[ ! -f "$HLA_CSV" ]     && echo "[ERROR] HLA CSV not found: $HLA_CSV"         && ERRORS=$((ERRORS+1))
[ ! -f "$FASTA" ]       && echo "[ERROR] Reference not found: $FASTA"         && ERRORS=$((ERRORS+1))
# IEDB is downloaded automatically by the Nextflow pipeline on first run
[ ! -d "$VEP_PLUGINS" ] && echo "[ERROR] VEP plugins not found: $VEP_PLUGINS" && ERRORS=$((ERRORS+1))

if [ $ERRORS -gt 0 ]; then
    echo "[ERROR] $ERRORS missing path(s). Fix config.sh and rerun."
    exit 1
fi
echo "[OK] All required paths exist"


# ── Step 1b: Configure environment variables ──────────────────────────────────
echo ""
echo "[STEP 1b] Configuring environment variables..."

# Redirect Nextflow home to scratch (avoids home quota issues)
if ! grep -q "NXF_HOME" $BASE/config.sh; then
    printf "\nexport NXF_HOME=$BASE/.nextflow\n" >> $BASE/config.sh
fi
# Redirect pip installs to scratch
if ! grep -q "PYTHONUSERBASE" $BASE/config.sh; then
    printf "\nexport PYTHONUSERBASE=$BASE/.local\n" >> $BASE/config.sh
fi
source $BASE/config.sh
export NXF_HOME=$BASE/.nextflow
export PYTHONUSERBASE=$BASE/.local
echo "[OK] Environment variables configured"
# ── Step 2: Install miniforge ─────────────────────────────────────────────────
echo ""
echo "[STEP 2] Installing miniforge..."

MINIFORGE_DIR=$BASE/miniforge3

if [ -f "$MINIFORGE_DIR/bin/mamba" ]; then
    echo "[OK] Miniforge already installed"
else
    echo "[INFO] Downloading and installing miniforge..."
    curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
        -o /tmp/miniforge_setup.sh
    bash /tmp/miniforge_setup.sh -b -p $MINIFORGE_DIR
    rm /tmp/miniforge_setup.sh
    echo "[OK] Miniforge installed at $MINIFORGE_DIR"
fi

source $MINIFORGE_DIR/etc/profile.d/conda.sh
conda activate base

# ── Step 3: Create nf_pvacseq environment ────────────────────────────────────
echo ""
echo "[STEP 3] Creating nf_pvacseq conda environment..."

NF_ENV_PATH=$BASE/conda_envs/nf_pvacseq
mkdir -p $BASE/conda_envs

if [ -d "$NF_ENV_PATH" ]; then
    echo "[OK] nf_pvacseq env already exists"
else
    echo "[INFO] Building nf_pvacseq env with mamba (Nextflow + Java)..."
    mamba create -y -p $NF_ENV_PATH \
        -c conda-forge -c bioconda \
        nextflow=25.10.* openjdk python=3.11 samtools
    echo "[OK] nf_pvacseq environment created"
fi

NF_BIN=$NF_ENV_PATH/bin/nextflow

# ── Index reference genome (required for TMZ signature analysis) ──────────────
echo ""
echo "[STEP 3b] Indexing reference genome for TMZ signature analysis..."
if [ ! -f "$FASTA" ]; then
    echo "[WARN]  FASTA not found — skipping indexing (TMZ signature unavailable)"
elif [ -f "${FASTA}.fai" ]; then
    echo "[OK] Reference index already exists"
else
    echo "[INFO] Running samtools faidx on reference FASTA..."
    $NF_ENV_PATH/bin/samtools faidx $FASTA
    echo "[OK] Reference genome indexed"
fi
echo "[INFO] Nextflow: $($NF_BIN -version 2>&1 | head -1)"

# ── Step 4: Build Nextflow pipeline conda environments ───────────────────────
echo ""
echo "[STEP 4] Building Nextflow pipeline conda environments..."

mkdir -p $NEXTFLOW_WORK/conda

NF_PIPELINE=$BASE/pipeline/nextflow-pvacseq

# VEP environment
VEP_ENV=$NEXTFLOW_WORK/conda/env-7093d00cc21f39945dcdbeec541c6ae5
if [ -d "$VEP_ENV" ] || [ -L "$VEP_ENV" ]; then
    echo "[OK] VEP env already exists"
else
    echo "[INFO] Building VEP environment (~20-30 min)..."
    mamba env create --yes -p $VEP_ENV \
        --file $NF_PIPELINE/modules/local/vep/environment.yml
    echo "[OK] VEP environment built"
fi

# pVACseq environment
PVACSEQ_ENV=$NEXTFLOW_WORK/conda/env-4a82d2b85d92a884189ec0284fcb54f7
if [ -d "$PVACSEQ_ENV" ] || [ -L "$PVACSEQ_ENV" ]; then
    echo "[OK] pVACseq env already exists"
else
    echo "[INFO] Building pVACseq environment (~20-30 min)..."
    mamba env create --yes -p $PVACSEQ_ENV \
        --file $NF_PIPELINE/modules/local/configure_pvacseq/environment.yml
    echo "[OK] pVACseq environment built"
fi

# MultiQC environment
MULTIQC_ENV=$NEXTFLOW_WORK/conda/env-55e377717f27765e46a811a08ed80f85
if [ -d "$MULTIQC_ENV" ] || [ -L "$MULTIQC_ENV" ]; then
    echo "[OK] MultiQC env already exists"
else
    echo "[INFO] Building MultiQC environment (~5 min)..."
    mamba env create --yes -p $MULTIQC_ENV \
        --file $NF_PIPELINE/modules/nf-core/multiqc/environment.yml
    echo "[OK] MultiQC environment built"
fi


# ── Step 4b: Patch pvactools to use sys.executable for MHCnuggetsI ───────────
echo ""
echo "[STEP 4b] Patching pvactools MHCnuggets to use correct Python..."
PRED_CLASS=$PVACSEQ_ENV/lib/python3.11/site-packages/pvactools/lib/prediction_class.py
if grep -q "sys.executable" $PRED_CLASS; then
    echo "[OK] pvactools already patched"
else
    sed -i 's/arguments = \["python", script,/arguments = [sys.executable, script,/' $PRED_CLASS
    echo "[OK] pvactools patched"
fi

# ── Step 4c: Download MHCflurry models ───────────────────────────────────────
echo ""
echo "[STEP 4c] Downloading MHCflurry models..."
MHCFLURRY_DIR=$BASE/resources/mhcflurry
if [ -d "$MHCFLURRY_DIR/4" ]; then
    echo "[OK] MHCflurry models already present"
else
    mkdir -p $MHCFLURRY_DIR
    # mhcflurry ignores MHCFLURRY_DATA_PATH during download — fetch to default then copy
    $PVACSEQ_ENV/bin/mhcflurry-downloads fetch models_class1_pan
    DEFAULT_MHCFLURRY=$(python3 -c "import appdirs; print(appdirs.user_data_dir('mhcflurry'))" 2>/dev/null || echo "$HOME/.local/share/mhcflurry")
    cp -r $DEFAULT_MHCFLURRY/4 $MHCFLURRY_DIR/
    echo "[OK] MHCflurry models copied to $MHCFLURRY_DIR"
fi
if ! grep -q "MHCFLURRY_DATA_PATH" $BASE/config.sh; then
    printf "\nexport MHCFLURRY_DATA_PATH=$MHCFLURRY_DIR\n" >> $BASE/config.sh
fi
# ── Step 5: Download VEP cache v113 ──────────────────────────────────────────
echo ""
echo "[STEP 5] Downloading VEP cache v113..."

VEP_CACHE_TARGET=$VEP_CACHE/homo_sapiens/113_GRCh38

if [ -d "$VEP_CACHE_TARGET" ]; then
    echo "[OK] VEP cache v113 already present"
else
    echo "[INFO] Downloading VEP cache (~15GB, ~30 min)..."
    mkdir -p $VEP_CACHE
    conda activate $VEP_ENV
    vep_install \
        --CACHEDIR $VEP_CACHE \
        --SPECIES homo_sapiens \
        --ASSEMBLY GRCh38 \
        --CACHE_VERSION 113 \
        --AUTO c \
        --CONVERT \
        --NO_BIOPERL \
        --NO_HTSLIB \
        --NO_TEST \
        --NO_UPDATE
    conda activate base
    echo "[OK] VEP cache downloaded"
fi

# ── Step 6: Patch vep/main.nf to unset PERL5LIB ──────────────────────────────
echo ""
echo "[STEP 6] Patching vep/main.nf for Perl compatibility..."

VEP_MODULE=$NF_PIPELINE/modules/local/vep/main.nf

python3 << PYEOF
content = open('$VEP_MODULE').read()
if 'unset PERL5LIB' not in content:
    content = content.replace(
        '    vep \\\\\\\\\n        -i',
        '    unset PERL5LIB\n    vep \\\\\\\\\n        -i'
    )
    open('$VEP_MODULE', 'w').write(content)
    print('[OK] vep/main.nf patched')
else:
    print('[OK] vep/main.nf already patched')
PYEOF

# ── Step 7: Create ENSEMBLVEP_DOWNLOAD symlink ────────────────────────────────
echo ""
echo "[STEP 7] Creating ENSEMBLVEP_DOWNLOAD symlink..."

DOWNLOAD_ENV=$NEXTFLOW_WORK/conda/env-2def1406b7b4183115ec06885ae604a4
if [ -L "$DOWNLOAD_ENV" ] || [ -d "$DOWNLOAD_ENV" ]; then
    echo "[OK] ENSEMBLVEP_DOWNLOAD symlink already exists"
else
    ln -s $VEP_ENV $DOWNLOAD_ENV
    echo "[OK] Symlink created"
fi

# ── Step 8: Update pipeline scripts with correct paths ───────────────────────
echo ""
echo "[STEP 8] Updating pipeline scripts..."

# Set conda cacheDir in slurm.config — always overwrite to ensure correct user path
SLURM_CFG=$NF_PIPELINE/slurm.config
sed -i '/conda.cacheDir/d' $SLURM_CFG
sed -i "1s|^|conda.cacheDir = \"$NEXTFLOW_WORK/conda\"\n\n|" $SLURM_CFG
echo "[OK] slurm.config: conda.cacheDir set to $NEXTFLOW_WORK/conda"

# Update slurm.config ENSEMBLVEP_DOWNLOAD conda path
python3 << PYEOF
import re
content = open('$NF_PIPELINE/slurm.config').read()
pattern = r"withName: 'ENSEMBLVEP_DOWNLOAD' \{[^}]*\}"
replacement = "withName: 'ENSEMBLVEP_DOWNLOAD' {\n        conda = '$DOWNLOAD_ENV'\n    }"
new_content = re.sub(pattern, replacement, content)
open('$NF_PIPELINE/slurm.config', 'w').write(new_content)
print('[OK] slurm.config updated with ENSEMBLVEP_DOWNLOAD path')
PYEOF

# Update slurm_pvacseq.sh to use local miniforge and nextflow
PVACSEQ_SLURM=$BASE/pipeline/slurm/slurm_pvacseq.sh

# Replace conda activation with miniforge-based activation
python3 << PYEOF
content = open('$PVACSEQ_SLURM').read()

# Fix conda activation
old = 'module load lsfm-init-miniconda/1.0.0\nconda activate \$CONDA_ENV_NEXTFLOW'
new = 'source $MINIFORGE_DIR/etc/profile.d/conda.sh\nconda activate $NF_ENV_PATH'
content = content.replace(old, new)

# Fix nextflow binary - replace any existing nextflow run with absolute path
import re
content = re.sub(r'(?m)^(\s*)(\S+nextflow) run', r'\1$NF_BIN run', content)

open('$PVACSEQ_SLURM', 'w').write(content)
print('[OK] slurm_pvacseq.sh updated')
PYEOF

# Expand variables in slurm_pvacseq.sh
sed -i "s|\$MINIFORGE_DIR|${MINIFORGE_DIR}|g" $PVACSEQ_SLURM
sed -i "s|\$NF_ENV_PATH|${NF_ENV_PATH}|g" $PVACSEQ_SLURM
sed -i "s|\$NF_BIN|${NF_BIN}|g" $PVACSEQ_SLURM

# Replace __BASE__ placeholder in all slurm scripts with absolute path
for f in $BASE/pipeline/slurm/*.sh; do
    sed -i "s|__BASE__|$BASE|g" "$f"
done
echo "[OK] All slurm scripts: __BASE__ replaced with $BASE"
echo "[OK] Pipeline scripts updated"


echo ""
echo "[STEP 9] Patching pVACseq Nextflow module (cluster-portability fix)..."
# The published nf-core pVACseq module reads IEDB data files from the work dir,
# which on shared parallel filesystems (BeeGFS, Lustre) can hit transient
# open() failures under high concurrent I/O. This patch makes the module copy
# IEDB to node-local scratch before invoking pvacseq run. Idempotent.
python3 "$BASE/pipeline/scripts/patch_pvacseq_module.py" \
    "$NF_PIPELINE/modules/local/pvacseq/main.nf"

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  Setup complete! $(date)"
echo ""
echo "  Run the pipeline:"
echo "    sbatch $BASE/pipeline/slurm/master_pipeline.sh"
echo "============================================================"
