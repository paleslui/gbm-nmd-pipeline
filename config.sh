#!/usr/bin/env bash
# =============================================================================
# GBM NMD-Neoantigen Pipeline — Configuration
# =============================================================================
# Thesis: Identifying NMD-Sensitive Immunogenic Mutations in Primary and
#         Recurrent Glioblastoma
# Author: Luigi Palese (ZHAW / UniBe)
#
# HOW TO USE:
#   1. Set BASE to the absolute path of this repository on your system
#   2. Set SLURM_PARTITION and SLURM_PARTITION_GPU to your HPC partitions
#   3. Set SLURM_CONSTRAINT to your node constraint (or leave empty)
#   4. All other paths are derived automatically
#
# This file is sourced by all SLURM scripts — do not execute it directly.
# =============================================================================

# ── USER SETTINGS — edit these for your HPC environment ─────────────────────

# Absolute path to the root of this repository
# EDIT THIS: absolute path to the cloned repo on your machine
BASE=/path/to/gbm_nmd_pipeline

# SLURM partition for CPU jobs (analysis, VEP, NMD scoring)
SLURM_PARTITION=earth-3

# SLURM partition for GPU jobs (pVACseq with MHCflurry)
SLURM_PARTITION_GPU=earth-4

# SLURM node constraint (leave empty "" if not needed)
SLURM_CONSTRAINT=rhel8

# Email for SLURM notifications
SLURM_MAIL=paleslui@students.zhaw.ch

# Conda environment for Nextflow
CONDA_ENV_NEXTFLOW=nf_pvacseq

# ── DERIVED PATHS — do not edit ──────────────────────────────────────────────

# Pipeline source code
PIPELINE=$BASE/pipeline
NF_PIPELINE=$PIPELINE/nextflow-pvacseq
SLURM_SCRIPTS=$PIPELINE/slurm

# Input data
DATA=$BASE/data
VCF_DIR=$DATA/vcf
HLA_DIR=$DATA/hla_typing
HLA_CSV=$DATA/hla_typing_classI.csv
FASTA=$DATA/reference/GRCh38.primary_assembly.fa

# Resources (downloaded on first run if absent)
RESOURCES=$BASE/resources
IEDB=$RESOURCES/iedb
VEP_PLUGINS=$RESOURCES/VEP_plugins
VEP_CACHE=$RESOURCES/vep_cache
VEP_CACHE_VERSION=113

# Nextflow work directory (large, not in repo)
NEXTFLOW_WORK=$BASE/nextflow_work

# Outputs and logs
RESULTS=$BASE/results
LOGS=$BASE/logs
TMP=$BASE/tmp

# pVACseq algorithms (Class I only — Class II requires separate HLA typing)
PVACSEQ_ALGORITHMS="MHCflurry MHCnuggetsI NetMHC NetMHCpan NetMHCpanEL PickPocket SMM SMMPMBEC"

# TMB normalisation: total callable region in megabases.
# Default 30 Mb assumes standard whole-exome capture.
# Replace with exact value when capture-kit BED is available:
#   bedtools genomecov -bg -i target.bed | awk '{s += $3-$2} END {print s/1e6}'
# Exported so the Stage 1 Python child (gbm_analysis.py via slurm_python.sh)
# can read it from the environment (os.environ["TMB_DENOMINATOR_MB"]).
# Empirical value from resources/S31285117_Padded.bed
# (Agilent SureSelect Exome V7 + 100bp padding, MD5 46e1a67055...)
# Sum of (end - start) across 215,154 intervals / 1e6
export TMB_DENOMINATOR_MB=92.699206


export NXF_HOME=$BASE/.nextflow

export PYTHONUSERBASE=$BASE/.local

export MHCFLURRY_DATA_PATH=$BASE/resources/mhcflurry
