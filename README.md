# GBM NMD-Neoantigen Pipeline

A computational pipeline for identifying NMD-sensitive immunogenic mutations in primary and recurrent Glioblastoma (GBM). Developed as part of a Master's thesis at ZHAW.

**Thesis:** Identifying NMD-Sensitive Immunogenic Mutations in Primary and Recurrent Glioblastoma  
**Author:** Luigi Palese  
**Supervisors:** Erik Vassella, Tugce Bilgin Sonay, Maria Anisimova  
**Collaborator:** Massimo Maiolo (Inselspital Bern)

---

## Background

Temozolomide (TMZ) chemotherapy induces frameshift mutations in GBM tumor cells that introduce premature termination codons (PTCs). Nonsense-mediated mRNA decay (NMD) degrades these transcripts before they can be translated and presented to the immune system. This pipeline identifies which TMZ-induced mutations produce immunogenic neoantigens that are being silenced by NMD — candidates that could become visible to the immune system upon NMD inhibition.

---

## Pipeline Overview

```
Stage 1-4  gbm_analysis.py      VCF parsing, mutation burden, paired T/M comparison,
                                  gene-level recurrence, HLA typing → HTML report
Stage 5-6  nextflow-pvacseq     VEP v113 annotation (NMD plugin) + pVACseq neoantigen
                                  prediction (MHC Class I, 8 algorithms)
Stage 7    nmd_scoring.py       NMD sensitivity scoring using ensemble of:
                                  - VEP NMD plugin (Ensembl rules)
                                  - Lindeboom et al. 2019 rule-based method
                                  → Priority tiers + HTML report
```

All three stages are submitted automatically by a single master SLURM script with `afterok` dependencies.

---

## Requirements

- SLURM HPC cluster (Linux, RHEL8 or compatible)
- ~50GB scratch space for conda environments + VEP cache
- ~200k inodes available
- Internet access from compute nodes (for miniforge and VEP cache download)
- Your own WES VCF files and HLA typing data (see Input Data below)

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/paleslui/gbm-nmd-pipeline.git
cd gbm-nmd-pipeline
```

### 2. Edit `config.sh`

Open `config.sh` and change **only these two settings**:

```bash
# Absolute path to where you cloned this repository
BASE=/path/to/gbm-nmd-pipeline

# Your HPC SLURM partition names
SLURM_PARTITION=your-cpu-partition
SLURM_PARTITION_GPU=your-gpu-partition
```

Everything else is derived automatically.

### 3. Add your data

```
data/
├── vcf/                    ← your VCF files (see format below)
├── hla_typing/             ← HLA typing files from Dragen
├── hla_typing_classI.csv   ← HLA Class I alleles CSV (see format below)
└── reference/
    └── GRCh38.primary_assembly.fa   ← download from Ensembl
```

### 4. Run setup (one time only, ~1.5 hours)

Setup runs on the login node and is interactive (mamba prompts for env builds). Use `tmux` so an SSH dropout doesn't kill it:

```bash
tmux new -s gbm-setup
bash setup.sh 2>&1 | tee logs/setup.log
# Detach: Ctrl+B then D
# Reattach: tmux attach -t gbm-setup
```

This installs miniforge, builds all conda environments, downloads the VEP v113 cache, and patches all pipeline scripts automatically. You only need to run this once.

### 5. Run the pipeline

```bash
sbatch pipeline/slurm/master_pipeline.sh
```

Monitor progress:
```bash
watch -n 30 'squeue -u $USER'
```

---

## Input Data Format

### VCF files

Place unzipped VCF files in `data/vcf/` using this naming convention:

```
{patient_id}{T|M}-ensemble-annotated.vcf
```

Where `T` = primary (pre-treatment) and `M` = recurrent (post-treatment):

```
11T-ensemble-annotated.vcf    ← patient 11, primary
11M-ensemble-annotated.vcf    ← patient 11, recurrent
45T-ensemble-annotated.vcf
45M-ensemble-annotated.vcf
...
```

Both SnpEff-annotated and VEP-annotated VCFs are supported. VEP-annotated VCFs from Dragen with `--canonical --mane --plugin NMD` are recommended for best results.

### HLA typing CSV

`data/hla_typing_classI.csv` must contain Class I HLA alleles in this format:

```csv
Sample_ID,HLA_Types
11_T,HLA-A*26:01;HLA-A*01:01;HLA-B*08:01;HLA-B*38:01;HLA-C*07:01;HLA-C*12:03
45_T,HLA-A*02:01;HLA-A*03:01;HLA-B*07:02;HLA-B*15:01;HLA-C*03:04;HLA-C*07:02
```

Note: Sample IDs use underscore (e.g. `11_T`), VCF filenames do not (e.g. `11T`).

Only samples present in this CSV will be processed by pVACseq.

---

## Output Structure

Each run produces a timestamped set of output directories:

```
results/
├── gbm_analysis_{RUN_TS}/          ← Stage 1-4 outputs
│   └── run_{timestamp}/
│       ├── report.html             ← interactive HTML report
│       ├── summary_mutation_burden.tsv
│       ├── all_fs_sg_variants.tsv
│       ├── paired_variant_overlap.tsv
│       ├── gene_recurrence.tsv
│       ├── hla_typing_summary.tsv
│       └── plot_*.png
│
├── nextflow_{RUN_TS}/              ← Stage 5-6 outputs
│   ├── ensemblvep/                 ← VEP-annotated VCFs
│   │   └── vep_cache/             ← reused on subsequent runs
│   └── pvactools/
│       └── {sample_id}/
│           └── MHC_Class_I/
│               └── {sample}.filtered.tsv   ← neoantigen candidates
│
└── nmd_scoring_{sample}_{RUN_TS}/ ← Stage 7 outputs
    ├── report_nmd.html             ← NMD scoring report
    └── nmd_scored_candidates.tsv
```

---

## NMD Scoring

Each neoantigen candidate from pVACseq is classified into one of:

| Classification | Meaning |
|---|---|
| **SENSITIVE** | Transcript predicted to be degraded by NMD — candidate exposed by NMD inhibition |
| **INSENSITIVE** | Transcript predicted to escape NMD (e.g. last exon rule) |
| **NOT_APPLICABLE** | Missense variant — NMD does not apply |
| **UNCERTAIN** | Insufficient transcript information |

Candidates are then prioritised into tiers:

| Tier | Criteria |
|---|---|
| **Tier 1** | NMD-SENSITIVE + IC50 < 50 nM — highest priority therapeutic targets |
| **Tier 2** | NMD-SENSITIVE + IC50 50-500 nM |
| **Tier 3** | NMD-INSENSITIVE + strong binder — controls |

Scoring uses an ensemble of the VEP NMD plugin (Ensembl) and the Lindeboom et al. 2019 rule-based method. Both methods must agree for high-confidence classification.

---

## Resources

All resources are downloaded automatically by `setup.sh`:

| Resource | Version | Location after setup |
|---|---|---|
| VEP cache | v113 GRCh38 | `resources/vep_cache/` |
| VEP plugins | current | `resources/VEP_plugins/` |
| IEDB mhc_i | current | `resources/iedb/mhc_i/` |

---

## Configuration Reference

All pipeline settings are in `config.sh`. The only required edits are `BASE` and the SLURM partition names. Other settings you may want to adjust:

| Setting | Default | Description |
|---|---|---|
| `VEP_CACHE_VERSION` | 113 | VEP cache version |
| `PVACSEQ_ALGORITHMS` | 8 Class I algorithms (MHCflurry, MHCnuggetsI, NetMHC, NetMHCpan, NetMHCpanEL, PickPocket, SMM, SMMPMBEC) | Prediction algorithms for pVACseq |
| `SLURM_MAIL` | — | Email for SLURM notifications |

---


## Troubleshooting

**Setup fails with `git clone ... EVP_KDF_ctrl undefined symbol` on a compute node**
Compute nodes have a libcrypto/openssl mismatch that breaks `git clone` for pip-from-git dependencies. Run `setup.sh` on the login node in tmux, not through `sbatch`.

**Pipeline fails at `CONFIGURE_PVACSEQ` with `Failed to create Conda environment`**
Nextflow couldn't find the pre-built conda envs in `nextflow_work/conda/`. Most often caused by deleting `nextflow_work/` between runs. Rerun `bash setup.sh` to rebuild and reset paths.

**Stage 7 fails with `EOFError: Compressed file ended before the end-of-stream marker`**
A stale, truncated VEP VCF from a previously killed run was picked up. The pipeline reads finalized VCFs from `results/nextflow_*/ensemblvep/` — if you still see this, delete the stale work dir under `nextflow_work/` and rerun.

**`module load lsfm-init-miniconda/1.0.0` fails on your cluster**
That module is ZHAW-specific. Setup.sh installs its own miniforge into `$BASE/miniforge3/` and uses that throughout — the `module load` lines in the slurm scripts only initialise the cluster's base conda before activating the project's own env. If your cluster uses a different module path, edit the `module load` lines in `pipeline/slurm/*.sh`.

**Inode quota exceeded (BeeGFS / Lustre)**
Conda envs are inode-heavy (~200k for the four envs combined). Check usage with `du --inodes -s ./*` and clear unused conda envs elsewhere if needed.

**MHC Class II algorithms fail with `not valid. Skipping.`**
This pipeline runs MHC Class I only. Class II HLA allele formats (HLA-DPA1, DPB1, DQA1, DQB1, DRB1) are rejected by pVACtools 5.3.1 and would require additional setup not covered here.
## Citation

If you use this pipeline, please cite:

- **pVACtools**: Hundal et al. (2020) *Cancer Immunology Research* — pVACtools: A Computational Toolkit to Identify and Visualize Cancer Neoantigens
- **VEP**: McLaren et al. (2016) *Genome Biology* — The Ensembl Variant Effect Predictor
- **NMD rules**: Lindeboom et al. (2019) *Nature Genetics* — The rules and impact of nonsense-mediated mRNA decay in human cancers

---

## License

MIT License — see LICENSE file.
