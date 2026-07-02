# GBM NMD-Neoantigen Pipeline

A computational pipeline for identifying NMD-sensitive immunogenic mutations in primary and recurrent Glioblastoma (GBM). Developed as part of a Master's thesis at ZHAW. While the focus is GBM NMD/neoantigen biology, the paired primary/recurrent design and the five analytical modules generalise to any paired-cohort study - any two matched timepoints or conditions with per-patient VCFs.

**Thesis:** Identifying NMD-Sensitive Immunogenic Mutations in Primary and Recurrent Glioblastoma  
**Author:** Luigi Palese  
**Supervisors:** Maria Anisimova, Tugce Bilgin Sonay, Erik Vassella, Massimo Maiolo

---

## Background

Temozolomide (TMZ) chemotherapy induces frameshift mutations in GBM tumor cells that introduce premature termination codons (PTCs). Nonsense-mediated mRNA decay (NMD) degrades these transcripts before they can be translated and presented to the immune system. This pipeline identifies which TMZ-induced mutations produce immunogenic neoantigens that are being silenced by NMD - candidates that could become visible to the immune system upon NMD inhibition.

---

## Pipeline Overview

The pipeline runs in three stages, all submitted by a single master SLURM script with `afterok` dependencies:

| Stage | Script | Function |
|-------|--------|----------|
| **1** | `gbm_analysis.py` | Mutation landscape: VCF parsing, mutation burden, paired T/M comparison, gene-level recurrence, HLA typing summary, SBS11 (temozolomide) mutational-signature analysis |
| **2** | `nextflow-pvacseq` | NMD-relevant variant filter → VEP v113 annotation (NMD plugin) + pVACseq neoantigen prediction (MHC Class I, 8 algorithms) |
| **3** | `nmd_scoring.py` + `nmd_cohort_summary.py` | NMD sensitivity scoring per sample (ensemble of VEP NMD plugin + Lindeboom rules) and cohort-level aggregate report |

A single `sbatch master_pipeline.sh` produces a complete, reproducible run with all three stages chained together.

---

## Analytical modules

The three stages implement five analytical modules. Modules 1–2 run in Stage 1 (`gbm_analysis.py`); Modules 3–5 run in Stage 3 (`nmd_scoring.py` per sample, aggregated by `nmd_cohort_summary.py`). Each module emits a per-sample summary table and a headline plot.

| # | Module | What it tests | Per-sample TSV | Headline plot |
|---|--------|---------------|----------------|---------------|
| **1** | Mutation burden analysis (Stage 1) | Whether tumour mutation burden (TMB) shifts between primary (T) and recurrent (M) timepoints across paired patients. | `1_gbm_analysis/sample_mutation_summary.tsv` | `plot_paired_tmb.png` |
| **2** | Truncating mutation enrichment (Stage 1) | Whether the fraction of truncating variants (frameshift / stop-gained / splice / start-lost) is enriched at recurrence - the PTC-generating events NMD acts on. | `1_gbm_analysis/sample_mutation_summary.tsv` | `plot_paired_truncating_fraction.png` |
| **3** | NMD-sensitive mutation analysis (Stage 3) | Whether the count of NMD-sensitive truncating variants (VEP NMD plugin + Lindeboom ensemble) changes between T and M. | `3_nmd_analysis/cohort/per_sample_nmd_summary.tsv` | `m3_nmd_sensitive_paired.png` |
| **4** | Neoepitope burden analysis (Stage 3) | Whether the per-sample neoepitope (neoantigen) burden from pVACseq differs between timepoints. | `3_nmd_analysis/cohort/per_sample_neoepitope_summary.tsv` | `m4_neoepitope_burden_paired.png` |
| **5** | Neoepitopes stratified by NMD status (Stage 3) | Whether neoepitopes arising from NMD-sensitive vs NMD-insensitive transcripts diverge between T and M - the core thesis hypothesis (epitopes silenced by NMD). | `3_nmd_analysis/cohort/per_sample_neoepitope_nmd_summary.tsv` | `m5_neoepitope_by_nmd_paired.png` |

---

## Requirements

- SLURM HPC cluster (Linux, RHEL8 or compatible)
- ~50 GB scratch space for conda environments + VEP cache + IEDB tools
- ~200 k inodes available (conda envs are inode-heavy)
- Internet access from the compute node running `setup.sh` (for miniforge, VEP cache, IEDB, mhcflurry, VEP plugins)
- Your own WES VCF files and HLA typing data (see Input Data below)

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/paleslui/gbm-nmd-pipeline.git
cd gbm-nmd-pipeline
```

### 2. Edit `config.sh`

Open `config.sh` and change **only these settings**:

```bash
# Absolute path to where you cloned this repository
BASE=/path/to/gbm-nmd-pipeline

# Your HPC SLURM partition + constraint + email
SLURM_PARTITION=earth-3
SLURM_CONSTRAINT=rhel8
SLURM_MAIL=you@example.com
```

Everything else is derived automatically.

### 3. Add your data

```
data/
├── vcf/                     ← your VCF files (see format below)
├── hla_typing/              ← HLA typing files (Dragen TSVs, optional)
├── hla_typing_classI.csv    ← HLA Class I alleles CSV (required)
└── reference/
    └── GRCh38.primary_assembly.fa   ← from Ensembl
```

### 4. Run setup (one time only, ~1h 20m)

Setup downloads miniforge, builds 5 conda envs, downloads the VEP v113 cache, downloads MHCflurry models, populates VEP plugins, and patches the Nextflow pipeline. It is submitted as a SLURM job that runs on a compute node:

```bash
sbatch pipeline/slurm/slurm_setup.sh
```

Monitor:
```bash
tail -f logs/Slurm-<job_id>.out
```

You only need to run this once per fresh clone. Setup is idempotent - re-running it on a partial install fixes anything missing without redoing existing steps.

### 5. Run the pipeline

```bash
sbatch pipeline/slurm/master_pipeline.sh
```

Three SLURM jobs are submitted with `afterok` dependencies (Stage 1 → Stage 2 → Stage 3). Total wall time: ~2.5 h on 56 paired GBM samples.

Monitor progress:
```bash
watch -n 30 'squeue -u $USER'
```

When complete, open the cohort report:
```
results/run_<TIMESTAMP>/3_nmd_analysis/cohort/cohort_report.html
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

VEP-annotated VCFs from Dragen with `--canonical --mane --plugin NMD` are recommended for best results. SnpEff-annotated input is also accepted.

### HLA typing CSV

`data/hla_typing_classI.csv` must contain Class I HLA alleles in this format:

```csv
Sample_ID,HLA_Types
11_T,HLA-A*26:01;HLA-A*01:01;HLA-B*08:01;HLA-B*38:01;HLA-C*07:01;HLA-C*12:03
45_T,HLA-A*02:01;HLA-A*03:01;HLA-B*07:02;HLA-B*15:01;HLA-C*03:04;HLA-C*07:02
...
```

Note: Sample IDs use underscore (`11_T`); VCF filenames do not (`11T`).
Only samples present in this CSV will be processed by pVACseq.

---

## Output Structure

Every run produces a single timestamped directory with three sub-folders, one per stage:

```
results/
└── run_<TIMESTAMP>/
    ├── 1_gbm_analysis/                          ← Stage 1
    │   ├── stage1_report.html                   ← interactive HTML report
    │   ├── sample_mutation_summary.tsv          ← per-sample burden, truncating fraction, SBS11
    │   ├── paired_variant_overlap.tsv           ← T/M shared-vs-private counts (extended)
    │   ├── paired_stats.tsv                     ← Module 1+2 paired Wilcoxon + BH-FDR
    │   ├── all_truncating_variants.tsv          ← every truncating variant (FS/SG/splice/start-lost)
    │   ├── gene_recurrence.tsv
    │   ├── hla_typing_summary.tsv
    │   └── plot_*.png                           ← 8 headline PNGs (300 dpi)
    │
    ├── 2_pvacseq/                               ← Stage 2 (Nextflow)
    │   ├── ensemblvep/                          ← VEP-annotated VCFs per sample
    │   ├── pvactools/<sample>/MHC_Class_I/      ← pVACseq output per sample
    │   │   └── <sample>.filtered.tsv            ← neoantigen candidates
    │   ├── multiqc/                             ← MultiQC report
    │   └── pipeline_info/                       ← Nextflow execution metadata
    │
    └── 3_nmd_analysis/                          ← Stage 3
        ├── per_sample/<sample>/                 ← per-sample NMD reports (one dir per sample)
        │   ├── nmd_scored_candidates.tsv
        │   ├── nmd_hla_breakdown.tsv
        │   └── report_nmd.html                  ← per-sample HTML
        └── cohort/                              ← cohort-level summary
            ├── cohort_candidates.tsv            ← all candidates (101 cols; adds binder_class + clonality)
            ├── cohort_summary.tsv
            ├── cohort_paired.tsv                ← T (primary) vs M (recurrent) per patient
            ├── cohort_tier1.tsv
            ├── per_sample_nmd_summary.tsv             ← Module 3 per-sample table
            ├── per_sample_neoepitope_summary.tsv      ← Module 4 per-sample table
            ├── per_sample_neoepitope_nmd_summary.tsv  ← Module 5 per-sample table
            ├── paired_stats_stage3.tsv                ← Module 3/4/5 paired Wilcoxon + BH-FDR
            ├── m3_*.png m4_*.png m5_*.png             ← 12 module plots (300 dpi)
            └── cohort_report.html               ← MAIN cohort HTML (~2.8 MB: 5 modules +
                                                   per-pair overlap widget + per-patient drill-down)
```

The cohort report (`3_nmd_analysis/cohort/cohort_report.html`) is the primary deliverable: a single HTML with overview cards, tier distributions, T-vs-M paired comparison per patient, top recurrent genes, HLA breakdown, and ranked TIER1 candidates across the cohort.

### Per-sample summary tables

The refactored pipeline writes one row per sample in each of these tables, the inputs to the paired statistics:

```
results/run_<TS>/1_gbm_analysis/sample_mutation_summary.tsv          # per-sample burden + truncating + SBS11
results/run_<TS>/1_gbm_analysis/paired_variant_overlap.tsv           # T/M shared-vs-private counts
results/run_<TS>/1_gbm_analysis/paired_stats.tsv                     # Module 1+2 paired Wilcoxon + BH-FDR
results/run_<TS>/1_gbm_analysis/all_truncating_variants.tsv          # every truncating variant (FS/SG/splice/start-lost)
results/run_<TS>/3_nmd_analysis/cohort/per_sample_nmd_summary.tsv             # NEW - Module 3
results/run_<TS>/3_nmd_analysis/cohort/per_sample_neoepitope_summary.tsv      # NEW - Module 4
results/run_<TS>/3_nmd_analysis/cohort/per_sample_neoepitope_nmd_summary.tsv  # NEW - Module 5
results/run_<TS>/3_nmd_analysis/cohort/paired_stats_stage3.tsv                # NEW - Module 3/4/5 Wilcoxon
```

The temozolomide (SBS11) mutational signature is reported as the `sbs11_count` / `sbs11_pct` columns of `sample_mutation_summary.tsv` (SBS11 = C>T at CpC, the COSMIC temozolomide signature).

### Statistical methodology

- **Paired test.** Every primary-vs-recurrent comparison uses the Wilcoxon signed-rank test (`scipy.stats.wilcoxon`) over patients with both timepoints.
- **Multiple-testing correction.** Benjamini–Hochberg FDR is applied *within each module family*: Modules 1+2 form one family in Stage 1 (`paired_stats.tsv`), and Modules 3+4+5 form one family in Stage 3 (`paired_stats_stage3.tsv`).
- **Usable pairs.** Reported `n_pairs` reflects only usable pairs - patients missing a timepoint, or with an undefined denominator for a given metric, are dropped from that test.
- **Effect sizes.** Each test reports `median_delta` and `mean_delta` (recurrent − primary) alongside the raw and BH-adjusted p-values, so direction and magnitude are visible independent of significance.

### TMB normalisation

- `TMB_per_Mb` is computed as truncating/total variant counts divided by a callable-region denominator taken from `config.sh`'s `TMB_DENOMINATOR_MB` (default **30 Mb**, assuming standard whole-exome capture).
- When a real capture-kit BED file is available, compute the empirical denominator:

  ```bash
  bedtools genomecov -bg -i target.bed | awk '{s += $3-$2} END {print s/1e6}'
  ```

- Set the result as `TMB_DENOMINATOR_MB` in `config.sh` and re-run Stage 1 (`sbatch pipeline/slurm/slurm_python.sh`, ~1 min) to refresh the TMB columns. The value is exported from `config.sh` so the Stage 1 Python child reads it from the environment.

---

## NMD-Relevant Variant Filter

Before pVACseq runs, raw VCFs are filtered to keep only variants whose VEP/SnpEff `Consequence` field includes one of the NMD-relevant categories:

- `frameshift_variant`
- `stop_gained`
- `stop_lost`
- `splice_donor_variant`, `splice_acceptor_variant`
- `splice_region_only` (with truncating downstream consequence)

Implemented in `pipeline/scripts/filter_nmd_relevant.sh` and integrated into `pipeline/slurm/slurm_pvacseq_filtered.sh`. The filter typically reduces variant counts ~28× (e.g. a sample with 1660 raw variants → 32 NMD-relevant variants), which dramatically speeds up Stage 2 (pVACseq runtime is dominated by IEDB binding-prediction calls).

The filter is idempotent - already-filtered VCFs in `data/vcf_filtered/` are skipped on re-runs.

---

## NMD Scoring

Each neoantigen candidate from pVACseq is classified by NMD sensitivity using an **ensemble of two methods**:

1. **VEP NMD plugin** (Ensembl) - reads the `NMD` field from the VEP CSQ annotation. An empty NMD field on a truncating variant means NMD is triggered (SENSITIVE). `NMD_escaping_variant` = INSENSITIVE.
2. **Lindeboom rules** ([Lindeboom et al. 2019, *Nat Genet*](https://www.nature.com/articles/s41588-019-0517-5)) - applied in priority order:

| Rule | Condition | Classification |
|------|-----------|----------------|
| Rule 4 | PTC within 150 nt of start codon | INSENSITIVE (pioneer round escape) |
| Rule 1 | PTC in last exon | INSENSITIVE (no downstream EJC) |
| Rule 3 | PTC in exon > 407 nt | INSENSITIVE (EJC density too low) |
| Rule 2 | PTC > 55 nt upstream of last EJC | SENSITIVE (canonical NMD) |

**Ensemble confidence score (0–3):** 3 = both methods agree, 2 = single method available, 1 = methods disagree, 0 = no data.

**Priority tiers:**

| Tier | Criteria | Interpretation |
|------|----------|----------------|
| **TIER1** | NMD-SENSITIVE + IC50 < 50 nM | Primary therapeutic targets - silenced by NMD, exposed by NMD inhibition |
| **TIER2** | NMD-SENSITIVE + IC50 50–500 nM | Moderate binders, potentially relevant after NMD inhibition |
| **TIER3** | NMD-INSENSITIVE + IC50 < 500 nM | Already expressed - controls for immune response without NMD inhibition |
| Unclassified | Missense variants or insufficient transcript info | NMD does not apply |

The cohort report breaks these out by timepoint (primary vs recurrent) and per patient, which is the key analysis for the thesis hypothesis.

---

## Resources

All resources are downloaded automatically by `setup.sh`. Total ~25 GB.

| Resource | Version | Location after setup |
|---|---|---|
| miniforge3 | latest | `miniforge3/` |
| Conda envs (orchestrator + Nextflow inner envs) | pinned via Nextflow conda specs | `conda_envs/`, `nextflow_work/conda/` |
| VEP cache | v113 GRCh38 | `resources/vep_cache/` |
| VEP plugins (Ensembl + pVACtools custom) | v4.0.7 | `resources/VEP_plugins/` |
| MHCflurry models | latest | `resources/mhcflurry/` |
| IEDB MHC-I binding tools | latest | `resources/iedb/mhc_i/` (downloaded by Nextflow on first pipeline run) |

---

## Configuration Reference

All pipeline settings live in `config.sh`. The only required edits are `BASE`, `SLURM_PARTITION`, `SLURM_CONSTRAINT`, and `SLURM_MAIL`. Other settings you may want to adjust:

| Setting | Default | Description |
|---|---|---|
| `VEP_CACHE_VERSION` | `113` | VEP cache version |
| `PVACSEQ_ALGORITHMS` | 8 Class I algorithms (MHCflurry, MHCnuggetsI, NetMHC, NetMHCpan, NetMHCpanEL, PickPocket, SMM, SMMPMBEC) | Prediction algorithms used by pVACseq (median across all 8 = the binding consensus filter) |

pVACseq's default `--tdna-vaf 0.25` is used for tumor-DNA VAF filtering. See the troubleshooting section for guidance on tuning this for low-purity samples.

---

## Reproducibility Testing

The pipeline has been validated by three escalating destructive reproducibility tests on the same input data:

| Level | Wipe scope | Setup runtime | Pipeline runtime | Result |
|-------|-----------|---------------|------------------|--------|
| 1 | Run state only (kept all envs and resources) | n/a (idempotent) | 2 h 11 m | 56/56 ✓, 0 retries |
| 2 | Conda envs + Nextflow inner envs | ~30 m | 2 h 14 m | 56/56 ✓, 1 transient retry auto-recovered |
| 3 | Full fresh-clone (miniforge3, all envs, all resources) | 1 h 18 m | 2 h 23 m | 56/56 ✓, 0 retries |

All three runs produced **bit-for-bit identical** results across all 56 samples (338 total epitopes), with identical gold-standard hits in 11_T (ITGA4 frameshift YCIKLIHIV at HLA-C*12:03, BRAT1 frameshift STMSFCGTL at HLA-C*12:03).

**Integrity gate.** The reproducibility guarantee is enforced, not assumed. `setup.sh` STEP 10 is a hard gate: after building the conda environments it imports every critical package (scipy, Bio, pandas, pvactools, ...) and verifies the pVACtools `sys.executable` patch, exiting non-zero and naming the offending package if anything is missing. This closes the failure mode behind the earlier `20260630` production run, where a silent environment breakage (`No module named 'scipy'` in Stage 1; `No module named 'Bio'` in Stage 2) was reported as a successful setup. With this gate in place a fresh `git clone` → `setup.sh` → `master_pipeline.sh` reliably reproduces the byte-equivalent baseline, and a corrupted environment fails loudly at setup time rather than silently mid-run.

---

## Troubleshooting

> See [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) for the reproducibility failure modes (env integrity gate, python hijacking, the env-4a82 MHCnuggets patch, empty-sample `KeyError`, disk quota, `__BASE__` placeholders) and how to diagnose them.

**Setup fails before any conda activation with a git/openssl error**
This is the libk5crypto/OPENSSL_1_1_1b symbol conflict that hits system git on RHEL8 compute nodes when `LD_LIBRARY_PATH` is inherited from a conda env. `setup.sh` STEP 1c pre-populates `resources/VEP_plugins/` with `LD_LIBRARY_PATH` unset before any conda activation, so this error should not appear during pipeline runs. If you see it from a manual `git clone`, run it from a clean shell (no conda env active).

**Stage 2 fails with `Failed to create Conda environment`**
Nextflow couldn't find or build envs in `nextflow_work/conda/`. Most often caused by deleting `nextflow_work/` between runs or a partial setup. Re-run `sbatch pipeline/slurm/slurm_setup.sh` - it is idempotent and rebuilds anything missing.

**Stage 3 fails with `EOFError: Compressed file ended before the end-of-stream marker`**
A stale, truncated VEP VCF from a previously killed run was picked up. Delete the offending work dir under `nextflow_work/` and re-run.

**A sample produces 0 NMD candidates**
This usually means pVACseq's coverage filter (`--tdna-vaf 0.25` by default) rejected all of the sample's variants - common in low-tumor-purity samples. The variants exist in `2_pvacseq/pvactools/<sample>/MHC_Class_I/<sample>.all_epitopes.tsv` but none pass the strict default. To recover them, lower `--tdna-vaf` (e.g. to 0.10). The 8-algorithm binding consensus (median IC50 < 500 nM) protects against false positives, so a relaxed VAF still produces high-confidence neoantigens.

**`module load lsfm-init-miniconda/1.0.0` fails on your cluster**
That module is ZHAW-specific. `setup.sh` installs its own miniforge into `$BASE/miniforge3/` and uses that throughout - the `module load` lines in the slurm scripts only initialise the cluster's base conda before activating the project's own env. If your cluster uses different module paths, edit the `module load` lines in `pipeline/slurm/*.sh`.

**Inode quota exceeded (BeeGFS / Lustre)**
Conda envs are inode-heavy (~200 k for the four envs combined). Check usage with `du --inodes -s ./*` and clear unused conda envs elsewhere if needed.

**MHC Class II algorithms fail with `not valid. Skipping.`**
This pipeline runs MHC Class I only. Class II HLA allele formats (HLA-DPA1, DPB1, DQA1, DQB1, DRB1) are rejected by pVACtools 5.3.1 and would require additional setup not covered here.

---

## Citation

If you use this pipeline, please cite:

- **pVACtools**: Hundal et al. (2020) *Cancer Immunology Research* - pVACtools: A Computational Toolkit to Identify and Visualize Cancer Neoantigens
- **VEP**: McLaren et al. (2016) *Genome Biology* - The Ensembl Variant Effect Predictor
- **VEP NMD plugin**: Coelho et al. (2020) *NAR Genomics & Bioinformatics* - Nonsense-mediated mRNA decay prediction in the Ensembl Variant Effect Predictor
- **NMD rules**: Lindeboom et al. (2019) *Nature Genetics* - The rules and impact of nonsense-mediated mRNA decay in human cancers
- **Mutational signatures (SBS11)**: Alexandrov et al. (2020) *Nature* - The repertoire of mutational signatures in human cancer (COSMIC SBS11 = temozolomide exposure)
- **Paired statistics**: Virtanen et al. (2020) *Nature Methods* - SciPy 1.0 (`scipy.stats.wilcoxon`, Wilcoxon signed-rank test for all paired primary-vs-recurrent comparisons)
- **Multiple-testing correction**: Benjamini & Hochberg (1995) *J. R. Stat. Soc. B* - Controlling the false discovery rate (BH-FDR applied within each module family)

---

## License

MIT License - see LICENSE file.
