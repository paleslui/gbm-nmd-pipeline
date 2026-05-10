# Pipeline Scripts — Technical Reference

This directory contains all analysis scripts and the Nextflow pipeline for VEP annotation and pVACseq neoantigen prediction. The flow corresponds to the three pipeline stages:

```
Stage 1 → gbm_analysis.py
Stage 2 → scripts/filter_nmd_relevant.sh + nextflow-pvacseq/
Stage 3 → nmd_scoring.py + nmd_cohort_summary.py
```

All three stages are submitted by `slurm/master_pipeline.sh`. Each writes into `$RESULTS/run_<TIMESTAMP>/{1_gbm_analysis,2_pvacseq,3_nmd_analysis}/`.

---

## gbm_analysis.py — Stage 1

Parses all SnpEff- or VEP-annotated VCFs, computes mutation burden, paired T/M comparison, gene-level recurrence, TMZ signature enrichment, and HLA typing summary. Generates a self-contained HTML report.

**Usage (standalone):**
```bash
python gbm_analysis.py \
    --vcf_dir  <path/to/vcf/> \
    --out_dir  <path/to/output/> \
    --hla_dir  <path/to/hla/>          # optional
    --fasta    <path/to/GRCh38.fa>     # optional, for TMZ signature
```

**Inputs:**

| Argument | Required | Description |
|---|---|---|
| `--vcf_dir` | Yes | Directory of `{id}{T\|M}-ensemble-annotated.vcf` files |
| `--out_dir` | Yes | Output directory |
| `--hla_dir` | No | Directory containing HLA typing CSV or Dragen TSV files |
| `--fasta` | No | GRCh38 FASTA for TMZ signature (auto-detected from `data/reference/`) |

**Outputs (in `$RUN_DIR/1_gbm_analysis/`):**

| File | Description |
|---|---|
| `report.html` | Self-contained interactive HTML report (10 sections) |
| `summary_mutation_burden.tsv` | Per-sample PASS variant counts |
| `all_fs_sg_variants.tsv` | All frameshift / stop-gained variants |
| `paired_variant_overlap.tsv` | T vs M variant overlap per patient |
| `gene_recurrence.tsv` | Gene-level FS/SG hit counts |
| `tmz_signature.tsv` | C>T@CpG fraction per recurrent sample |
| `hla_typing_summary.tsv` | HLA alleles per sample (if `--hla_dir` provided) |
| `plot_*.png` | Individual plot files |

**Report sections:**
1. Dataset overview
2. Total somatic mutation burden
3. High-impact truncating variants (FS / SG)
4. Paired primary vs recurrent comparison
5. SNV vs indel breakdown
6. Variant overlap — shared vs timepoint-specific
7. Recurrence-acquired FS/SG — TMZ candidate pool
8. Gene-level recurrence
9. TMZ mutational signature — SBS11 enrichment
10. NMD neoantigen candidate prioritisation

---

## scripts/filter_nmd_relevant.sh — Stage 2 pre-filter

Filters raw VEP/SnpEff-annotated VCFs to NMD-relevant consequences only. Used to reduce variant counts ~28× before running pVACseq, which dramatically speeds up Stage 2.

**Usage (standalone):**
```bash
bash scripts/filter_nmd_relevant.sh <input.vcf> <output.vcf.gz>
```

**Kept consequences:** `frameshift_variant`, `stop_gained`, `stop_lost`, `splice_donor_variant`, `splice_acceptor_variant`, `splice_region_only` (with truncating downstream consequence).

**Note:** Requires `bgzip` to write the gzip-indexed output. Used inside `slurm_pvacseq_filtered.sh` from within the `nf_pvacseq` conda env (which has `bgzip` via `samtools`). When run standalone, ensure `bgzip` is on the PATH.

The integrated wrapper in `slurm_pvacseq_filtered.sh` is idempotent — already-filtered VCFs in `data/vcf_filtered/` are skipped.

---

## nextflow-pvacseq/ — Stage 2

Nextflow pipeline for VEP v113 annotation (with NMD plugin) and pVACseq MHC Class I neoantigen prediction. Submitted via `slurm/slurm_pvacseq_filtered.sh`.

**Key modules:**

- `subworkflows/local/setup_vep_env/` — VEP environment setup (cache + plugins)
- `modules/local/vep/` — VEP v113 annotation with NMD, Frameshift, Wildtype plugins
- `modules/local/pvacseq/` — patched to copy IEDB to node-local SSD before invoking pvacseq run (avoids BeeGFS file-open contention under concurrent pVACseq tasks)
- `subworkflows/local/configure_pvacseq_iedb/` — pVACseq setup
- `modules/nf-core/multiqc/` — QC report

**Patches applied by `setup.sh`:**

- STEP 1c — Pre-populates `resources/VEP_plugins/` from GitHub + pVACtools archive before any conda activation. Avoids the libk5crypto/OPENSSL symbol-lookup error that breaks the Nextflow `DOWNLOAD_VEP_PLUGINS` task on RHEL8 compute nodes.
- STEP 6 — `modules/local/vep/main.nf` — `unset PERL5LIB` before VEP call (Perl conflict fix).
- STEP 7 — `slurm.config` — `ENSEMBLVEP_DOWNLOAD` conda path.
- STEP 9 — `modules/local/pvacseq/main.nf` — node-local IEDB copy patch (BeeGFS contention fix).

**Outputs (in `$RUN_DIR/2_pvacseq/`):**

| Path | Description |
|---|---|
| `ensemblvep/` | VEP-annotated VCFs per sample |
| `pvactools/<sample>/MHC_Class_I/<sample>.filtered.tsv` | pVACseq neoantigen candidates that pass binding + coverage filters |
| `pvactools/<sample>/MHC_Class_I/<sample>.all_epitopes.tsv` | All epitope predictions before filtering (useful for tuning thresholds) |
| `multiqc/` | MultiQC QC report |
| `pipeline_info/` | Nextflow execution trace, timeline, report |

---

## nmd_scoring.py — Stage 3 (per sample)

Scores each pVACseq candidate for NMD sensitivity using an ensemble of the VEP NMD plugin and Lindeboom et al. 2019 rule-based method.

**Usage (standalone):**
```bash
python nmd_scoring.py \
    --pvacseq_tsv <path/to/<sample>.filtered.tsv> \
    --vep_vcf     <path/to/<sample>-annotated_vep.vcf.gz> \   # optional
    --out_dir     <path/to/output/>
```

**Inputs:**

| Argument | Required | Description |
|---|---|---|
| `--pvacseq_tsv` | Yes | Filtered TSV from pVACseq |
| `--vep_vcf` | No | VEP-annotated VCF (gzipped). If absent, Lindeboom rule-based scoring is used alone |
| `--out_dir` | Yes | Output directory |

**NMD scoring methods:**

*Method 1 — VEP NMD plugin*  
Reads the `NMD` field from the VEP CSQ annotation. Empty field on a truncating variant = NMD-SENSITIVE. `NMD_escaping_variant` = NMD-INSENSITIVE.

*Method 2 — Lindeboom rules (priority order):*

| Rule | Condition | Classification |
|------|-----------|----------------|
| Rule 4 | PTC within 150 nt of start codon | INSENSITIVE |
| Rule 1 | PTC in last exon | INSENSITIVE |
| Rule 3 | PTC in exon > 407 nt | INSENSITIVE |
| Rule 2 | PTC > 55 nt upstream of last EJC | SENSITIVE |

*Ensemble confidence (0–3):*

| Score | Meaning |
|---|---|
| 3 | Both methods agree — high confidence |
| 2 | Single method available — medium |
| 1 | Methods disagree — flag for review |
| 0 | No transcript data — unclassifiable |

**Priority tiers:**

| Tier | Criteria |
|------|----------|
| TIER1 | NMD-SENSITIVE + IC50 < 50 nM — primary therapeutic targets |
| TIER2 | NMD-SENSITIVE + IC50 50–500 nM — moderate binders |
| TIER3 | NMD-INSENSITIVE + IC50 < 500 nM — controls |
| Unclassified | Missense variants or no transcript data |

**Outputs (in `$RUN_DIR/3_nmd_analysis/per_sample/<sample>/`):**

| File | Description |
|---|---|
| `report_nmd.html` | Self-contained per-sample HTML (4 plots + tables) |
| `nmd_scored_candidates.tsv` | All candidates with NMD scores and tiers |
| `nmd_hla_breakdown.tsv` | TIER1/2 candidates summarised by HLA allele |

---

## nmd_cohort_summary.py — Stage 3 (cohort-level)

Aggregates per-sample NMD outputs across the cohort. Imports plotting functions from `nmd_scoring.py` so the cohort report's visual style matches the per-sample reports exactly.

**Usage (standalone):**
```bash
python nmd_cohort_summary.py \
    --input_dir <path/to/3_nmd_analysis/per_sample/> \
    --out_dir   <path/to/3_nmd_analysis/cohort/>
```

**Inputs:**

| Argument | Required | Description |
|---|---|---|
| `--input_dir` | Yes | Directory with one subdirectory per sample, each containing `nmd_scored_candidates.tsv` |
| `--out_dir` | Yes | Where the cohort outputs are written |

**Outputs (in `$RUN_DIR/3_nmd_analysis/cohort/`):**

| File | Description |
|---|---|
| `cohort_report.html` | Single HTML report with 7 plots covering tier distribution (overall and split by timepoint), per-patient T-vs-M paired bar chart, IC50 distribution, top recurrent genes, HLA breakdown, NMD method comparison, and confidence distribution |
| `cohort_candidates.tsv` | All candidates concatenated, with `sample`/`patient`/`timepoint` columns prepended |
| `cohort_summary.tsv` | Cohort-level counts: total, by tier, by NMD class, by confidence |
| `cohort_paired.tsv` | Per-patient T (primary) vs M (recurrent) comparison: `T_total`, `M_total`, `T_tier1`, `M_tier1`, `delta_total`, `delta_tier1` |
| `cohort_tier1.tsv` | TIER1-only candidates across the cohort, ranked by IC50 |

The cohort report is the primary deliverable of a pipeline run — it surfaces the recurrent-vs-primary expansion of NMD-sensitive neoantigens that is the core thesis hypothesis.

---

## slurm/ — SLURM Job Scripts

| Script | Purpose |
|--------|---------|
| `slurm_setup.sh` | One-time setup: installs miniforge3, builds 5 conda envs, downloads VEP cache + MHCflurry, populates VEP plugins, applies pipeline patches. Submit once after cloning. |
| `master_pipeline.sh` | Orchestrator: submits all 3 stages with `afterok` dependencies, exports `RUN_TS` and `RUN_DIR` to children. This is the single command for an end-to-end run. |
| `slurm_python.sh` | Stage 1: runs `gbm_analysis.py` with `RUN_DIR/1_gbm_analysis/` as output |
| `slurm_pvacseq_filtered.sh` | Stage 2: integrates the NMD-relevant filter, then runs Nextflow with `RUN_DIR/2_pvacseq/` as output |
| `slurm_pvacseq.sh` | Legacy: Stage 2 without the pre-filter. Kept for reproducibility comparisons against the unfiltered path. Not used by `master_pipeline.sh` |
| `slurm_nmd.sh` | Stage 3: runs `nmd_scoring.py` per sample into `RUN_DIR/3_nmd_analysis/per_sample/<sample>/`, then runs `nmd_cohort_summary.py` into `RUN_DIR/3_nmd_analysis/cohort/` |

All slurm scripts can be submitted standalone (they auto-default `RUN_TS` and `RUN_DIR` to a fresh timestamped directory if not exported by a parent job).
