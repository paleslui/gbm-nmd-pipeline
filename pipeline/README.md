# Pipeline Scripts — Technical Reference

This directory contains all analysis scripts and the Nextflow pipeline for
VEP annotation and pVACseq neoantigen prediction.

---

## gbm_analysis.py — Stages 1-4

Parses all SnpEff-annotated VCFs, computes mutation burden, paired T/M
comparison, gene-level recurrence, TMZ signature enrichment, and NMD candidate
prioritisation. Generates a self-contained HTML report.

**Usage:**
```bash
python gbm_analysis.py \
    --vcf_dir  <path/to/vcf/> \
    --out_dir  <path/to/output/> \
    --hla_dir  <path/to/hla/>    # optional
    --fasta    <path/to/GRCh38.fa>  # optional, for TMZ signature
```

**Inputs:**
| Argument | Required | Description |
|---|---|---|
| `--vcf_dir` | Yes | Directory of `{id}{T\|M}-ensemble-annotated.vcf` files |
| `--out_dir` | Yes | Output directory (timestamped subdirectory created) |
| `--hla_dir` | No | Directory containing HLA typing CSV or Dragen TSV files |
| `--fasta` | No | GRCh38 FASTA for TMZ signature (auto-detected from `data/reference/`) |

**Outputs (in `out_dir/run_{timestamp}/`):**
| File | Description |
|---|---|
| `report.html` | Self-contained interactive HTML report (10 sections) |
| `summary_mutation_burden.tsv` | Per-sample PASS variant counts |
| `all_fs_sg_variants.tsv` | All frameshift/stop-gained variants |
| `paired_variant_overlap.tsv` | T vs M variant overlap per patient |
| `gene_recurrence.tsv` | Gene-level FS/SG hit counts |
| `tmz_signature.tsv` | C>T@CpG fraction per recurrent sample |
| `hla_typing_summary.tsv` | HLA alleles per sample (if --hla_dir provided) |
| `plot_*.png` | Individual plot files |

**Report sections:**
1. Dataset overview
2. Total somatic mutation burden
3. High-impact truncating variants (FS/SG)
4. Paired primary vs recurrent comparison
5. SNV vs indel breakdown
6. Variant overlap — shared vs timepoint-specific
7. Recurrence-acquired FS/SG — TMZ candidate pool
8. Gene-level recurrence
9. TMZ mutational signature — SBS11 enrichment
10. NMD neoantigen candidate prioritisation

---

## nextflow-pvacseq/ — Stages 5-6

Nextflow pipeline for VEP v113 annotation (with NMD plugin) and pVACseq
MHC Class I neoantigen prediction. Submitted via `slurm/slurm_pvacseq.sh`.

**Key modules:**
- `modules/local/vep/` — VEP v113 annotation with NMD, Frameshift, Wildtype plugins
- `modules/local/configure_pvacseq/` — pVACseq setup and execution
- `modules/nf-core/multiqc/` — QC report

**Important patches applied by setup.sh:**
- `modules/local/vep/main.nf` — `unset PERL5LIB` before VEP call (Perl conflict fix)
- `slurm.config` — `ENSEMBLVEP_DOWNLOAD` conda path (OOM fix)
- `pvactools/lib/prediction_class.py` — `sys.executable` for MHCnuggetsI (Python path fix)

**Outputs (in `results/nextflow_{RUN_TS}/`):**
| Path | Description |
|---|---|
| `ensemblvep/` | VEP-annotated VCFs per sample |
| `pvactools/{sample}/MHC_Class_I/{sample}.filtered.tsv` | pVACseq neoantigen candidates |
| `multiqc/` | MultiQC QC report |

---

## nmd_scoring.py — Stage 7

Scores each pVACseq candidate for NMD sensitivity using an ensemble of
the VEP NMD plugin and Lindeboom et al. 2019 rule-based method.

**Usage:**
```bash
python nmd_scoring.py \
    --pvacseq_tsv <path/to/11_T.filtered.tsv> \
    --vep_vcf     <path/to/11T-annotated_vep.vcf.gz> \  # optional
    --out_dir     <path/to/output/>
```

**Inputs:**
| Argument | Required | Description |
|---|---|---|
| `--pvacseq_tsv` | Yes | Filtered TSV from pVACseq |
| `--vep_vcf` | No | VEP-annotated VCF (gzipped). Falls back to rule-based only if absent |
| `--out_dir` | Yes | Output directory |

**NMD scoring methods:**

*Method 1 — VEP NMD plugin:*
Reads the `NMD` field from the VEP CSQ annotation. Empty field on a truncating
variant = NMD-SENSITIVE. `NMD_escaping_variant` = NMD-INSENSITIVE.

*Method 2 — Lindeboom rules (applied in priority order):*
| Rule | Condition | Classification |
|---|---|---|
| Rule 4 | PTC within 150nt of start codon | INSENSITIVE |
| Rule 1 | PTC in last exon | INSENSITIVE |
| Rule 3 | PTC in exon >407nt | INSENSITIVE |
| Rule 2 | PTC >55nt upstream of last EJC | SENSITIVE |

*Ensemble confidence score (0-3):*
| Score | Meaning |
|---|---|
| 3 | Both methods agree — high confidence |
| 2 | Single method available — medium confidence |
| 1 | Methods disagree — flag for review |
| 0 | No transcript data — unclassifiable |

**Priority tiers:**
| Tier | Criteria |
|---|---|
| Tier 1 | NMD-SENSITIVE + IC50 < 50 nM — primary therapeutic targets |
| Tier 2 | NMD-SENSITIVE + IC50 50-500 nM — moderate binders |
| Tier 3 | NMD-INSENSITIVE + IC50 < 500 nM — controls |
| Unclassified | Missense variants or no transcript data |

**Outputs (in `out_dir/`):**
| File | Description |
|---|---|
| `report_nmd.html` | Self-contained HTML report with 4 plots |
| `nmd_scored_candidates.tsv` | All candidates with NMD scores and tiers |
| `nmd_hla_breakdown.tsv` | Tier 1/2 candidates summarised by HLA allele |

---

## slurm/ — SLURM Job Scripts

| Script | Description |
|---|---|
| `master_pipeline.sh` | Submits all 3 stages with `afterok` dependencies |
| `slurm_python.sh` | Stage 1-4: runs gbm_analysis.py |
| `slurm_pvacseq.sh` | Stage 5-6: runs Nextflow VEP + pVACseq |
| `slurm_nmd.sh` | Stage 7: runs nmd_scoring.py for all samples |