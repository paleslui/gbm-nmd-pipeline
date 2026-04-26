# Input Data — Format Guide

This directory contains all input data required by the pipeline.
Only `hla_typing_classI.csv` and the VCF files need to be provided by the user.
The reference genome and HLA typing files from Dragen are placed here manually.

---

## vcf/ — Somatic Variant Call Files

Place unzipped SnpEff-annotated VCF files here using this naming convention:

```
{patient_id}{T|M}-ensemble-annotated.vcf
```

- `T` = primary tumour (pre-treatment)
- `M` = recurrent tumour (post-TMZ)

**Example:**
```
vcf/
├── 11T-ensemble-annotated.vcf
├── 11M-ensemble-annotated.vcf
├── 45T-ensemble-annotated.vcf
├── 45M-ensemble-annotated.vcf
└── ...
```

**Requirements:**
- SnpEff annotation with `ANN` INFO field
- `PASS` filter tag on variant records
- Both T and M files required per patient for paired analysis
- Patients with only T or only M are parsed but excluded from paired comparison

**VEP-annotated VCFs (Dragen):**
VEP-annotated VCFs from Dragen with `--canonical --mane --plugin NMD` are
automatically detected via the `CSQ` INFO field and used for more accurate
NMD scoring in Stage 7. SnpEff VCFs (ANN field) are used for Stages 1-4.

---

## hla_typing_classI.csv — HLA Class I Alleles

Required for pVACseq (Stage 6). Contains MHC Class I alleles per sample.

**Format:**
```csv
Sample_ID,HLA_Types
11_T,HLA-A*26:01;HLA-A*01:01;HLA-B*08:01;HLA-B*38:01;HLA-C*07:01;HLA-C*12:03
45_T,HLA-A*02:01;HLA-A*03:01;HLA-B*07:02;HLA-B*15:01;HLA-C*03:04;HLA-C*07:02
```

**Notes:**
- `Sample_ID` uses underscore format (`11_T`), while VCF filenames do not (`11T`)
- Alleles separated by semicolons
- Only Class I alleles (HLA-A, HLA-B, HLA-C) — Class II requires separate typing
- Only samples present in this CSV will be processed by pVACseq

---

## hla_typing/ — Dragen HLA Typing Files (optional)

Raw Dragen HLA typing output for Stage 4 integration. If provided via
`--hla_dir`, `gbm_analysis.py` loads and summarises allele calls.

**Expected file format:**
```
S{patient}_{T|M}.hla.tsv         — allele calls from matched normal
```

**Alternative:** Provide a pre-formatted `hla_typing_classI.csv` directly
(see above) and skip `--hla_dir`.

---

## reference/ — Reference Genome

The GRCh38 reference genome is required for TMZ mutational signature analysis
(Stage 3D). It is placed here manually and indexed automatically by `setup.sh`.

**Required file:**
```
reference/
└── GRCh38.primary_assembly.fa          ← download from Ensembl
    GRCh38.primary_assembly.fa.fai      ← created automatically by setup.sh
```

**Download:**
```bash
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/\
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa \
   data/reference/GRCh38.primary_assembly.fa
```

The `.fai` index is created automatically by `setup.sh` using samtools.
If the FASTA is absent, the TMZ signature section in the report will show
"not available" but all other pipeline stages will run normally.