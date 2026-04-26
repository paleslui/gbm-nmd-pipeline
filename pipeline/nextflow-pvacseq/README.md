# nextflow_pvacseq

## Introduction

**nextflow_pvacseq** is a custom [Nextflow](https://www.nextflow.io/) pipeline that transforms MAF files into VCF, annotates them with [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), and analyzes them with [pVACseq](https://pvactools.readthedocs.io/en/latest/tools/pvacseq.html) to facilitate the investigation of tumor neoantigens.
It supports inputs in both MAF and VCF formats.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="img/nf_diagram3_dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="img/nf_diagram3.svg">
  <img alt="Project diagram" src="img/nf_diagram3_dark.png">
</picture>

## Pipeline Summary

The pipeline performs the following steps:

1. **Input Preprocessing**

    - Accepts MAF or VCF files as input.

    - Converts MAF to VCF (if required) using [maf2vcf](https://github.com/mskcc/vcf2maf/tree/main).

2. **Variant Annotation**

    - Annotates variants using [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), configured for pVACseq requirements.

3. **Loading HLA**

    - Reads and processes HLA typing information from a user-provided CSV file.

4. **pVACseq Setup**

    - Configures and downloads MHC class I and II reference files required by pVACseq if not provided.

5. **pVACseq Execution**

    - Runs [pVACseq](https://pvactools.readthedocs.io/en/latest/tools/pvacseq.html) to predict tumor neoantigens.

6. **MultiQC**

    - Aggregates results with [MultiQC](http://multiqc.info/).



## Usage

### Required Inputs

1. **Input Files**

    - A directory containing `*.maf` or `*.vcf` files.

2. **HLA Typing Information**

    - A CSV file (`--hla_csv`) in the format:

        ```
        Sample_ID,HLA_Types
        TCGA-G4-6310-01A-11D-1719-10,HLA-C05:01;HLA-C06:02;HLA-B45:01;HLA-A29:02;HLA-B44:02;HLA-A02:01
        ```

3. **Reference Genome**

    - A FASTA file (`--fasta`).

4. **VEP Requirements** (Optional)

    - Pre-downloaded **VEP cache** (`--vep_cache`) and/or **VEP plugins** (`--vep_plugins`).

    - If not provided, the pipeline will download the required files automatically.

5. **pVACseq Requirements** (Optional)

    - Pre-installed IEDB directory (`--pvacseq_iedb`).

    - If not provided, the pipeline will download and configure IEDB automatically.


### Running the Pipeline

```bash
nextflow run main.nf \
   -profile <docker|conda> \
   --input <INPUT DIRECTORY> \
   --hla_csv <HLA CSV FILE> \
   --fasta <REFERENCE FASTA> \
   --outdir <OUTPUT DIRECTORY>
```

### Testing the Pipeline

A minimal test dataset is included with the pipeline to verify that installation and execution work correctly.
The test profile uses **online VEP annotation**, so no local VEP cache is required.

⚠️ This mode is intended **only for testing**. It is **not recommended** for real analyses, as online VEP annotation is slower and less reproducible.

```bash
nextflow run main.nf -profile test,<docker|conda>
```

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

- [pVACseq](https://pvactools.readthedocs.io/en/latest/index.html)

> Jasreet Hundal+, Susanna Kiwala+, Joshua McMichael, Christopher A Miller, Alexander T Wollam, Huiming Xia, Connor J Liu, Sidi Zhao, Yang-Yang Feng, Aaron P Graubert, Amber Z Wollam, Jonas Neichin, Megan Neveau, Jason Walker, William E Gillanders, Elaine R Mardis, Obi L Griffith, Malachi Griffith. pVACtools: a computational toolkit to select and visualize cancer neoantigens. Cancer Immunology Research. 2020 Mar;8(3):409-420. DOI: 10.1158/2326-6066.CIR-19-0401. PMID: 31907209. (+) equal contribution.

> Jasreet Hundal, Susanna Kiwala, Yang-Yang Feng, Connor J. Liu, Ramaswamy Govindan, William C. Chapman, Ravindra Uppaluri, S. Joshua Swamidass, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. Accounting for proximal variants improves neoantigen prediction. Nature Genetics. 2018, DOI: 10.1038/s41588-018-0283-9. PMID: 30510237.

> Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. pVACseq: A genome-guided in silico approach to identifying tumor neoantigens. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: 26825632.

- [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)

> McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P Cunningham F. The Ensembl Variant Effect Predictor. Genome Biology Jun 6;17(1):122. (2016) doi:10.1186/s13059-016-0974-4

- [vcf2maf](https://github.com/mskcc/vcf2maf)

> Cyriac Kandoth. mskcc/vcf2maf: vcf2maf v1.6. (2020). doi:10.5281/zenodo.593251

- [nf-core template](https://nf-co.re/)

> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen. **The nf-core framework for community-curated bioinformatics pipelines.** _Nature Biotechnology._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
