# pvacseq: Usage

## Introduction

This pipeline is built using the **nf-core** template, providing a standardized structure. However, the pipeline is not published in the nf-core repository.

The pipeline is designed to run the **pVAC-seq** tool on multiple tumor samples in different environments.

## Input / Output Options

|Parameter|Description|Required|
|---|---|---|
|`--input`|Path to directory with VCF/MAF files. Both formats are supported; MAF will be converted.|yes|
|`--hla_csv`|CSV file with HLA types per sample (`Sample_ID`, `HLA_Types`).|yes|
|`--outdir`|Directory for pipeline results. |yes|



### Directory Structure

- The input directory can contain both VCF and MAF files.
- **VCF Files**: Must include sample genotype information (`GT` field) as described in the [pVACtools documentation](https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/gt.html).
- **MAF Files**: Will be transformed to VCF format before processing.

### Example Directory Structure

```plaintext
input_vcf_maf/
├── sample1.vcf
├── sample2.maf
├── sample3.vcf
```

### HLA Input

The HLA input file must be a comma-separated file (CSV) with the following columns:

| Column      | Description                                                                                                            |
|-------------|------------------------------------------------------------------------------------------------------------------------|
| `Sample_ID` | Unique sample identifier. This must match the sample name inside the input VCF or MAF file, **not** the file name.     |
| `HLA_Types` | Semicolon-separated list of HLA alleles in the format `HLA-[Gene][Allele group]:[Protein]`.                            |

- **Header Row**: The first row must contain column names (`Sample_ID`, `HLA_Types`).
- **Unique Samples**: Each row corresponds to one unique sample.

### HLA Format

HLA alleles must use the `HLA-[Gene][Allele group]:[Protein]` format. Examples:
- `HLA-A02:01`
- `HLA-B15:01`
- `HLA-C07:02`

Alleles in the `HLA_Types` column must be separated by semicolons `;`.

### Example HLA Input File

```csv title="hla_input.csv"
Sample_ID,HLA_Types
SAMPLE_1,HLA-C06:02;HLA-B45:01
SAMPLE_2,HLA-A29:02;HLA-B44:02;HLA-A02:01
SAMPLE_3,HLA-A11:01;HLA-B35:01;HLA-C04:01
```

## Reference Genome Options

|Parameter|Description|Required|
|---|---|---|
|`--fasta`|Path to reference genome FASTA file. Must be uncompressed and match the input variant coordinates.|yes|


### Requirements

- The FASTA file must be **unzipped**. Compressed versions (e.g., `.fa.gz`) are not currently supported.
- Ensure the FASTA file corresponds to the correct reference genome version used for the input VCF/MAF files.

## VEP Annotation

The pVACseq pipeline uses VEP (Variant Effect Predictor) for annotating input variants. Below are the parameters required for the tool and their behavior in the pipeline:

### VEP Options

|Parameter|Description|Notes|
|---|---|---|
|`--vep_plugins`|Path to installed VEP plugins (e.g. Wildtypes, Frameshift).|If not set, plugins are downloaded automatically.|
|`--vep_cache`|Path to local VEP cache directory.|If not set, cache can be auto-downloaded.|
|`--vep_cache_version`|Version of VEP cache.|Must match cache used. Default: `102`.|
|`--vep_genome`|Genome identifier used by VEP (e.g. `GRCh38`).|Must match cache used. Default: `GRCh38`.|
|`--vep_species`|Species name used by VEP (e.g. `homo_sapiens`).|Must match cache used. Default: `homo_sapiens`.|
|`--extra_vep_args`|Extra arguments to pass to VEP.|Optional.|

### Automatic download

#### **`vep_plugins`**
- If `vep_plugins` is not provided, the pipeline will download the required plugins automatically.
- Specify the downloaded plugin directory in subsequent runs using the `vep_plugins` parameter.

#### Notes

> - Always specify the `vep_cache` and `vep_plugins` directories after the first run to avoid unnecessary downloads.
> - Ensure the `vep_cache_version`, `vep_genome`, `vep_species` matche the version of the `vep_cache` directory provided.


## pVACseq Parameters

The pVACseq pipeline provides a range of configurable options for neoantigen prediction. Below are the key parameters and their behavior:

## pVACseq Options

|Parameter|Description|Required|
|---|---|---|
|`--pvacseq_algorithm`|Epitope prediction algorithms to use (e.g. `NetMHCpan`, `MHCflurry`, etc.). Multiple indicated with "," delimiter.|yes|
|`--pvacseq_iedb`|Path to local IEDB installation. If not set, required parts are downloaded automatically.|no|
|`--blastp_path`|Path to BLASTP binary, if reference proteome similarity is used.|no|
|`--genes_of_interest`|File listing genes of interest (one per line)|no|
|`--peptide_fasta`|Custom peptide FASTA for similarity searches instead of BLASTP.|no|
|`--ph_proximal_variants_vcf`|VCF with phased proximal variants (gzipped + tabix indexed).|no|
|`--extra_pvacseq_args`|Extra arguments to pass to `pvacseq`.|no|

#### **`pvacseq_iedb`**

- Path to the IEDB installation directory.

- Behavior depends on what you provide:

    - **If not provided or folder is empty**:
        The pipeline will automatically download the required IEDB components (`mhc_i`, `mhc_ii`, or both) based on the algorithms selected.

    - **If provided and the folder contains files**:
        The pipeline will assume it already contains a valid IEDB installation with `mhc_i` and/or `mhc_ii` subdirectories.

⚠️ **Important limitation**:
The full path to the IEDB directory **must be shorter than 57 characters**.

- If the path exceeds this limit, the pipeline will try to create a hard link in a shorter temporary path.

- If hard linking is not possible, the files will be copied instead.

- After pipeline execution, the temporary IEDB directory is automatically removed.

## Running the Pipeline

To run the pipeline, provide all required parameters in a configuration file and execute the pipeline using the following command:

```bash
nextflow run main.nf -profile <conda/docker>
```

### Test Profile

A `test` profile is available for running the pipeline with test data. The test dataset includes a MAF file derived from the TCGA dataset of a human tumor.

```bash
nextflow run main.nf -profile test,<conda/docker>
```
