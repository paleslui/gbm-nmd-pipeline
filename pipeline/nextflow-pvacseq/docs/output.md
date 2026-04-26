# pvacseq: Output

## Introduction

This document describes the output produced by the pVACseq pipeline.

## Output overview

The pVACseq pipeline output is grouped into three categories:

- **pVACseq Analysis** – Neoantigen prediction results.

- **MultiQC** – Summary report of VEP statistics.

- **Pipeline Information** – Metadata and logs from the pipeline run.


---

### pVACseq Analysis

The pVACseq outputs follow the structure and file formats described in the [official pVACseq documentation](https://pvactools.readthedocs.io/en/latest/pvacseq/output_files.html).

- `MHC_I/`: Results for MHC class I predictions.

- `MHC_II/`: Results for MHC class II predictions.


Each directory may contain:

- **`.all_epitopes.tsv`** – Comprehensive list of all predicted epitopes with scoring details.

- **`.filtered.tsv`** – Filtered list of candidate neoantigens passing the thresholds.


---

### MultiQC

- `multiqc/`

    - `multiqc_report.html`: standalone HTML file with VEP summary statistics.

    - `multiqc_data/`: parsed statistics from VEP runs.

    - `multiqc_plots/`: static plots included in the report.


[MultiQC](http://multiqc.info/) is used to provide a summary of VEP annotation results across all samples.
The report is limited to **VEP statistics only** in this pipeline.

---

### Pipeline Information

- `pipeline_info/`

    - Nextflow reports: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt`, `pipeline_dag.dot`/`pipeline_dag.svg`.

    - Pipeline reports: `pipeline_report.html`, `pipeline_report.txt`, `software_versions.yml`.

    - Run parameters: `params.json`.


[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides reports on execution details such as launch commands, runtime, and resource usage. These outputs can be used to troubleshoot or document pipeline runs.
