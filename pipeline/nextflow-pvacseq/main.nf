#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PVACSEQ_PIPELINE        } from './workflows/pvacseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_pvacseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_pvacseq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NF_PVACSEQ {

    take:
    maf_files // channel: directory with maf files read in from --input
    vcf_files // channel: directory with vcf files read in from --input
    main:

    //
    // WORKFLOW: Run pipeline
    //
    PVACSEQ_PIPELINE (
        maf_files,
        vcf_files,
        params.fasta,
        params.hla_csv
    )

    emit:
    multiqc_report = PVACSEQ_PIPELINE.out.multiqc_report // channel: /path/to/multiqc_report.html
    iedb_dir = PVACSEQ_PIPELINE.out.iedb_dir
    mode = PVACSEQ_PIPELINE.out.mode
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir
    )


    //
    // Generate input channels dynamically for MAF and VCF files
    //
    ch_maf_files = Channel
        .fromPath(params.input + "/*.maf")
        .map { file ->
            [ [id: file.baseName], file ] // Create tuples with metadata (id) and file path
        }

    ch_vcf_files = Channel
        .fromPath(params.input + "/*.vcf")
        .map { file ->
            [ [id: file.baseName], file ] // Create tuples with metadata (id) and file path
        }

    //
    // WORKFLOW: Run main workflow
    //
    NF_PVACSEQ (
        ch_maf_files,
        ch_vcf_files
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.monochrome_logs,
        NF_PVACSEQ.out.iedb_dir,
        NF_PVACSEQ.out.mode
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
