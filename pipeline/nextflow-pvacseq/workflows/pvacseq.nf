/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAF2VCF                } from '../modules/local/maf2vcf/main'
include { ENSEMBLVEP_VEP         } from '../modules/local/vep/main'
include { PVACTOOLS_PVACSEQ      } from '../modules/local/pvacseq/main'
include { CONFIGURE_PVACSEQ      } from '../modules/local/configure_pvacseq/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pvacseq_pipeline'
include { CONFIGURE_PVACSEQ_IEDB } from '../subworkflows/local/configure_pvacseq_iedb'
include { SETUP_VEP_ENVIRONMENT  } from '../subworkflows/local/setup_vep_env'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



workflow PVACSEQ_PIPELINE {

    take:
    ch_maf_files   // Channel: directory with input files (MAF)
    ch_vcf_files   // Channel: directory with input files (VCF)
    fasta          // Path to reference genome
    hla_csv        // Preprocessed HLA CSV file path

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Process MAF files (if provided)

    //
    // MODULE: Run MAF2VCF
    //
    MAF2VCF (
        ch_maf_files,
        fasta
    )

    ch_vcf_files = ch_vcf_files.mix(MAF2VCF.out.vcf)
    ch_versions = ch_versions.mix(MAF2VCF.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_vcf_files.map { tuple ->
            return tuple[1] // Pass only path to VCF
        }
    )

    // Check and Install ENSEMBLVEP_VEP Parameters
    SETUP_VEP_ENVIRONMENT (
        params.vep_cache ?: [],
        params.vep_cache_version ?: '',
        params.vep_genome ?: '',
        params.vep_species ?: '',
        params.vep_plugins ?: []
    )

    //
    // MODULE: Run ENSEMBLVEP_VEP
    //
    ENSEMBLVEP_VEP (
        ch_vcf_files.map { tuple ->
            return [tuple[0], tuple[1], []] // Pass only metadata and path to VCF
        },
        SETUP_VEP_ENVIRONMENT.out.vep_genome,
        SETUP_VEP_ENVIRONMENT.out.vep_species,
        SETUP_VEP_ENVIRONMENT.out.vep_cache_version,
        SETUP_VEP_ENVIRONMENT.out.vep_cache,
        [ [id:"reference_genome"], fasta ],
        [],
        SETUP_VEP_ENVIRONMENT.out.vep_plugins
    )

    ch_multiqc_files = ch_multiqc_files.mix(ENSEMBLVEP_VEP.out.report)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())

    // Load and Normalize HLA Data
    hla_ch = Channel
        .fromPath(hla_csv)
        .splitCsv(header: true)
        .map { row ->
            // Extract sample ID and raw HLA string
            def sample_id = row.Sample_ID
            def hla_string = row.HLA_Types.replace(';', ',').trim()

            // Normalize HLA string to consistent format
            def normalized_hla = hla_string.split(',')
                .collect { hla ->
                    hla.trim().replaceAll(/HLA-([A-Z])([0-9]+:[0-9]+)/, 'HLA-$1*$2')
                }
                .join(',')
            return [["id": sample_id], normalized_hla]
        }


    tumor_pvacseq_ch = ENSEMBLVEP_VEP.out.vcf
        .map { tuple ->
            // Extract tumor and normal samples
            // Parse the VCF file to extract sample names
            def vcf_file = file(tuple[1])
            def header_line = readVcfChromLine(vcf_file)
            def columns = header_line.split('\t')
            def tumor_sample = columns[-2] // Second-to-last column is typically the tumor sample
            def normal_sample = columns[-1] // Last column is typically the normal sample

            // Extract tumor and normal sample from files and add them
            return [
                // meta
                [
                    "id": tumor_sample,
                    "normal_sample": normal_sample
                ],
                // sample_name
                tumor_sample,
                // vcf
                tuple[1]
            ]
        }
    // Reorder just so that normal sample is a key now
    normal_pvacseq_ch = tumor_pvacseq_ch
        .map { meta, sample_name, vcf ->
            return [
                // meta
                [
                    "id": meta.normal_sample,
                    "normal_sample": sample_name
                ],
                meta.normal_sample,
                // vcf
                vcf
            ]
        }

    // Merge HLA info in case we have tumor sample id
    tumor_pvacseq_ch_hla = tumor_pvacseq_ch
        .map{meta, sample_name, vcf -> [meta.subMap(["id"]), meta, sample_name, vcf]}
        .join(hla_ch)
        .map { id, meta, sample_name, vcf, hla ->
            return [
                meta,
                sample_name,
                hla,
                vcf
            ]
        }


    // Merge HLA info in case we have normal sample id
    normal_pvacseq_ch_hla = normal_pvacseq_ch
        .map{meta, sample_name, vcf -> [meta.subMap(["id"]), meta, sample_name, vcf]}
        .join(hla_ch)
        .map { id, meta, sample_name, vcf, hla ->
            return [
                meta,
                sample_name,
                hla,
                vcf
            ]
        }

    pvacseq_ch = tumor_pvacseq_ch_hla.mix(normal_pvacseq_ch_hla)

    // Download mhc_i and mhc_ii iedb if required
    CONFIGURE_PVACSEQ_IEDB (
        params.pvacseq_iedb ?: [],
        params.pvacseq_algorithm ?: ''
    )

    // Configure pvacseq before running it
    CONFIGURE_PVACSEQ (
        CONFIGURE_PVACSEQ_IEDB.out.iedb_dir,
        CONFIGURE_PVACSEQ_IEDB.out.iedb_mhc_i,
        CONFIGURE_PVACSEQ_IEDB.out.iedb_mhc_ii
    )

    config_done_ch = CONFIGURE_PVACSEQ.out.config_file

    // Ensure PVACSEQ starts only after CONFIGURE_PVACSEQ has completed.
    // CONFIGURE_PVACSEQ emits 'env_config_done.txt' as a completion signal.
    // We combine pvacseq_ch with this signal so that PVACSEQ will not run prematurely,
    // preventing race conditions where the IEDB setup is not yet available.
    // The config_file is only used to enforce this dependency; its value is ignored.
    pvacseq_ready_ch = pvacseq_ch.map { tuple ->
        config_done_ch.first().get()
        return tuple
    }

    //
    // MODULE: Run pVAcseq tool
    //
    PVACTOOLS_PVACSEQ (
        pvacseq_ready_ch,
        params.pvacseq_algorithm,
        CONFIGURE_PVACSEQ.out.iedb_dir,
        params.blastp_path ?: [],
        params.genes_of_interest ?: [],
        params.peptide_fasta ?: [],
        params.ph_proximal_variants_vcf ?: []
    )

    ch_multiqc_files = ch_multiqc_files.mix(PVACTOOLS_PVACSEQ.out.mhc_i_filtered.map { meta,file -> [file]})
    ch_multiqc_files = ch_multiqc_files.mix(PVACTOOLS_PVACSEQ.out.mhc_ii_filtered.map { meta,file -> [file]})

    ch_versions = ch_versions.mix(PVACTOOLS_PVACSEQ.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    iedb_dir       = CONFIGURE_PVACSEQ_IEDB.out.iedb_dir
    mode           = CONFIGURE_PVACSEQ_IEDB.out.mode
}

// Get the #CHROM header line from a gzipped VCF file in Nextflow DSL2
def readVcfChromLine(Path path) {
    def line = null
    try {
        path.withInputStream { stream ->
            def reader = new BufferedReader(
                new InputStreamReader(new java.util.zip.GZIPInputStream(stream), 'ASCII')
            )
            // Use `find` to locate the #CHROM line
            line = reader.lines()
                         .takeWhile { it.startsWith('#') }
                         .findAll { it.startsWith('#CHROM') }
                         .first()

            if (!line) {
                throw new IllegalArgumentException("VCF file (${path}) does not contain a #CHROM header line.")
            }
        }
    } catch (Exception e) {
        log.warn "VCF file (${path}): Error while reading header"
        log.warn "${e.message}"
    }

    return line
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
