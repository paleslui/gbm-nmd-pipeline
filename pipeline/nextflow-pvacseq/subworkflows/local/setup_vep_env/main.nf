 include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'

//
// Process: Download VEP Plugins
//
process DOWNLOAD_VEP_PLUGINS {
    tag "$meta.id"
    label 'process_single'

    input:
    val(meta)

    output:
    tuple val(meta), path(prefix), emit: vep_plugins

    when:
    task.ext.when == null || task.ext.when

    script:
    def pvacseq_version = task.ext.args?['pvacseq-version'] ?: '4.0.7'
    prefix = task.ext.prefix ?: 'VEP_plugins'

    """
    mkdir -p "$prefix"
    git clone https://github.com/Ensembl/VEP_plugins.git "$prefix"
    wget "https://github.com/griffithlab/pVACtools/archive/refs/tags/v${pvacseq_version}.zip"
    unzip v${pvacseq_version}.zip "pVACtools-${pvacseq_version}/pvactools/tools/pvacseq/VEP_plugins/*"
    mv pVACtools-${pvacseq_version}/pvactools/tools/pvacseq/VEP_plugins/* "$prefix"
    rm -rf pVACtools-${pvacseq_version}
    rm v${pvacseq_version}.zip
    """
}


workflow SETUP_VEP_ENVIRONMENT {
    take:
    vep_cache
    vep_cache_version
    vep_genome
    vep_species
    vep_plugins


    main:
    // We cant have no version, but if vep_cache is provided we can't download the new one eather
    if (vep_cache && !vep_cache_version) {
        throw new Exception("VEP cache provided but version not indicated. Please provide 'vep_cache_version'.")
    }

    vep_cache_version = vep_cache_version ?: '102'
    vep_genome = vep_genome ?: 'GRCh38'
    vep_species = vep_species ?: 'homo_sapiens'

    println "Using VEP cache version: $vep_cache_version for assembly $vep_genome and species $vep_species"

    def use_database_mode =  params.extra_vep_args.contains('--database')

    if (use_database_mode) {
        println "Detected '--database' mode: skipping VEP cache download."
        cache_dir = []
    } else if (!(vep_cache && file(vep_cache).exists())) {
        println "VEP cache directory not found. Will download and extract."
        ENSEMBLVEP_DOWNLOAD(
            Channel.of([
                ["id": "${vep_species}_${vep_cache_version}"],
                vep_genome,
                vep_species,
                vep_cache_version
            ])
        )
        cache_dir = ENSEMBLVEP_DOWNLOAD.out.cache.map { meta, cache -> [cache] }
    } else {
        println "VEP cache found at $vep_cache"
        cache_dir = file("$vep_cache")
    }

    // Check VEP plugins directory
    if (!(vep_plugins && file(vep_plugins).exists())) {
        println "VEP plugins directory not found. Will clone from GitHub."
        DOWNLOAD_VEP_PLUGINS(["id": "VEP_plugins"])
        plugins_dir = DOWNLOAD_VEP_PLUGINS.out.vep_plugins.map{ meta, plugins -> [plugins] }
    } else {
        println "VEP plugins found at $vep_plugins"
        plugins_dir = file("$vep_plugins")
    }

    emit:
    vep_cache = cache_dir
    vep_plugins = plugins_dir

    vep_cache_version = vep_cache_version
    vep_genome = vep_genome
    vep_species = vep_species
}
