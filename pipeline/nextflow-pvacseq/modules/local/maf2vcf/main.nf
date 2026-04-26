// Based on nf-core module vcf2maf https://nf-co.re/modules/vcf2maf
process MAF2VCF {

    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b6fc09bed47d0dc4d8384ce9e04af5806f2cc91b:305092c6f8420acd17377d2cc8b96e1c3ccb7d26-0':
        'biocontainers/mulled-v2-b6fc09bed47d0dc4d8384ce9e04af5806f2cc91b:305092c6f8420acd17377d2cc8b96e1c3ccb7d26-0' }"

    input:
    tuple val(meta), path(maf) 
    path fasta                 // Required

    output:
    tuple val(meta), path("*.vcf"), path("*.pairs.tsv")   , emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    def VERSION = '1.6.21' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    maf2vcf.pl \\
        $args \\
        --ref-fasta $fasta \\
        --input-maf $maf \\
        --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maf2vcf: $VERSION
    END_VERSIONS
    """
}