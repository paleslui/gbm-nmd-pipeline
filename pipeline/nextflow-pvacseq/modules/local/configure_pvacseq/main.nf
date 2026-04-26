//
// MODULE: Configure pVACseq tools
//
process CONFIGURE_PVACSEQ {
    tag "configure_pvacseq"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/griffithlab/pvactools:5.3.1"

    input:
    path iedb_dir
    path iedb_mhc_i
    path iedb_mhc_ii

    output:
    path iedb_dir, emit: iedb_dir
    path 'env_config_done.txt', emit: config_file

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!(iedb_mhc_i || iedb_mhc_ii)) {
        throw new IllegalArgumentException("Error: No iedb_mhc_i or iedb_mhc_ii provided. At least one is required.")
    }

    // Check and log for MHC Class I
    if (iedb_mhc_i) {
        if (!iedb_mhc_i.isDirectory()) {
            println "Warning: Provided iedb_mhc_i as '${iedb_mhc_i}' is not a directory."
        } else {
            println "MHC Class I will be configured with IEDB in '${iedb_mhc_i}'"
        }
    } else {
        println "MHC Class I configuration skipped."
    }


    // Check and log for MHC Class II
    if (iedb_mhc_ii) {
        if (!iedb_mhc_ii.isDirectory()) {
            println "Warning: Provided iedb_mhc_ii as '${iedb_mhc_ii}' is not a directory."
        } else {
            println "MHC Class II will be configured with IEDB in '${iedb_mhc_ii}'"
        }
    } else {
        println "MHC Class II configuration skipped."
    }

    """
    mhcflurry-downloads fetch

    if [ -d "${iedb_mhc_i}" ]; then
        cd "${iedb_mhc_i}"
        ./configure -k "1"
        cd -
    fi

    if [ -d "${iedb_mhc_ii}" ]; then
        cd "${iedb_mhc_ii}"
        # Run configure.py but ignore failure
        if ! ./configure.py; then
            echo "Warning: configure.py failed"
        fi
        cd -
    fi

    # Generate the stats file
    touch env_config_done.txt
    echo "Configuration details:" >> env_config_done.txt
    if [ -d "${iedb_mhc_i}" ]; then
        echo "- MHC Class I configured at: ${iedb_mhc_i}" >> env_config_done.txt
    else
        echo "- MHC Class I: Not configured." >> env_config_done.txt
    fi
    if [ -d "${iedb_mhc_ii}" ]; then
        echo "- MHC Class II configured at: ${iedb_mhc_ii}" >> env_config_done.txt
    else
        echo "- MHC Class II: Not configured." >> env_config_done.txt
    fi
    """
}
