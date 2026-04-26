process SMART_LINK_IEDB {
    tag "smart_link_iedb"
    label 'process_single'

    input:
    path iedb_source
    path iedb_mhc_i
    path iedb_mhc_ii

    output:
    stdout emit: iedb_stdout

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!(iedb_mhc_i || iedb_mhc_ii)) {
        throw new IllegalArgumentException("Error: No iedb_mhc_i or iedb_mhc_ii provided. At least one is required.")
    }
    """
    smart_link_iedb.py --src "${iedb_source}"
    """
}
