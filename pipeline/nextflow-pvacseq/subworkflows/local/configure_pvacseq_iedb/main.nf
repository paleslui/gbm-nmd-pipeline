include { SMART_LINK_IEDB        } from '../../../modules/local/smart_link/main'

process DOWNLOAD_MHC_I {
    tag "download_mhc_i"
    label 'process_single'
    publishDir "${params.outdir}/download", mode: 'symlink'

    input:
    path pvacseq_iedb_dir

    output:
    path "$pvacseq_iedb_dir/mhc_i", emit: iedb_mhc_i

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p "$pvacseq_iedb_dir" 2>/dev/null || true
    if [ -d "$pvacseq_iedb_dir/mhc_i" ]; then
        # already present — just emit it
        :
    else
        wget -q https://downloads.iedb.org/tools/mhci/3.1.6/IEDB_MHC_I-3.1.6.tar.gz -O IEDB_MHC_I-3.1.6.tar.gz
        tar -zxf IEDB_MHC_I-3.1.6.tar.gz
        mv mhc_i "$pvacseq_iedb_dir/mhc_i"
        rm -f IEDB_MHC_I-3.1.6.tar.gz
    fi
    """
}

process DOWNLOAD_MHC_II {
    tag "download_mhc_ii"
    label 'process_single'
    publishDir "${params.outdir}/download", mode: 'symlink'

    input:
    path pvacseq_iedb_dir

    output:
    path "$pvacseq_iedb_dir/mhc_ii", emit: iedb_mhc_ii

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p "$pvacseq_iedb_dir" 2>/dev/null || true
    if [ -d "$pvacseq_iedb_dir/mhc_ii" ]; then
        :
    else
        wget -q https://downloads.iedb.org/tools/mhcii/3.1.12/IEDB_MHC_II-3.1.12.tar.gz -O IEDB_MHC_II-3.1.12.tar.gz
        tar -zxf IEDB_MHC_II-3.1.12.tar.gz
        mv mhc_ii "$pvacseq_iedb_dir/mhc_ii"
        rm -f IEDB_MHC_II-3.1.12.tar.gz
    fi
    """
}

//
// Workflow: Configure pVACseq tools
//
workflow CONFIGURE_PVACSEQ_IEDB {
    take:
    pvacseq_iedb_dir         // directory with IEDB for pVACseq
    pvacseq_algorithm        // string: algorithms for pVACseq

    main:

    def valid_mhc_i_algorithms = [
        'BigMHC_EL', 'BigMHC_IM', 'DeepImmuno', 'MHCflurry', 'MHCflurryEL',
        'MHCnuggetsI', 'NNalign', 'NetMHC', 'NetMHCpan', 'NetMHCpanEL',
        'PickPocket', 'SMM', 'SMMPMBEC', 'SMMalign', 'all', 'all_class_i'
    ]

    def valid_mhc_ii_algorithms = [
        'MHCnuggetsII', 'NetMHCIIpan', 'NetMHCIIpanEL', 'NetMHCcons',
        'all', 'all_class_ii'
    ]

    def requires_mhc_i = false
    def requires_mhc_ii = false

    // Determine which MHC classes are required
    pvacseq_algorithm.split(' ').each { algorithm ->
        if (valid_mhc_i_algorithms.contains(algorithm)) {
            requires_mhc_i = true
        }
        if (valid_mhc_ii_algorithms.contains(algorithm)) {
            requires_mhc_ii = true
        }
    }

    // If we dont have any requirements, the algorithm string is wrong
    if (!requires_mhc_i && !requires_mhc_ii) {
        throw new IllegalArgumentException("Invalid algorithm string: '${pvacseq_algorithm}'. It must match at least one valid MHC class I or II algorithm.")
    }

    // Create iedb directory if needed
    if (!(pvacseq_iedb_dir)) {
        println "IEDB is not indicated"
        iedb_dir = file("$params.outdir/iedb")
        iedb_dir.mkdirs()
    } else {
        iedb_dir = file("$pvacseq_iedb_dir")
        if (!iedb_dir.exists()) {
            println "IEDB indicated, but does not exists"
            iedb_dir.mkdirs()
        }
    }
    println "IEDB folder will be $iedb_dir"

    def classI = [
        'BigMHC_EL','BigMHC_IM','DeepImmuno','MHCflurry','MHCflurryEL',
        'MHCnuggetsI','NNalign','NetMHC','NetMHCpan','NetMHCpanEL',
        'PickPocket','SMM','SMMPMBEC','SMMalign','all','all_class_i'
    ]
    def classII = [
        'MHCnuggetsII','NetMHCIIpan','NetMHCIIpanEL','NetMHCcons',
        'all','all_class_ii'
    ]

    def needs_i  = pvacseq_algorithm.tokenize(' ').any { it in classI }
    def needs_ii = pvacseq_algorithm.tokenize(' ').any { it in classII }

    // only dereference `.out` when the process actually runs; otherwise -> empty channel
    def mhc_i  = []
    def mhc_ii = []
    if (needs_i) {
        DOWNLOAD_MHC_I(iedb_dir)
        mhc_i = DOWNLOAD_MHC_I.out.iedb_mhc_i
    }
    if (needs_ii) {
        DOWNLOAD_MHC_II(iedb_dir)
        mhc_ii = DOWNLOAD_MHC_II.out.iedb_mhc_ii
    }

    // smart-link waits for whatever is required
    SMART_LINK_IEDB(iedb_dir, mhc_i, mhc_ii)

    // Parse stdout to two channels
    lines_ch       = SMART_LINK_IEDB.out.iedb_stdout.map { s -> s.readLines().findAll { it?.trim() } }
    iedb_dir_short = lines_ch.map { ls -> file(ls[0].trim()) }                    // Path channel
    link_mode      = lines_ch.map { ls -> ls.size() > 1 ? ls[1].trim() : 'original' }

    // Derive mhc_i/mhc_ii paths as channels
    mhc_i_path  = iedb_dir_short.map { p -> p.resolve('mhc_i') }
    mhc_ii_path = iedb_dir_short.map { p -> p.resolve('mhc_ii') }

    emit:
    iedb_dir = iedb_dir_short
    iedb_mhc_i = mhc_i_path
    iedb_mhc_ii = mhc_ii_path
    mode = link_mode
}
