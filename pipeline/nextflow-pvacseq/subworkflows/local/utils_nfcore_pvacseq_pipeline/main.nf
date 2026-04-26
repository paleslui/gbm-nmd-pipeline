//
// Subworkflow with functionality specific to the nextflow_pvacseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin/main'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    emit:
    versions       = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/


workflow PIPELINE_COMPLETION {

    take:
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    iedb_dir        // String absolute path
    mode            // "original" | "hardlink" | "copy"

    main:
    workflow.onComplete {
        completionSummary(monochrome_logs)
        def status = unlinkIedb(iedb_dir as String, mode as String)
        log.info "Unlink iedb ${status} (dir=${iedb_dir}, mode=${mode})"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "pVACtools (Hundal, Jasreet, et al. 2020)",
            "VEP (McLaren W et al. 2016)",
            "vcf2maf (Cyriac Kandoth, 2020)",
            "MultiQC (Ewels et al. 2016)",
            "nf-core template (Philip Ewels et al., 2020)"
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
            "<li>Jasreet Hundal+, Susanna Kiwala+, Joshua McMichael, Christopher A Miller, Alexander T Wollam, Huiming Xia, Connor J Liu, Sidi Zhao, Yang-Yang Feng, Aaron P Graubert, Amber Z Wollam, Jonas Neichin, Megan Neveau, Jason Walker, William E Gillanders, Elaine R Mardis, Obi L Griffith, Malachi Griffith. pVACtools: a computational toolkit to select and visualize cancer neoantigens. Cancer Immunology Research. 2020 Mar;8(3):409-420. DOI: 10.1158/2326-6066.CIR-19-0401. PMID: 31907209. (+) equal contribution. </li>",
            "<li>Jasreet Hundal, Susanna Kiwala, Yang-Yang Feng, Connor J. Liu, Ramaswamy Govindan, William C. Chapman, Ravindra Uppaluri, S. Joshua Swamidass, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. Accounting for proximal variants improves neoantigen prediction. Nature Genetics. 2018, DOI: 10.1038/s41588-018-0283-9. PMID: 30510237.</li>",
            "<li>Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. pVACseq: A genome-guided in silico approach to identifying tumor neoantigens. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: 26825632.</li>",
            "<li>McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P Cunningham F. The Ensembl Variant Effect Predictor. Genome Biology Jun 6;17(1):122. (2016) doi:10.1186/s13059-016-0974-4</li>",
            "<li>Cyriac Kandoth. mskcc/vcf2maf: vcf2maf v1.6. (2020). doi:10.5281/zenodo.593251</li>",
            "<li>Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen. The nf-core framework for community-curated bioinformatics pipelines.Nature Biotechnology. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = """
<ul>

  <li><a href="https://www.ensembl.org/info/docs/tools/vep/index.html">VEP</a>:
      McLaren W, Gil L, Hunt SE, et&nbsp;al. <em>Genome Biology</em> (2016) 17:122.
      doi: <a href="https://doi.org/10.1186/s13059-016-0974-4">10.1186/s13059-016-0974-4</a></li>
  <li><a href="https://github.com/mskcc/vcf2maf">vcf2maf</a>:
      Kandoth C. <em>Zenodo</em> (2020).
      doi: <a href="https://doi.org/10.5281/zenodo.593251">10.5281/zenodo.593251</a></li>
  <li><a href="https://nf-co.re/">nf-core template</a>:
      Ewels P, Peltzer A, Fillinger S, et&nbsp;al. <em>Nature Biotechnology</em> (2020).
      doi: <a href="https://doi.org/10.1038/s41587-020-0439-x">10.1038/s41587-020-0439-x</a></li>
</ul>
""".stripIndent().trim()
    meta["tool_bibliography"] = """
<p><strong>pVACseq / pVACtools</strong><br/>
Hundal J, Kiwala S, McMichael J, Miller CA, Wollam AT, Xia H, Liu CJ, Zhao S, Feng Y-Y, Graubert AP, Wollam AZ, Neichin J, Neveau M, Walker J, Gillanders WE, Mardis ER, Griffith OL, Griffith M.
pVACtools: a computational toolkit to select and visualize cancer neoantigens. <em>Cancer Immunology Research</em>. 2020;8(3):409–420.
doi: <a href="https://doi.org/10.1158/2326-6066.CIR-19-0401">10.1158/2326-6066.CIR-19-0401</a>.<br/>
Hundal J, Kiwala S, Feng Y-Y, Liu CJ, Govindan R, Chapman WC, Uppaluri R, Swamidass SJ, Griffith OL, Mardis ER, Griffith M.
Accounting for proximal variants improves neoantigen prediction. <em>Nature Genetics</em>. 2018.
doi: <a href="https://doi.org/10.1038/s41588-018-0283-9">10.1038/s41588-018-0283-9</a>.<br/>
Hundal J, Carreno BM, Petti AA, Linette GP, Griffith OL, Mardis ER, Griffith M.
pVACseq: A genome-guided <em>in silico</em> approach to identifying tumor neoantigens. <em>Genome Medicine</em>. 2016;8:11.
doi: <a href="https://doi.org/10.1186/s13073-016-0264-5">10.1186/s13073-016-0264-5</a>.
</p>

<p><strong>VEP</strong><br/>
McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F.
The Ensembl Variant Effect Predictor. <em>Genome Biology</em>. 2016;17:122.
doi: <a href="https://doi.org/10.1186/s13059-016-0974-4">10.1186/s13059-016-0974-4</a>.
</p>

<p><strong>vcf2maf</strong><br/>
Kandoth C. mskcc/vcf2maf: vcf2maf v1.6. <em>Zenodo</em>. 2020.
doi: <a href="https://doi.org/10.5281/zenodo.593251">10.5281/zenodo.593251</a>.
</p>

<p><strong>nf-core template (used to scaffold this pipeline)</strong><br/>
Ewels P, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S.
The nf-core framework for community-curated bioinformatics pipelines. <em>Nature Biotechnology</em>. 2020 Feb 13.
doi: <a href="https://doi.org/10.1038/s41587-020-0439-x">10.1038/s41587-020-0439-x</a>.
</p>
""".stripIndent().trim()

    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

def unlinkIedb(String chosenDir, String mode) {
    try {
        // basic arg checks
        if (!chosenDir || chosenDir.trim().isEmpty())
            return 'no_path'
        if (!mode || mode == 'original')
            return 'no_cleanup'   // only remove temp dirs created by hardlink/copy

        // resolve path (don’t follow symlinks for existence check)
        def p = java.nio.file.Paths.get(chosenDir.trim())
        if (!java.nio.file.Files.exists(p, java.nio.file.LinkOption.NOFOLLOW_LINKS))
            return 'already_missing'

        try {
            // canonical without following symlinks
            p = p.toRealPath(java.nio.file.LinkOption.NOFOLLOW_LINKS)
        }
        catch (java.nio.file.NoSuchFileException ignored) {
            return 'already_missing'
        }

        // safety rails
        if (!java.nio.file.Files.isDirectory(p, java.nio.file.LinkOption.NOFOLLOW_LINKS))
            return "error: refusing to remove non-directory: ${p}"
        if (p.getParent() == null) // root like "/"
            return "error: refusing to remove root: ${p}"

        // optional extra guard against nuking generic mount points
        def name = p.getFileName()?.toString() ?: ""
        if (name in ['tmp','temp','scratch','mnt','local','fsx','lustre','gpfs'])
            return "error: refusing to remove generic directory name: ${p}"

        // walk depth-first and delete
        def stream = java.nio.file.Files.walk(p)
        def paths  = stream.sorted(java.util.Comparator.reverseOrder())
                           .collect(java.util.stream.Collectors.toList())
        stream.close()

        for (q in paths) {
            java.nio.file.Files.deleteIfExists(q)
        }

        return 'removed'
    }
    catch (Exception e) {
        return "error: ${e.class.simpleName}: ${e.message}"
    }
}
