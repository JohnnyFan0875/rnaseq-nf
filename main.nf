nextflow.enable.dsl=2

import groovy.json.JsonOutput

// Include modules
include { fastqc as fastqc_pre } from './modules/fastqc.nf'
include { fastqc as fastqc_post } from './modules/fastqc.nf'
include { fastp }  from './modules/fastp.nf'
include { multiqc } from './modules/multiqc.nf'
include { kallisto } from './modules/kallisto.nf'
include { de_analysis } from './modules/de_analysis.nf'
include { gsea } from './modules/gsea.nf'
include { ora } from './modules/ora.nf'

// Define parameter
if (!params.project_name) {
    error "Please provide a project name using --project_name"
}

def parseYamlFile(path) {
    def yamlText = file(path).text
    def slurper = new groovy.yaml.YamlSlurper()
    return slurper.parseText(yamlText)
}

def project_dir = "./data/${params.project_name}"
def config_path = "${project_dir}/config/project_config.yaml"
def project_cfg = parseYamlFile(config_path)
def metadata_file = project_cfg.metadata_file
def kallisto_index = project_cfg.kallisto_index
def kallisto_threads = project_cfg.kallisto_threads ?: 4
def de_method = project_cfg.de_method
def comparisons = project_cfg.comparisons
def result_dir = "${project_dir}/results"
def comparisons_json = JsonOutput.toJson(comparisons)

// main workflow
workflow {

    samples = Channel
        .fromPath("${metadata_file}")
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            tuple(row.sample_id, 
                  file("${project_dir}/raw/${row.fastq1}"), 
                  file("${project_dir}/raw/${row.fastq2}")
            )
        }

    // Pre-trim QC
    fastqc_pre_out = fastqc_pre(samples,
        Channel.value("pre"),
        Channel.value(result_dir)
    )

    // Fastp trimming
    trimmed = fastp(samples,
        Channel.value(result_dir)
    )

    // Post-trim QC
    fastqc_post_out = fastqc_post(trimmed.trimmed_reads,
        Channel.value("post"),
        Channel.value(result_dir)
    )

    // Kallisto quantification
    kallisto_out = kallisto(trimmed.trimmed_reads, 
        file(kallisto_index), 
        Channel.value(result_dir),
        Channel.value(kallisto_threads)
    )

    quant_results_all = kallisto_out.quant_results.collect()

    // multiqc (without post-trimming fastqc)
    all_qc_reports = fastqc_pre_out.fastqc_reports
        .mix(trimmed.html_reports)
        .mix(trimmed.json_reports)
        .mix(kallisto_out.quant_log)
        .collect()

    multiqc(all_qc_reports, Channel.value(result_dir))

    // Differential expression analysis
    de_results = de_analysis(
        quant_results_all,
        file(metadata_file),
        Channel.value(comparisons_json),
        Channel.value(de_method),
        Channel.value(result_dir),
        file("scripts/"),
        file("reference/")
    )

    // Define enrichment analysis inputs
    enrichment_inputs = de_results.deg_results
        .map { deg_dir -> tuple(
            deg_dir,
            file(metadata_file),
            comparisons_json,
            result_dir,
            de_method,
            file("scripts/")
        )}
 
    // GSEA analysis
    gsea_results = gsea(enrichment_inputs)
 
    // ORA analysis
    ora_results = ora(enrichment_inputs)

}
