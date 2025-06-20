nextflow.enable.dsl=2

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

def project_dir = "./data/${params.project_name}"
def config_path = "${project_dir}/config/project_config.yaml"

def parseYamlFile(path) {
    def yamlText = file(path).text
    def slurper = new groovy.yaml.YamlSlurper()
    return slurper.parseText(yamlText)
}

def project_cfg = parseYamlFile(config_path)

import groovy.json.JsonOutput

params.project_dir = project_dir
params.metadata_file = project_cfg.metadata_file
params.kallisto_index = project_cfg.kallisto_index
params.kallisto_threads = project_cfg.kallisto_threads ?: 4
params.de_method = project_cfg.de_method
params.comparisons = project_cfg.comparisons
params.result_dir = "${params.project_dir}/results"
params.comparisons_json = JsonOutput.toJson(project_cfg.comparisons)

// main workflow
workflow {

    // Create a channel for the R script file(s)
    def edgeR_script = file('scripts/edgeR.R')

    samples = Channel
        .fromPath("${params.metadata_file}")
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.sample_id, 
                  file("${params.project_dir}/raw/${row.sample_id}_R1.fastq.gz"), 
                  file("${params.project_dir}/raw/${row.sample_id}_R2.fastq.gz")
            )
        }

    // Pre-trim QC
    fastqc_pre_out = fastqc_pre(samples, Channel.value("pre"), Channel.value(params.result_dir))

    // Fastp trimming
    trimmed = fastp(samples, Channel.value(params.result_dir))

    // Post-trim QC
    fastqc_post_out = fastqc_post(trimmed.trimmed_reads, Channel.value("post"), Channel.value(params.result_dir))

    // Kallisto quantification
    kallisto_out = kallisto(trimmed.trimmed_reads, 
        file(params.kallisto_index), 
        Channel.value(params.result_dir),
        Channel.value(params.kallisto_threads)
    )

    // multiqc (without post-trimming fastqc)
    all_qc_reports = fastqc_pre_out.fastqc_reports
        .mix(trimmed.html_reports)
        .mix(trimmed.json_reports)
        .mix(kallisto_out.quant_log)
        .collect()

    multiqc(all_qc_reports, Channel.value(params.result_dir))

    // Differential expression analysis
    de_results = de_analysis(
        file("${params.result_dir}/kallisto_quant"),
        file(params.metadata_file),
        Channel.value(params.comparisons_json),
        Channel.value(params.de_method),
        Channel.value(params.result_dir),
        file("scripts/"),
        file("reference/")
    )

    // Define enrichment analysis inputs
    enrichment_inputs = de_results.deg_results
        .map { deg_dir -> tuple(
            deg_dir,
            file(params.metadata_file),
            params.comparisons_json,
            params.result_dir,
            params.de_method,
            file("scripts/")
        )}

    // GSEA analysis
    gsea_results = gsea(enrichment_inputs)

    // ORA analysis
    ora_results = ora(enrichment_inputs)

}
