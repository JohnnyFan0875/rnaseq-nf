nextflow.enable.dsl=2

import groovy.json.JsonOutput

// ── Include modules ───────────────────────────────────────────────────────────
include { fastqc as fastqc_pre  } from './modules/fastqc.nf'
include { fastqc as fastqc_post } from './modules/fastqc.nf'
include { fastp        }          from './modules/fastp.nf'
include { multiqc      }          from './modules/multiqc.nf'
include { kallisto     }          from './modules/kallisto.nf'
include { de_analysis  }          from './modules/de_analysis.nf'
include { gsea         }          from './modules/gsea.nf'
include { ora          }          from './modules/ora.nf'

// ── Parameter validation ──────────────────────────────────────────────────────
if (!params.project_name) {
    error """
    =========================================
    ERROR: --project_name is required.

    Example:
      nextflow run main.nf \\
        -profile docker \\
        -params-file params/my_project.yaml
    =========================================
    """
}

if (!params.outdir) {
    error """
    =========================================
    ERROR: --outdir is required.

    Set it in your params YAML:
      outdir: results/my_project
    or pass it on the command line:
      --outdir results/my_project
    =========================================
    """
}

// ── Resolve metadata_file default ────────────────────────────────────────────
def resolved_metadata = params.metadata_file
    ?: "./data/${params.project_name}/metadata/sample_metadata.csv"

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    if (!params.comparisons) {
        error "params.comparisons is required. Define it in your -params-file YAML."
    }

    def comparisons_json = JsonOutput.toJson(params.comparisons)

    // Build sample channel from metadata CSV
    samples = Channel
        .fromPath(resolved_metadata, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            tuple(
                row.sample_id,
                file("./data/${params.project_name}/raw/${row.fastq1}", checkIfExists: true),
                file("./data/${params.project_name}/raw/${row.fastq2}", checkIfExists: true)
            )
        }

    // ── Pre-trim QC ───────────────────────────────────────────────────────────
    fastqc_pre_out = fastqc_pre(
        samples,
        Channel.value("pre")
    )

    // ── Adapter trimming ──────────────────────────────────────────────────────
    if (!params.skip_trimming) {
        trimmed = fastp(samples)

        fastqc_post_out = fastqc_post(
            trimmed.trimmed_reads,
            Channel.value("post")
        )

        reads_for_quant    = trimmed.trimmed_reads
        fastp_qc_logs      = trimmed.html_reports.mix(trimmed.json_reports)
        post_fastqc_reports = fastqc_post_out.fastqc_reports

    } else {
        fastqc_post_out    = Channel.empty()
        reads_for_quant    = samples
        fastp_qc_logs      = Channel.empty()
        post_fastqc_reports = Channel.empty()
    }

    // ── Kallisto quantification ───────────────────────────────────────────────
    kallisto_out = kallisto(
        reads_for_quant,
        file(params.kallisto_index, checkIfExists: true)
    )

    // ── MultiQC ───────────────────────────────────────────────────────────────
    // FIX: post_fastqc_reports now mixed in — was previously omitted, causing
    // post-trim FastQC data to be absent from the MultiQC report.
    all_qc_reports = fastqc_pre_out.fastqc_reports
        .mix(post_fastqc_reports)
        .mix(fastp_qc_logs)
        .mix(kallisto_out.quant_log)
        .collect()

    multiqc(all_qc_reports)

    // ── Differential expression ───────────────────────────────────────────────
    quant_sample_ids = kallisto_out.quant_results.map { id, _dir -> id  }.collect()
    quant_abundance  = kallisto_out.quant_results.map { _id, dir -> dir }.collect()

    de_results = de_analysis(
        quant_sample_ids,
        quant_abundance,
        file(resolved_metadata),
        Channel.value(comparisons_json),
        Channel.value(params.de_method),
        file("${projectDir}/scripts/"),
        file(params.ref_dir)
    )

    // ── Enrichment analysis ───────────────────────────────────────────────────
    if (!params.skip_enrichment) {

        // FIX: Removed file(resolved_metadata) from enrichment_inputs tuple —
        // gsea.nf and ora.nf no longer declare metadata_file in their input,
        // since it was staged but never consumed by the R scripts.
        enrichment_inputs = de_results.deg_results
            .map { deg_dir -> tuple(
                deg_dir,
                comparisons_json,
                params.de_method,
                file("${projectDir}/scripts/")
            )}

        gsea_results = gsea(enrichment_inputs)
        ora_results  = ora(enrichment_inputs)
    }
}
