nextflow.enable.dsl=2

import groovy.json.JsonOutput

// ── Include modules ───────────────────────────────────────────────────────────
include { FASTQ_TRIM_FASTP_FASTQC as fastq_preprocess } from './subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { MULTIQC                 as multiqc          } from './modules/nf-core/multiqc/main'
include { KALLISTO_QUANT          as kallisto         } from './modules/nf-core/kallisto/quant/main'
include { PREPARE_REFERENCE       as prepare_reference } from './subworkflows/local/prepare_reference/main'
include { de_analysis                                 } from './modules/de_analysis.nf'
include { gsea                                        } from './modules/gsea.nf'
include { ora                                         } from './modules/ora.nf'

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

if (!params.metadata_file) {
    error """
    =========================================
    ERROR: --metadata_file is required.

    Set it in your params YAML:
      metadata_file: data/my_project/metadata/sample_metadata.csv
    or pass it on the command line:
      --metadata_file data/my_project/metadata/sample_metadata.csv
    =========================================
    """
}

// ── Resolve paths ────────────────────────────────────────────────────────────
def resolved_metadata = params.metadata_file
def resolved_ref_dir = params.ref_dir
    ?: "${projectDir}/reference"
def resolved_gtf = "${resolved_ref_dir}/Homo_sapiens.GRCh38.113.gtf"
def resolved_genome_fasta = "${resolved_ref_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
def resolved_transcripts_fasta = "${resolved_ref_dir}/transcripts.GRCh38.113.fa"
def resolved_kallisto_index = "${resolved_ref_dir}/transcripts.GRCh38.113.fa.idx"
def resolved_tx2gene = "${resolved_ref_dir}/tx2gene.tsv"
def resolved_gtf_url = params.gtf_url
    ?: 'https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz'
def resolved_genome_fasta_url = params.genome_fasta_url
    ?: 'https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    if (!params.comparisons) {
        error "params.comparisons is required. Define it in your -params-file YAML."
    }

    def comparisons_json = JsonOutput.toJson(params.comparisons)

    prepared_reference = prepare_reference(
        resolved_ref_dir,
        resolved_gtf,
        resolved_genome_fasta,
        resolved_transcripts_fasta,
        resolved_kallisto_index,
        resolved_tx2gene,
        resolved_gtf_url,
        resolved_genome_fasta_url,
        params.reference_auto_download,
        file("${projectDir}/scripts/")
    )

    // Build sample channel from metadata CSV
    samples = Channel
        .fromPath(resolved_metadata, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [
                id        : row.sample_id,
                single_end: false
            ]

            if (row.containsKey('strandedness') && row.strandedness) {
                meta.strandedness = row.strandedness
            }

            tuple(
                meta,
                [
                    file(row.fastq1, checkIfExists: true),
                    file(row.fastq2, checkIfExists: true)
                ]
            )
        }

    // ── Read QC and trimming ──────────────────────────────────────────────────
    qc_trimmed_reads = fastq_preprocess(
        samples.map { meta, reads -> tuple(meta, reads, []) },
        false,
        false,
        false,
        params.skip_trimming,
        false
    )

    reads_for_quant = qc_trimmed_reads.reads

    // ── Kallisto quantification ───────────────────────────────────────────────
    kallisto_out = kallisto(
        reads_for_quant,
        prepared_reference.kallisto_index.map { idx -> tuple([id: 'reference'], idx) },
        [],
        [],
        [],
        []
    )

    // ── MultiQC ───────────────────────────────────────────────────────────────
    all_qc_reports = qc_trimmed_reads.fastqc_raw_zip
        .map { _meta, report -> report }
        .mix(qc_trimmed_reads.fastqc_trim_zip.map { _meta, report -> report })
        .mix(qc_trimmed_reads.trim_json.map { _meta, report -> report })
        .mix(qc_trimmed_reads.trim_html.map { _meta, report -> report })
        .mix(qc_trimmed_reads.trim_log.map { _meta, report -> report })
        .mix(kallisto_out.log.map { _meta, log_file -> log_file })
        .collect()

    multiqc(
        all_qc_reports.map { reports ->
            tuple([id: 'multiqc'], reports, [], [], [], [])
        }
    )

    // ── Differential expression ───────────────────────────────────────────────
    quant_sample_ids = kallisto_out.results.map { meta, _dir -> meta.id }.collect()
    quant_abundance  = kallisto_out.results.map { _meta, dir -> dir }.collect()

    de_results = de_analysis(
        quant_sample_ids,
        quant_abundance,
        file(resolved_metadata),
        Channel.value(comparisons_json),
        Channel.value(params.de_method),
        file("${projectDir}/scripts/"),
        prepared_reference.tx2gene
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
