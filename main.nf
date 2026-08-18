nextflow.enable.dsl=2

import groovy.json.JsonOutput

// ── Include modules ───────────────────────────────────────────────────────────
include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc/main'
include { FASTP } from './modules/nf-core/fastp/main'
include { MULTIQC as multiqc } from './modules/nf-core/multiqc/main'
include { KALLISTO_QUANT as kallisto } from './modules/nf-core/kallisto/quant/main'
include { download_reference_asset as download_gtf_asset } from './modules/local/download_reference_asset/main'
include { download_reference_asset as download_genome_fasta_asset } from './modules/local/download_reference_asset/main'
include { gffread_transcripts } from './modules/local/gffread_transcripts/main'
include { kallisto_index_build } from './modules/local/kallisto_index_build/main'
include { salmon_index_build } from './modules/local/salmon_index_build/main'
include { salmon_quant } from './modules/local/salmon_quant/main'
include { deg_analysis } from './modules/local/deg_analysis/main'
include { gsea } from './modules/local/gsea/main'
include { ora } from './modules/local/ora/main'

def getFastpReadsAfterFiltering(json_file) {
    if (workflow.stubRun) {
        return 1
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toLong()
}

def resolveInputPath(pathString) {
    if (!pathString) {
        return null
    }

    def rawPath = pathString.toString()
    def direct = file(rawPath)
    if (direct.exists()) {
        return direct
    }

    if (rawPath.startsWith('/mnt/d/rna_seq_pipeline/')) {
        def remapped = rawPath.replace('/mnt/d/rna_seq_pipeline', projectDir.toString())
        def remappedFile = file(remapped)
        if (remappedFile.exists()) {
            return remappedFile
        }
    }

    return file(rawPath, checkIfExists: true)
}

def existingFile(String candidate) {
    def persisted = file(candidate)
    return persisted.exists() && persisted.size() > 0
}

def existingDirectory(String candidate, String sentinel) {
    def persisted = new File(candidate, sentinel)
    return persisted.exists() && persisted.length() > 0
}

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

def supported_quant_strategies = ['kallisto', 'salmon']
if (!supported_quant_strategies.contains(params.quant_strategy)) {
    error "Unsupported --quant_strategy '${params.quant_strategy}'. Supported values: ${supported_quant_strategies.join(', ')}"
}

// ── Resolve paths ────────────────────────────────────────────────────────────
def resolved_metadata = params.metadata_file
def resolved_ref_dir = params.ref_dir ?: "${projectDir}/reference"
def resolved_gtf = "${resolved_ref_dir}/Homo_sapiens.GRCh38.113.gtf"
def resolved_genome_fasta = "${resolved_ref_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
def resolved_transcripts_fasta = "${resolved_ref_dir}/transcripts.GRCh38.113.fa"
def resolved_kallisto_index = "${resolved_ref_dir}/transcripts.GRCh38.113.kallisto.idx"
def resolved_salmon_index = "${resolved_ref_dir}/salmon_index"
def resolved_gtf_url = params.gtf_url ?: 'https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz'
def resolved_genome_fasta_url = params.genome_fasta_url ?: 'https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

def needs_transcriptome = true
def de_input_mode = 'tximport'
def quant_type = params.quant_strategy

def has_gtf = existingFile(resolved_gtf)
def has_genome_fasta = existingFile(resolved_genome_fasta)
def has_transcripts_fasta = existingFile(resolved_transcripts_fasta)
def has_kallisto_index = existingFile(resolved_kallisto_index)
def has_salmon_index = existingDirectory(resolved_salmon_index, 'versionInfo.json')

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    if (!params.comparisons) {
        error "params.comparisons is required. Define it in your -params-file YAML."
    }

    def comparisons_json = JsonOutput.toJson(params.comparisons)

    gtf_channel = has_gtf
        ? Channel.value(file(resolved_gtf, checkIfExists: true))
        : download_gtf_asset('Homo_sapiens.GRCh38.113.gtf', resolved_gtf_url, params.reference_auto_download).gtf

    genome_fasta_channel = has_genome_fasta
        ? Channel.value(file(resolved_genome_fasta, checkIfExists: true))
        : download_genome_fasta_asset('Homo_sapiens.GRCh38.dna.primary_assembly.fa', resolved_genome_fasta_url, params.reference_auto_download).fasta

    transcripts_fasta_channel = needs_transcriptome
        ? (
            has_transcripts_fasta
                ? Channel.value(file(resolved_transcripts_fasta, checkIfExists: true))
                : gffread_transcripts(gtf_channel, genome_fasta_channel, 'transcripts.GRCh38.113.fa').transcripts_fasta
        )
        : Channel.empty()

    kallisto_index_channel = params.quant_strategy == 'kallisto'
        ? (
            has_kallisto_index
                ? Channel.value(file(resolved_kallisto_index, checkIfExists: true))
                : kallisto_index_build(transcripts_fasta_channel, 'transcripts.GRCh38.113.kallisto.idx').kallisto_index
        )
        : Channel.empty()

    salmon_index_channel = params.quant_strategy == 'salmon'
        ? (
            has_salmon_index
                ? Channel.value(file(resolved_salmon_index, checkIfExists: true))
                : salmon_index_build(transcripts_fasta_channel).salmon_index
        )
        : Channel.empty()

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
                    resolveInputPath(row.fastq1),
                    resolveInputPath(row.fastq2)
                ]
            )
        }

    // ── Read QC and trimming ──────────────────────────────────────────────────
    reads_with_adapters = samples.map { meta, reads -> tuple(meta, reads, []) }
    reads_for_fastqc_raw = reads_with_adapters.map { meta, reads, _adapter_fasta -> tuple(meta, reads) }

    FASTQC_RAW(reads_for_fastqc_raw)

    trim_json = Channel.empty()
    trim_html = Channel.empty()
    trim_log = Channel.empty()
    fastqc_trim_zip = Channel.empty()

    if (params.skip_trimming) {
        reads_for_quant = reads_for_fastqc_raw
    } else {
        FASTP(
            reads_with_adapters,
            false,
            false,
            false
        )

        trim_json = FASTP.out.json
        trim_html = FASTP.out.html
        trim_log = FASTP.out.log

        reads_for_quant = FASTP.out.reads
            .join(FASTP.out.json)
            .map { meta, reads, json ->
                if (getFastpReadsAfterFiltering(json) > 0) {
                    tuple(meta, reads)
                }
            }

        FASTQC_TRIM(reads_for_quant)
        fastqc_trim_zip = FASTQC_TRIM.out.zip
    }

    quant_result_dirs = Channel.empty()
    quant_aux_reports = Channel.empty()
    if (params.quant_strategy == 'kallisto') {
        kallisto_out = kallisto(
            reads_for_quant,
            kallisto_index_channel.map { idx -> tuple([id: 'reference'], idx) },
            [],
            [],
            [],
            []
        )

        quant_result_dirs = kallisto_out.results.map { _meta, dir -> dir }
        quant_aux_reports = kallisto_out.log.map { _meta, log_file -> log_file }
    } else if (params.quant_strategy == 'salmon') {
        salmon_out = salmon_quant(
            reads_for_quant,
            salmon_index_channel,
            transcripts_fasta_channel
        )

        quant_result_dirs = salmon_out.results.map { _meta, dir -> dir }
        quant_aux_reports = salmon_out.meta_info
            .map { _meta, report -> report }
            .mix(salmon_out.lib_format_counts.map { _meta, report -> report })
            .mix(salmon_out.log.map { _meta, report -> report })
    }

    // ── MultiQC ───────────────────────────────────────────────────────────────
    all_qc_reports = FASTQC_RAW.out.zip
        .map { _meta, report -> report }
        .mix(fastqc_trim_zip.map { _meta, report -> report })
        .mix(trim_json.map { _meta, report -> report })
        .mix(trim_html.map { _meta, report -> report })
        .mix(trim_log.map { _meta, report -> report })
        .mix(quant_aux_reports)
        .collect()

    multiqc(
        all_qc_reports.map { reports ->
            tuple([id: 'multiqc'], reports, [], [], [], [])
        }
    )

    // ── Differential expression ───────────────────────────────────────────────
    expression_inputs = quant_result_dirs.collect()

    de_results = deg_analysis(
        expression_inputs,
        file(resolved_metadata),
        Channel.value(comparisons_json),
        Channel.value(params.de_method),
        file("${projectDir}/scripts/"),
        gtf_channel,
        Channel.value(de_input_mode),
        Channel.value(quant_type)
    )

    // ── Enrichment analysis ───────────────────────────────────────────────────
    if (!params.skip_enrichment) {
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
