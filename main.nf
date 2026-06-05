nextflow.enable.dsl=2

import groovy.json.JsonOutput

// ── Include modules ───────────────────────────────────────────────────────────
include { FASTQC as FASTQC_RAW                           } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM                          } from './modules/nf-core/fastqc/main'
include { FASTP                                          } from './modules/nf-core/fastp/main'
include { MULTIQC                 as multiqc           } from './modules/nf-core/multiqc/main'
include { KALLISTO_QUANT          as kallisto          } from './modules/nf-core/kallisto/quant/main'
include { prepare_reference_assets                        } from './modules/local/prepare_reference_assets/main'
include { deg_analysis                                 } from './modules/local/deg_analysis/main'
include { gsea                                         } from './modules/local/gsea/main'
include { ora                                          } from './modules/local/ora/main'

def getFastpReadsAfterFiltering(json_file) {
    if (workflow.stubRun) {
        return 1
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toLong()
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

// ── Resolve paths ────────────────────────────────────────────────────────────
def resolved_metadata = params.metadata_file
def resolved_ref_dir = params.ref_dir
    ?: "${projectDir}/reference"
def resolved_gtf = "${resolved_ref_dir}/Homo_sapiens.GRCh38.113.gtf"
def resolved_genome_fasta = "${resolved_ref_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
def resolved_transcripts_fasta = "${resolved_ref_dir}/transcripts.GRCh38.113.fa"
def resolved_kallisto_index = "${resolved_ref_dir}/transcripts.GRCh38.113.kallisto.idx"
def resolved_gtf_url = params.gtf_url
    ?: 'https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz'
def resolved_genome_fasta_url = params.genome_fasta_url
    ?: 'https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
def persisted_reference_candidates = [
    resolved_gtf,
    resolved_genome_fasta,
    resolved_transcripts_fasta,
    resolved_kallisto_index
]
def has_persisted_reference_assets = [
    resolved_gtf,
    resolved_genome_fasta,
    resolved_transcripts_fasta,
    resolved_kallisto_index
].every { candidate ->
    def persisted_file = file(candidate)
    persisted_file.exists() && persisted_file.size() > 0
}
def existing_reference_files = persisted_reference_candidates
    .collect { file(it) }
    .findAll { persisted_file -> persisted_file.exists() && persisted_file.size() > 0 }

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    if (!params.comparisons) {
        error "params.comparisons is required. Define it in your -params-file YAML."
    }

    def comparisons_json = JsonOutput.toJson(params.comparisons)

    prepared_reference = has_persisted_reference_assets
        ? [
            gtf              : Channel.value(file(resolved_gtf, checkIfExists: true)),
            genome_fasta     : Channel.value(file(resolved_genome_fasta, checkIfExists: true)),
            transcripts_fasta: Channel.value(file(resolved_transcripts_fasta, checkIfExists: true)),
            kallisto_index   : Channel.value(file(resolved_kallisto_index, checkIfExists: true))
        ]
        : prepare_reference_assets(
            resolved_ref_dir,
            existing_reference_files,
            resolved_gtf,
            resolved_genome_fasta,
            resolved_transcripts_fasta,
            resolved_kallisto_index,
            resolved_gtf_url,
            resolved_genome_fasta_url,
            params.reference_auto_download
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
    all_qc_reports = FASTQC_RAW.out.zip
        .map { _meta, report -> report }
        .mix(fastqc_trim_zip.map { _meta, report -> report })
        .mix(trim_json.map { _meta, report -> report })
        .mix(trim_html.map { _meta, report -> report })
        .mix(trim_log.map { _meta, report -> report })
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

    de_results = deg_analysis(
        quant_sample_ids,
        quant_abundance,
        file(resolved_metadata),
        Channel.value(comparisons_json),
        Channel.value(params.de_method),
        file("${projectDir}/scripts/"),
        prepared_reference.gtf
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
