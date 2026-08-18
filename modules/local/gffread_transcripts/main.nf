process gffread_transcripts {

    tag "${output_name}"

    container 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4'

    input:
    path gtf
    path genome_fasta
    val output_name

    output:
    path "${output_name}", emit: transcripts_fasta

    script:
    """
    gffread -w "${output_name}" -g "${genome_fasta}" "${gtf}"
    """
}
