process fastp {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)
    val result_dir

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: trimmed_reads
    path("${sample_id}.html"), emit: html_reports
    path("${sample_id}.json"), emit: json_reports

    publishDir { "${result_dir}/fastp_trim" }, mode: 'copy'

    script:
    """
    fastp \
        -i $read1 \
        -I $read2 \
        -o ${sample_id}_R1.trimmed.fastq.gz \
        -O ${sample_id}_R2.trimmed.fastq.gz \
        --html ${sample_id}.html \
        --json ${sample_id}.json \
        --thread 4
    """
}
