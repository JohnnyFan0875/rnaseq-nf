process fastqc {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)
    val mode
    val result_dir

    output:
    path "*_fastqc.*", emit: fastqc_reports

    publishDir { "${result_dir}/fastqc_${mode}" }, mode: 'copy'

    script:
    """
    fastqc $read1 $read2 --threads 4
    """
}
