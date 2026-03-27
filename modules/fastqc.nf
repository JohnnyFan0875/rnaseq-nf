process fastqc {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)
    val mode        // "pre" or "post" — kept for tag/log purposes only

    output:
    path "*_fastqc.*", emit: fastqc_reports

    script:
    """
    fastqc $read1 $read2 --threads ${task.cpus}
    """
}
