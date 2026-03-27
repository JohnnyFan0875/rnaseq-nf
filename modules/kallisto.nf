process kallisto {

    tag "$sample_id"

    // publishDir is managed centrally in conf/modules.config

    input:
    tuple val(sample_id), path(read1), path(read2)
    path kallisto_index

    output:
    tuple val(sample_id), path("${sample_id}"),                    emit: quant_results
    path("${sample_id}/kallisto_${sample_id}.log"),                emit: quant_log

    script:
    """
    mkdir -p ${sample_id}
    kallisto quant \\
        -i ${kallisto_index} \\
        -o ${sample_id} \\
        -t ${task.cpus} \\
        ${read1} \\
        ${read2} 2>&1 | tee ${sample_id}/kallisto_${sample_id}.log
    """
}
