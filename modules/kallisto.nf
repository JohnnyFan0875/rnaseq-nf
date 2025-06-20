process kallisto {
    
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)
    path kallisto_index
    val result_dir
    val kallisto_threads

    output:
    tuple val(sample_id), path("${sample_id}/abundance.tsv"), emit: quant_results
    path("${sample_id}/kallisto_${sample_id}.log"), emit: quant_log

    publishDir { "${result_dir}/kallisto_quant/" }, mode: 'copy'

    script:
    """
    mkdir -p $sample_id
    kallisto quant \
        -i $kallisto_index \
        -o $sample_id \
        -t $kallisto_threads \
        $read1 \
        $read2 2>&1 | tee ${sample_id}/kallisto_${sample_id}.log

    """
}
