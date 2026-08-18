process salmon_quant {

    tag "${meta.id}"

    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'

    input:
    tuple val(meta), path(reads)
    path index
    path transcripts_fasta

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val(meta), path("${prefix}_meta_info.json"), emit: meta_info, optional: true
    tuple val(meta), path("${prefix}_lib_format_counts.json"), emit: lib_format_counts, optional: true
    tuple val(meta), path("${prefix}.salmon.log"), emit: log

    script:
    prefix = "${meta.id}"
    def lib_type = meta.strandedness == 'forward' ? 'ISF' : meta.strandedness == 'reverse' ? 'ISR' : 'IU'

    """
    salmon quant \
        --threads ${task.cpus} \
        --index "${index}" \
        --libType ${lib_type} \
        -1 "${reads[0]}" \
        -2 "${reads[1]}" \
        -o "${prefix}" \
        2> >(tee "${prefix}.salmon.log" >&2)

    cp "${prefix}/aux_info/meta_info.json" "${prefix}_meta_info.json"
    cp "${prefix}/lib_format_counts.json" "${prefix}_lib_format_counts.json"
    """

    stub:
    prefix = "${meta.id}"
    """
    mkdir -p "${prefix}/aux_info"
    touch "${prefix}/quant.sf"
    touch "${prefix}/cmd_info.json"
    touch "${prefix}/lib_format_counts.json"
    touch "${prefix}/aux_info/meta_info.json"
    touch "${prefix}.salmon.log"
    cp "${prefix}/aux_info/meta_info.json" "${prefix}_meta_info.json"
    cp "${prefix}/lib_format_counts.json" "${prefix}_lib_format_counts.json"
    """
}
