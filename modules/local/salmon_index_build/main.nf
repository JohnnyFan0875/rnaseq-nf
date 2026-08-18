process salmon_index_build {

    tag "salmon_index"

    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'

    input:
    path transcripts_fasta

    output:
    path 'salmon_index', emit: salmon_index

    script:
    """
    salmon index \
        --threads ${task.cpus} \
        -t "${transcripts_fasta}" \
        -i salmon_index
    """
}
