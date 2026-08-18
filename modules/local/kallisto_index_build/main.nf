process kallisto_index_build {

    tag "${output_name}"

    container 'quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0'

    input:
    path transcripts_fasta
    val output_name

    output:
    path "${output_name}", emit: kallisto_index

    script:
    """
    kallisto index -i "${output_name}" "${transcripts_fasta}"
    """
}
