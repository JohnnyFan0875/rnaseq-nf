process ora {

    tag "ORA"

    input:
    tuple path(deg_dir), val(comparisons_json), val(de_method), path(script_dir)

    output:
    path "ora_results", emit: ora_out

    script:
    """
    mkdir -p ora_results
    Rscript ${script_dir}/ora.R \\
        --deg_dir ${deg_dir} \\
        --out_dir ora_results/ \\
        --script_dir ${script_dir} \\
        --comparisons '${comparisons_json}' \\
        --de_method '${de_method}' \\
        --organism '${params.organism}'
    """

    stub:
    """
    mkdir -p ora_results
    touch ora_results/.stub
    """
}
