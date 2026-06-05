process gsea {

    tag "GSEA"

    input:
    tuple path(deg_dir), val(comparisons_json), val(de_method), path(script_dir)

    output:
    path "gsea_results", emit: gsea_out

    script:
    """
    mkdir -p gsea_results
    Rscript ${script_dir}/gsea.R \\
        --deg_dir ${deg_dir} \\
        --out_dir gsea_results/ \\
        --script_dir ${script_dir} \\
        --comparisons '${comparisons_json}' \\
        --de_method '${de_method}' \\
        --organism '${params.organism}'
    """

    stub:
    """
    mkdir -p gsea_results
    touch gsea_results/.stub
    """
}
