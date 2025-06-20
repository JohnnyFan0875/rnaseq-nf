process gsea {

    tag "GSEA"

    input:
    tuple (
        path(deg_dir), 
        path(metadata_file), 
        val(comparisons_json), 
        val(result_dir), 
        val(de_method), 
        path(script_dir) 
    )

    output:
    path "gsea_results", emit: gsea_out

    publishDir "${result_dir}/", mode: 'copy'

    script:
    """
    mkdir -p gsea_results
    Rscript ${script_dir}/gsea.R \\
        --deg_dir ${deg_dir} \\
        --metadata ${metadata_file} \\
        --out_dir gsea_results/ \\
        --script_dir ${script_dir} \\
        --comparisons '${comparisons_json}' \\
        --de_method '${de_method}'
    """
}
