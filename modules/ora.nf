process ora {

    tag "ORA"

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
    path "ora_results", emit: gsea_out

    publishDir "${result_dir}/", mode: 'copy'

    script:
    """
    mkdir -p ora_results
    Rscript ${script_dir}/ora.R \\
        --deg_dir ${deg_dir} \\
        --metadata ${metadata_file} \\
        --out_dir ora_results/ \\
        --script_dir ${script_dir} \\
        --comparisons '${comparisons_json}' \\
        --de_method '${de_method}'
    """
}
