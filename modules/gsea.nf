process gsea {

    tag "GSEA"

    // publishDir is managed centrally in conf/modules.config

    input:
    // FIX: Removed path(metadata_file) from input tuple — it was staged into
    // the work directory on every run but never passed to gsea.R nor consumed
    // by the script. Removing it eliminates unnecessary I/O and avoids the
    // misleading impression that metadata drives GSEA logic.
    // If metadata is needed in future, re-add it here and in gsea.R.
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
}
