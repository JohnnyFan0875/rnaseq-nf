process de_analysis {

    tag "$method"

    input:
    val  sample_ids
    path quant_dirs
    path metadata_file
    val  comparisons_json
    val  method
    val  result_dir
    path script_dir
    path ref_dir

    output:
    path "DEG_${method}", emit: deg_results

    publishDir { "${result_dir}/de_analysis/" }, mode: 'copy'

    script:
    def script_name = method == 'edgeR' ? 'edgeR.R' : 'deseq2.R'
    def ids = sample_ids instanceof List ? sample_ids : [sample_ids]

    """
    mkdir -p DEG_${method}
    QUANT_FILES=\$(for id in ${ids.join(' ')}; do echo "\$(pwd)/\${id}/abundance.tsv"; done | paste -sd ',')

    Rscript ${script_dir}/${script_name} \\
        --quant_files \$QUANT_FILES \\
        --metadata ${metadata_file} \\
        --comparisons '${comparisons_json}' \\
        --out_dir DEG_${method} \\
        --ref_dir ${ref_dir} \\
        --script_dir ${script_dir}
    """
}