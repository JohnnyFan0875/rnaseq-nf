process de_analysis {

    tag "$method"

    input:
    val quant_files
    path metadata_file
    val comparisons_json
    val method
    val result_dir
    path script_dir
    path ref_dir

    output:
    path "DEG_${method}", emit: deg_results

    publishDir { "${result_dir}/de_analysis/" }, mode: 'copy'

    script:
    def script_name = method == 'edgeR' ? 'edgeR.R' : 'deseq2.R'
    def quant_files_str = quant_files.join(',')

    """
    mkdir -p DEG_${method}
    Rscript ${script_dir}/${script_name} \\
        --quant_files ${quant_files_str} \\
        --metadata ${metadata_file} \\
        --comparisons '${comparisons_json}' \\
        --out_dir DEG_${method} \\
        --ref_dir ${ref_dir} \\
        --script_dir ${script_dir}
    """

}
