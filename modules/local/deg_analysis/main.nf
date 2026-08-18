process deg_analysis {

    tag "$method"

    input:
    path expression_inputs
    path metadata_file
    val  comparisons_json
    val  method
    path script_dir
    path gtf_file
    val  input_mode
    val  quant_type

    output:
    path "DEG_${method}", emit: deg_results

    script:
    def script_name = method == 'edgeR' ? 'edgeR.R' : 'deseq2.R'
    def expression_pattern = quant_type == 'salmon' ? 'quant.sf' : 'abundance.tsv'

    """
    mkdir -p DEG_${method}
    Rscript ${script_dir}/${script_name} \\
        --input_mode ${input_mode} \\
        --quant_type ${quant_type} \\
        --expression_inputs "\$(find -L . -name '${expression_pattern}' | sort | paste -sd ',')" \\
        --metadata ${metadata_file} \\
        --comparisons '${comparisons_json}' \\
        --out_dir DEG_${method} \\
        --gtf ${gtf_file} \\
        --script_dir ${script_dir}
    """

    stub:
    """
    mkdir -p DEG_${method}
    touch DEG_${method}/stub.txt
    """
}
