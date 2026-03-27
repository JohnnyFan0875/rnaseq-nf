process de_analysis {

    tag "$method"

    input:
    val  sample_ids
    path quant_dirs
    path metadata_file
    val  comparisons_json
    val  method
    path script_dir
    path ref_dir

    output:
    path "DEG_${method}", emit: deg_results

    script:
    def script_name = method == 'edgeR' ? 'edgeR.R' : 'deseq2.R'
    def ids = sample_ids instanceof List ? sample_ids : [sample_ids]

    """
    mkdir -p DEG_${method}

    QUANT_FILES=\$(find -L . -name "abundance.tsv" | sort | paste -sd ',')

    if [[ -z "\$QUANT_FILES" ]]; then
        echo "ERROR: No abundance.tsv files found in staged quant directories." >&2
        exit 1
    fi

    Rscript ${script_dir}/${script_name} \\
        --quant_files \$QUANT_FILES \\
        --metadata ${metadata_file} \\
        --comparisons '${comparisons_json}' \\
        --out_dir DEG_${method} \\
        --ref_dir ${ref_dir} \\
        --script_dir ${script_dir}
    """
}
