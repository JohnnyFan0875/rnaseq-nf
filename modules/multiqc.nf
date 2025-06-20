process multiqc {

    tag "multiqc"

    input:
    path qc_reports
    val result_dir

    output:
    path "multiqc_report.html", emit: multiqc_report
    path "multiqc_data", emit: multiqc_data

    publishDir  { "${result_dir}/multiqc" }, mode: 'copy'

    script:
    """
    TMPDIR=\$(mktemp -d)
    multiqc $qc_reports -o \$TMPDIR --force
    cp -r \$TMPDIR/* .
    """

}
