process multiqc {

    tag "multiqc"

    input:
    path qc_reports

    output:
    path "multiqc_report.html", emit: multiqc_report
    path "multiqc_data",        emit: multiqc_data

    script:
    """
    multiqc $qc_reports -o . --force
    """
}
