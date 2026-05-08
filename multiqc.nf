process MULTIQC {
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path qc_files

    output:
    path "*multiqc_report.html", emit: report

    script:
    """
    multiqc --force --title "FH vs Healthy RNA-seq QC" -o . .

    [ ! -f *multiqc_report.html ] && echo "CONTRACT FAIL: MultiQC report missing" && exit 2
    echo "CONTRACT PASS: MultiQC report generated"
    """
}
