process FASTQC {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads}

    # Output contract
    ZIP_COUNT=\$(ls *.zip 2>/dev/null | wc -l)
    [ "\$ZIP_COUNT" -eq 0 ] && echo "CONTRACT FAIL: no FastQC zip produced" && exit 2
    echo "CONTRACT PASS: \$ZIP_COUNT FastQC zip(s) produced"
    """
}
