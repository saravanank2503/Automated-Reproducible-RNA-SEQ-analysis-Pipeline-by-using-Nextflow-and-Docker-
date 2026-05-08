process FEATURECOUNTS {
    tag "${meta.id}"
    publishDir "${params.outdir}/featurecounts", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path  gtf

    output:
    tuple val(meta), path("${meta.id}.counts.txt"),  emit: counts
    path "${meta.id}.counts.txt.summary",             emit: summary

    script:
    def paired = meta.paired ? "-p" : ""
    """
    featureCounts \
        -T ${task.cpus} \
        -a ${gtf} \
        -o ${meta.id}.counts.txt \
        ${paired} \
        -s 2 \
        -t exon \
        -g gene_id \
        --minOverlap 10 \
        -Q 10 \
        ${bam}

    ASSIGNED=\$(grep "Assigned" ${meta.id}.counts.txt.summary | awk '{print \$2}')
    echo "Assigned reads: \$ASSIGNED"
    [ ! -s "${meta.id}.counts.txt" ] && echo "CONTRACT FAIL: counts empty" && exit 2
    echo "CONTRACT PASS: \$ASSIGNED reads assigned"
    """
}
