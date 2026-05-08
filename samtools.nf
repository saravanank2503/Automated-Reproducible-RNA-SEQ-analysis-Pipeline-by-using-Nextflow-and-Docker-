process SAMTOOLS_SORT_INDEX {
    tag "${meta.id}"
    publishDir "${params.outdir}/bam", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path "${meta.id}.flagstat.txt", emit: flagstat

    script:
    """
    samtools sort  -@ ${task.cpus} -o ${meta.id}.sorted.bam ${bam}
    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    samtools flagstat ${meta.id}.sorted.bam > ${meta.id}.flagstat.txt

    # Output contract
    [ ! -f "${meta.id}.sorted.bam.bai" ] && echo "CONTRACT FAIL: BAI index missing" && exit 2
    MAPPED=\$(samtools view -c -F 4 ${meta.id}.sorted.bam)
    [ "\$MAPPED" -lt 10000 ] && echo "CONTRACT FAIL: only \$MAPPED mapped reads" && exit 2
    echo "CONTRACT PASS: \$MAPPED mapped reads, BAI index present"
    """
}
