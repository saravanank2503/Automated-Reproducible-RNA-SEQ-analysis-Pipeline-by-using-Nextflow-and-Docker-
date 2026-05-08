process STAR_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/star/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)
    path  star_index
    path  gtf

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${meta.id}.Log.final.out",                                    emit: log

    script:
    def reads_in = meta.paired ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    """
    STAR \
        --runThreadN       ${task.cpus} \
        --genomeDir        ${star_index} \
        --sjdbGTFfile      ${gtf} \
        --readFilesIn      ${reads_in} \
        --readFilesCommand zcat \
        --outSAMtype       BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFileNamePrefix ${meta.id}. \
        --outFilterMismatchNmax 2 \
        --quantMode        GeneCounts \
        --twopassMode      Basic \
        --limitBAMsortRAM  30000000000

    BAM="${meta.id}.Aligned.sortedByCoord.out.bam"
    [ ! -f "\$BAM" ] && echo "CONTRACT FAIL: BAM not produced" && exit 2
    [ ! -s "\$BAM" ] && echo "CONTRACT FAIL: BAM is empty"     && exit 2

    RATE=\$(grep "Uniquely mapped reads %" ${meta.id}.Log.final.out | awk '{print \$NF}')
    echo "Unique mapping rate: \$RATE"
    echo "CONTRACT PASS: BAM produced, mapping rate \$RATE"
    """
}
