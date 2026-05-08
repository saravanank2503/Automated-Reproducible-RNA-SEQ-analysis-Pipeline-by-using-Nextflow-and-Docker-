process TRIMGALORE {
    tag "${meta.id}"
    publishDir "${params.outdir}/trimgalore", mode: params.publish_dir_mode,
        saveAs: { f -> f.endsWith('.txt') ? "logs/$f" : params.save_trimmed ? "trimmed/$f" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fq.gz"), emit: reads
    path "*_trimming_report.txt",              emit: log

    script:
    def cores = Math.min(task.cpus, 4)
    if (meta.paired)
    """
    trim_galore \
        --paired \
        --quality 20 \
        --length 36 \
        --adapter  AGATCGGAAGAGC \
        --adapter2 AGATCGGAAGAGC \
        --cores ${cores} \
        --gzip \
        --basename ${meta.id} \
        ${reads[0]} ${reads[1]}

    mv ${meta.id}_val_1.fq.gz ${meta.id}_R1_trimmed.fq.gz
    mv ${meta.id}_val_2.fq.gz ${meta.id}_R2_trimmed.fq.gz

    COUNT=\$(gzip -dc ${meta.id}_R1_trimmed.fq.gz | awk 'NR%4==1' | wc -l)
    [ "\$COUNT" -lt 1000 ] && echo "CONTRACT FAIL: only \$COUNT reads" && exit 2
    echo "CONTRACT PASS: \$COUNT pairs after trimming"
    """
    else
    """
    trim_galore \
        --quality 20 \
        --length 36 \
        --adapter AGATCGGAAGAGC \
        --cores ${cores} \
        --gzip \
        --basename ${meta.id} \
        ${reads[0]}

    mv ${meta.id}_trimmed.fq.gz ${meta.id}_R1_trimmed.fq.gz

    COUNT=\$(gzip -dc ${meta.id}_R1_trimmed.fq.gz | awk 'NR%4==1' | wc -l)
    [ "\$COUNT" -lt 1000 ] && echo "CONTRACT FAIL: only \$COUNT reads" && exit 2
    echo "CONTRACT PASS: \$COUNT reads after trimming"
    """
}
