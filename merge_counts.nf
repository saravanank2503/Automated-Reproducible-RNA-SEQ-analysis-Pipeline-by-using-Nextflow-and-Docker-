process MERGE_COUNTS {
    tag "merging counts"
    publishDir "${params.outdir}/counts", mode: params.publish_dir_mode

    input:
    path count_files

    output:
    path "merged_counts.txt", emit: matrix

    script:
    """
    Rscript ${projectDir}/bin/merge_counts.R

    [ ! -s merged_counts.txt ] && echo "CONTRACT FAIL: merged matrix empty" && exit 2
    GENES=\$(wc -l < merged_counts.txt)
    echo "CONTRACT PASS: \$GENES genes in merged matrix"
    """
}
