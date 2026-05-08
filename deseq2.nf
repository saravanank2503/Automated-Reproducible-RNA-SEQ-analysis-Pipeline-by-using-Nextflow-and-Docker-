process DESEQ2_ANALYSIS {
    tag "DESeq2"
    publishDir "${params.outdir}/deseq2", mode: params.publish_dir_mode

    input:
    path count_matrix
    path samplesheet
    val  condition_col
    val  fdr
    val  lfc

    output:
    path "deseq2_results*.csv",           emit: results
    path "deseq2_normalized_counts.csv",  emit: normalized
    path "plots/*.pdf",                   emit: plots
    path "deseq2_session_info.txt",       emit: session

    script:
    """
    mkdir -p plots

    Rscript ${projectDir}/bin/deseq2_analysis.R \
        --counts      ${count_matrix} \
        --samplesheet ${samplesheet} \
        --condition   ${condition_col} \
        --fdr         ${fdr} \
        --lfc         ${lfc} \
        --outdir      ./

    [ ! -f deseq2_session_info.txt ] && \
        echo "CONTRACT FAIL: DESeq2 did not complete" && exit 2

    echo "CONTRACT PASS: DESeq2 complete"
    echo "Plots generated:"
    ls plots/*.pdf
    """
}
