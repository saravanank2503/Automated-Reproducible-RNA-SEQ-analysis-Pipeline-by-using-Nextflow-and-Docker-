suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
    library(dplyr)
    library(optparse)
})

# ── Parse arguments ───────────────────────────────────────────────────────────
opts <- parse_args(OptionParser(option_list = list(
    make_option("--counts",      type = "character"),
    make_option("--samplesheet", type = "character"),
    make_option("--condition",   type = "character", default = "condition"),
    make_option("--fdr",         type = "double",    default = 0.05),
    make_option("--lfc",         type = "double",    default = 1.0),
    make_option("--outdir",      type = "character", default = "./")
)))

cat("=== DESeq2 Analysis ===\n")
cat("Counts:    ", opts$counts,      "\n")
cat("Metadata:  ", opts$samplesheet, "\n")
cat("Condition: ", opts$condition,   "\n")
cat("FDR:       ", opts$fdr,         "\n")
cat("LFC:       ", opts$lfc,         "\n")

dir.create(file.path(opts$outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
counts <- read.table(opts$counts, header = TRUE, row.names = 1,
                     sep = "\t", check.names = FALSE)
meta   <- read.csv(opts$samplesheet, header = TRUE, stringsAsFactors = FALSE)
COND   <- opts$condition

cat("Counts matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ── Align samples ─────────────────────────────────────────────────────────────
meta         <- meta[meta$sample %in% colnames(counts), ]
counts       <- counts[, meta$sample, drop = FALSE]
meta[[COND]] <- factor(meta[[COND]])

cat("Conditions:", paste(levels(meta[[COND]]), collapse = " vs "), "\n")
cat("Sample-condition map:\n")
print(meta[, c("sample", COND)])

# ── DESeq2 object ─────────────────────────────────────────────────────────────
dds  <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = meta,
    design    = as.formula(paste0("~", COND))
)
keep <- rowSums(counts(dds) >= 5) >= 1
dds  <- dds[keep, ]
cat("Genes after filter:", nrow(dds), "\n")
dds  <- DESeq(dds)

# ── All pairwise comparisons ──────────────────────────────────────────────────
conditions  <- levels(meta[[COND]])
comparisons <- combn(conditions, 2, simplify = FALSE)

for (comp in comparisons) {
    name <- paste(comp[1], "vs", comp[2], sep = "_")
    cat("\nComparison:", name, "\n")

    res <- results(dds,
                   contrast = c(COND, comp[1], comp[2]),
                   alpha    = opts$fdr)

    tryCatch({
        res <- lfcShrink(dds,
                         contrast = c(COND, comp[1], comp[2]),
                         res      = res,
                         type     = "ashr",
                         quiet    = TRUE)
        cat("LFC shrinkage applied\n")
    }, error = function(e) {
        cat("LFC shrinkage skipped:", e$message, "\n")
    })

    res_df <- as.data.frame(res) %>%
        tibble::rownames_to_column("gene_id") %>%
        arrange(padj) %>%
        mutate(
            significant = !is.na(padj) & padj < opts$fdr &
                          abs(log2FoldChange) >= opts$lfc,
            direction   = case_when(
                log2FoldChange > 0 ~ "UP",
                log2FoldChange < 0 ~ "DOWN",
                TRUE               ~ "NS"
            )
        )

    n_sig <- sum(res_df$significant, na.rm = TRUE)
    cat("Significant DEGs:", n_sig, "\n")

    write.csv(res_df,
              file.path(opts$outdir, paste0("deseq2_results_", name, ".csv")),
              row.names = FALSE)

    # ── Volcano plot ──────────────────────────────────────────────────────────
    vp <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue),
                              color = significant)) +
        geom_point(alpha = 0.5, size = 1.5) +
        scale_color_manual(values = c("grey70", "#e74c3c")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                   color = "grey40") +
        geom_vline(xintercept = c(-opts$lfc, opts$lfc), linetype = "dashed",
                   color = "grey40") +
        labs(title = paste("Volcano:", gsub("_", " ", name)),
             x = "log2 Fold Change", y = "-log10(p-value)") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    if (n_sig > 0) {
        top_g <- head(res_df[res_df$significant, ], 15)
        vp    <- vp + geom_label_repel(data = top_g,
                                        aes(label = gene_id),
                                        size = 2.5, max.overlaps = 15)
    }
    ggsave(file.path(opts$outdir, "plots",
                     paste0("volcano_", name, ".pdf")),
           vp, width = 10, height = 8)

    # ── MA plot ───────────────────────────────────────────────────────────────
    res_raw <- results(dds, contrast = c(COND, comp[1], comp[2]))
    pdf(file.path(opts$outdir, "plots", paste0("maplot_", name, ".pdf")),
        width = 8, height = 6)
    plotMA(res_raw, alpha = opts$fdr,
           main = paste("MA Plot:", gsub("_", " ", name)),
           ylim = c(-6, 6))
    abline(h = c(-opts$lfc, opts$lfc), lty = 2, col = "blue")
    dev.off()
    cat("MA plot saved:", name, "\n")
}

# ── Normalized counts ─────────────────────────────────────────────────────────
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts,
          file.path(opts$outdir, "deseq2_normalized_counts.csv"))

# ── PCA plot ──────────────────────────────────────────────────────────────────
tryCatch({
    vsd      <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = COND, returnData = TRUE)
    pvar     <- round(100 * attr(pca_data, "percentVar"))

    pca_p <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
        geom_point(size = 4) +
        geom_text_repel(size = 3) +
        xlab(paste0("PC1: ", pvar[1], "% variance")) +
        ylab(paste0("PC2: ", pvar[2], "% variance")) +
        labs(title = "PCA — All Samples") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(file.path(opts$outdir, "plots", "pca.pdf"),
           pca_p, width = 8, height = 6)
    cat("PCA plot saved\n")
}, error = function(e) {
    cat("PCA skipped:", e$message, "\n")
})

# ── Sample distance heatmap ───────────────────────────────────────────────────
tryCatch({
    vsd       <- vst(dds, blind = FALSE)
    samp_dist <- dist(t(assay(vsd)))
    samp_mat  <- as.matrix(samp_dist)
    annot_df  <- data.frame(row.names = colnames(vsd),
                             condition = vsd[[COND]])
    pdf(file.path(opts$outdir, "plots", "sample_distance_heatmap.pdf"),
        width = 10, height = 9)
    pheatmap(samp_mat,
             clustering_distance_rows = samp_dist,
             clustering_distance_cols = samp_dist,
             annotation_col = annot_df,
             color = colorRampPalette(
                 rev(RColorBrewer::brewer.pal(9, "Blues")))(255),
             main = "Sample-to-Sample Distances")
    dev.off()
    cat("Sample distance heatmap saved\n")
}, error = function(e) cat("Sample distance error:", e$message, "\n"))

# ── Top DEG heatmap ───────────────────────────────────────────────────────────
tryCatch({
    vsd     <- vst(dds, blind = FALSE)
    all_sig <- do.call(rbind, lapply(comparisons, function(comp) {
        name <- paste(comp[1], "vs", comp[2], sep = "_")
        f    <- file.path(opts$outdir,
                          paste0("deseq2_results_", name, ".csv"))
        d    <- read.csv(f)
        d[!is.na(d$significant) & d$significant == TRUE, ]
    }))
    top50 <- head(unique(all_sig$gene_id), 50)
    if (length(top50) >= 2) {
        mat      <- assay(vsd)[top50, ]
        mat      <- t(scale(t(mat)))
        annot_df <- data.frame(row.names = colnames(vsd),
                                condition = vsd[[COND]])
        pdf(file.path(opts$outdir, "plots", "top50_degs_heatmap.pdf"),
            width = 12, height = 14)
        pheatmap(mat,
                 annotation_col = annot_df,
                 fontsize_row   = 7,
                 color = colorRampPalette(
                     c("#2166ac", "white", "#d6604d"))(100),
                 main = paste("Top", length(top50), "DEGs (Z-score)"))
        dev.off()
        cat("Top DEG heatmap saved\n")
    }
}, error = function(e) cat("Heatmap error:", e$message, "\n"))

# ── Dispersion plot ───────────────────────────────────────────────────────────
pdf(file.path(opts$outdir, "plots", "dispersion.pdf"), width = 8, height = 6)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# ── Session info ──────────────────────────────────────────────────────────────
sink(file.path(opts$outdir, "deseq2_session_info.txt"))
cat("DESeq2 Analysis Complete\n")
cat("FDR:", opts$fdr, "\n")
cat("LFC:", opts$lfc, "\n")
cat("Genes tested:", nrow(dds), "\n")
sessionInfo()
sink()

cat("\nDESeq2 complete. Output:", opts$outdir, "\n")
