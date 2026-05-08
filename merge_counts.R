# bin/merge_counts.R — pure R, no Nextflow syntax
# Merges all featureCounts output files into one matrix

files <- list.files(pattern = "\\.counts\\.txt$", full.names = TRUE)
cat("Merging", length(files), "count files\n")

read_fc <- function(f) {
    d      <- read.table(f, header = TRUE, skip = 1, sep = "\t", check.names = FALSE)
    sample <- gsub("\\.counts\\.txt$", "", basename(f))
    data.frame(gene_id = d[["Geneid"]], count = d[[ncol(d)]],
               stringsAsFactors = FALSE) |>
        setNames(c("gene_id", sample))
}

merged <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE),
                 lapply(files, read_fc))

merged[is.na(merged)] <- 0
rownames(merged)       <- merged$gene_id
merged$gene_id         <- NULL

keep   <- rowSums(merged) > 0
merged <- merged[keep, ]

write.table(merged, "merged_counts.txt", sep = "\t",
            quote = FALSE, col.names = NA)

cat("Matrix:", nrow(merged), "genes x", ncol(merged), "samples\n")
cat("Samples:", paste(colnames(merged), collapse = ", "), "\n")
