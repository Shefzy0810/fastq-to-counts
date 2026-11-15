# Phase 6: Salmon -> tximport -> DESeq2 (dex vs untreated, paired by donor)
suppressMessages({library(tximport); library(DESeq2)})

samples <- read.delim("data/samples.tsv")
samples$donor <- factor(samples$donor)
samples$treatment <- factor(samples$treatment, levels = c("untreated", "dex"))
rownames(samples) <- samples$sample

files <- file.path("quant", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read.delim("ref/tx2gene.tsv")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2],
                ignoreAfterBar = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = samples,
                                design = ~ donor + treatment)
keep <- rowSums(counts(dds) >= 10) >= 4   # expressed in at least 4 samples
dds <- dds[keep, ]
cat("Genes kept after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment", "dex", "untreated"), alpha = 0.05)
summary(res)

genes <- unique(tx2gene[, c("gene_id", "gene_name")])
out <- as.data.frame(res)
out$gene_id <- rownames(out)
out$gene_name <- genes$gene_name[match(out$gene_id, genes$gene_id)]
out <- out[order(out$padj), c("gene_id", "gene_name", "baseMean",
                              "log2FoldChange", "lfcSE", "pvalue", "padj")]

dir.create("results", showWarnings = FALSE)
write.csv(out, "results/deseq2_dex_vs_untreated.csv", row.names = FALSE)
saveRDS(dds, "results/dds.rds")

cat("\nTop 10 genes:\n"); print(head(out, 10), row.names = FALSE)
cat("\nKnown dex-response genes:\n")
print(out[out$gene_name %in% c("DUSP1", "KLF15", "PER1", "FKBP5", "TSC22D3"),
          c("gene_name", "log2FoldChange", "padj")], row.names = FALSE)
