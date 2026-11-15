# Phase 7: PCA, volcano and MA plots
suppressMessages({library(DESeq2); library(ggplot2)})
dds <- readRDS("results/dds.rds")
res <- read.csv("results/deseq2_dex_vs_untreated.csv")

# PCA on variance-stabilised counts
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("treatment", "donor"), returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
p1 <- ggplot(pca, aes(PC1, PC2, colour = treatment, shape = donor)) +
  geom_point(size = 4) + theme_bw() +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
       title = "PCA of airway samples (VST counts)")
ggsave("results/pca.png", p1, width = 6, height = 4.5, dpi = 150)

# Volcano plot
res <- res[!is.na(res$padj), ]
res$sig <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "padj<0.05 & |LFC|>1", "not sig")
known <- res[res$gene_name %in% c("DUSP1", "KLF15", "PER1", "FKBP5", "TSC22D3"), ]
p2 <- ggplot(res, aes(log2FoldChange, -log10(padj), colour = sig)) +
  geom_point(size = 0.8, alpha = 0.6) + theme_bw() +
  scale_colour_manual(values = c("grey60", "firebrick")) +
  geom_text(data = known, aes(label = gene_name), colour = "black",
            vjust = -0.7, size = 3.5) +
  labs(title = "Dex vs untreated", colour = NULL)
ggsave("results/volcano.png", p2, width = 6.5, height = 5, dpi = 150)

# MA plot
png("results/ma_plot.png", width = 900, height = 700, res = 150)
plotMA(results(dds, contrast = c("treatment", "dex", "untreated")), ylim = c(-5, 5))
dev.off()
cat("Saved: results/pca.png, results/volcano.png, results/ma_plot.png\n")
