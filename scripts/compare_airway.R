# Phase 7: compare our Salmon/DESeq2 results with the published airway counts
suppressMessages({library(DESeq2); library(airway); library(ggplot2)})

data("airway")
airway$dex <- relevel(airway$dex, "untrt")
dds_pub <- DESeqDataSet(airway, design = ~ cell + dex)
dds_pub <- dds_pub[rowSums(counts(dds_pub) >= 10) >= 4, ]
dds_pub <- DESeq(dds_pub, quiet = TRUE)
pub <- as.data.frame(results(dds_pub, contrast = c("dex", "trt", "untrt"), alpha = 0.05))
pub$gene <- rownames(pub)

ours <- read.csv("results/deseq2_dex_vs_untreated.csv")
ours$gene <- sub("\\..*", "", ours$gene_id)   # drop GENCODE version (.14)
ours <- ours[!duplicated(ours$gene), ]

m <- merge(ours, pub, by = "gene", suffixes = c(".ours", ".pub"))
m <- m[!is.na(m$padj.ours) & !is.na(m$padj.pub), ]
so <- m$padj.ours < 0.05; sp <- m$padj.pub < 0.05

out <- c(
  paste("Genes tested in both:", nrow(m)),
  paste("Significant (padj<0.05) ours:", sum(so), "| published:", sum(sp)),
  paste("Significant in both:", sum(so & sp)),
  paste0("Our hits confirmed by published: ", round(100 * sum(so & sp) / sum(so), 1), "%"),
  paste0("Published hits recovered by ours: ", round(100 * sum(so & sp) / sum(sp), 1), "%"),
  paste("Same direction among shared hits:",
        sum(sign(m$log2FoldChange.ours[so & sp]) == sign(m$log2FoldChange.pub[so & sp])),
        "of", sum(so & sp)),
  paste("Pearson r of log2FC (all shared genes):",
        round(cor(m$log2FoldChange.ours, m$log2FoldChange.pub), 3)))
cat(out, sep = "\n"); writeLines(out, "results/compare_airway.txt")

p <- ggplot(m, aes(log2FoldChange.pub, log2FoldChange.ours)) +
  geom_point(size = 0.6, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, colour = "firebrick") + theme_bw() +
  labs(x = "log2FC, published airway counts (full depth)",
       y = "log2FC, this pipeline (Salmon, 2M reads)",
       title = "Our results vs published airway")
ggsave("results/compare_airway.png", p, width = 6, height = 5, dpi = 150)
