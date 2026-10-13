# ================================
# 拟南芥 WT DMSO vs NAE 18:3 差异分析
# ================================

library(DESeq2)

# 1. 读 featureCounts 输出
counts_raw <- read.table("~/AA/counts/gene_counts.txt",
                         header=TRUE, row.names=1, skip=1, sep="\t",
                         check.names=FALSE)

# 只保留样本列（去掉 Chr/Start/End/Strand/Length）
counts <- counts_raw[, 6:ncol(counts_raw)]

# 重命名列
colnames(counts) <- c("SRR10435452","SRR10435453","SRR10435454",
                      "SRR10435455","SRR10435456","SRR10435457")

# 2. 分组信息
coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(c("DMSO","DMSO","DMSO",
                       "NAE","NAE","NAE"),
                     levels = c("DMSO","NAE"))
)
print(coldata)

# 3. 构建 DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# 4. 预过滤：至少 3 个样本 counts >= 10
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("过滤后基因数:", nrow(dds), "\n")

# 5. 差异分析
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition","NAE","DMSO"))
res <- res[order(res$padj), ]

# 6. 导出结果
write.csv(as.data.frame(res), "~/AA/results/DESeq2_results.csv")

# 7. 显著差异基因
res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
cat("显著差异基因数 (padj<0.05 & |log2FC|>1):", nrow(res_sig), "\n")
write.csv(as.data.frame(res_sig), "~/AA/results/DESeq2_sig_genes.csv")

# 8. 标准化表达矩阵
norm_counts <- counts(dds, normalized=TRUE)
write.csv(norm_counts, "~/AA/results/normalized_counts.csv")

# 9. 样本聚类 / PCA 图
vsd <- vst(dds, blind=FALSE)

pdf("~/AA/results/PCA.pdf", width=6, height=5)
plotPCA(vsd, intgroup="condition") + ggplot2::theme_bw()
dev.off()

# --- 用 pheatmap 替换 heatmap ---
library(pheatmap)

pdf("~/AA/results/sample_distance_heatmap.pdf", width=9, height=7)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, colnames(vsd), sep="-")
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,
         main = "Sample distance",
         fontsize = 10,
         fontsize_row = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
dev.off()
# --- 替换结束 ---

# 10. MA 图、火山图
pdf("~/AA/results/MA_plot.pdf", width=6, height=5)
plotMA(res, ylim=c(-5,5), main="MA plot NAE vs DMSO")
dev.off()

pdf("~/AA/results/volcano.pdf", width=6, height=5)
plot(res$log2FoldChange, -log10(res$padj),
     pch=20, cex=0.5,
     xlab="log2 Fold Change (NAE/DMSO)",
     ylab="-log10(padj)",
     main="Volcano plot",
     col=ifelse(res$padj<0.05 & abs(res$log2FoldChange)>1, "red", "grey"))
abline(h=-log10(0.05), v=c(-1,1), lty=2, col="blue")
dev.off()

cat("分析完成，结果在 ~/AA/results/\n")
