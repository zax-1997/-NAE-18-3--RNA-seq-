library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
library(ggplot2)
library(stringr)  # 新增：用于长文字换行

setwd("/home/zax/AA/results")

sig <- read.csv("DESeq2_sig_genes.csv", row.names = 1)
gene_list <- rownames(sig)
cat("显著基因数:", length(gene_list), "\n")
head(gene_list)

outdir <- "/home/zax/AA/results/enrich"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (ont in c("BP","MF","CC")) {
  ego <- enrichGO(gene = gene_list,
                  OrgDb = org.At.tair.db,
                  keyType = "TAIR",
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file.path(outdir, paste0("GO_", ont, ".csv")), row.names = FALSE)
    
    # 修改1：加入 str_wrap 对长文字自动换行（每40个字符换行），并缩小字体大小
    p <- dotplot(ego, showCategory = 20, title = paste0("GO ", ont)) +
         scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
         theme(axis.text.y = element_text(size = 9))
    
    # 修改2：大幅增加高度（height）到 15，宽度（width）到 10，解决重叠
    ggsave(file.path(outdir, paste0("GO_", ont, "_dotplot.pdf")), p, 
           width = 10, height = 15, limitsize = FALSE)
    cat("GO", ont, "完成:", nrow(as.data.frame(ego)), "条\n")
  } else {
    cat("GO", ont, "无显著条目\n")
  }
}

kk <- enrichKEGG(gene = gene_list, organism = "ath",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05)
if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
  write.csv(as.data.frame(kk), file.path(outdir, "KEGG.csv"), row.names = FALSE)
  
  # 修改3：KEGG 图同样处理，换行 + 增加高度
  p <- dotplot(kk, showCategory = 20, title = "KEGG Pathway") +
       scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
       theme(axis.text.y = element_text(size = 9))
       
  ggsave(file.path(outdir, "KEGG_dotplot.pdf"), p, 
         width = 10, height = 15, limitsize = FALSE)
  cat("KEGG 完成:", nrow(as.data.frame(kk)), "条\n")
} else {
  cat("KEGG 无显著条目\n")
}

cat("富集分析完成\n")
