# Arabidopsis thaliana RNA-seq Analysis: Response to NAE 18:3

> 拟南芥幼苗响应 NAE 18:3（3h）的 RNA-seq 有参转录组分析。包含从 NCBI SRA 数据下载、质控、比对、定量、差异分析到功能富集的完整流程。

## 📌 实验设计
- **物种**：拟南芥 (*Arabidopsis thaliana*)
- **处理条件**：Control (DMSO, 3h) vs Treatment (NAE 18:3, 3h)
- **生物学重复**：3 vs 3
- **数据来源**：NCBI SRA (SRR10435452 - SRR10435457)

## 🔄 分析流程
1. 数据下载：SRA Toolkit (fasterq-dump)
2. 质控清洗：FastQC → fastp → FastQC → MultiQC
3. 序列比对：HISAT2 → TAIR10 参考基因组
4. 基因定量：featureCounts
5. 差异分析：DESeq2
6. 功能富集：clusterProfiler (GO + KEGG)

## 📊 主要结果
- 共鉴定出 **5567** 个显著差异表达基因（|log2FC| > 1, padj < 0.05）。
- **GO 富集**：显著富集于缺氧响应、免疫响应、茉莉酸响应等生物学过程。
- **KEGG 富集**：植物激素信号转导、植物-病原体相互作用为核心通路。

## 🛠️ 依赖工具
**上游**：WSL2 / Conda / SRA Toolkit / FastQC / fastp / MultiQC / HISAT2 / featureCounts
**下游**：R (4.5.3) / DESeq2 / clusterProfiler / org.At.tair.db / enrichplot / ggplot2

## 📂 目录结构
- `counts/`：featureCounts 定量结果
- `results/`：差异分析结果、PCA/火山图/热图、富集气泡图
- `scripts/`：Shell 与 R 分析脚本
