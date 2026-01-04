#===============================================================================
# 步骤2：差异分析
# 分析MDH2、SIRT5、铁死亡基因在心衰中的表达变化
#===============================================================================

library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)

#======================== 参数设置（需要修改）========================
geoID <- "GSE59867"
workDir <- "./02_diff"
dataDir <- "./01_data"

# 差异分析阈值
logFCfilter <- 0.585   # logFC阈值（0.585=1.5倍，1=2倍）
pvalFilter <- 0.05     # 校正P值阈值

# 目标基因（你关注的基因）
targetGenes <- c("MDH2", "SIRT5")

# 铁死亡相关基因（核心基因）
ferroptosisGenes <- c(
  # 铁死亡促进基因
  "ACSL4", "LPCAT3", "ALOX15", "TFRC", "SLC7A11", "GPX4",
  "NCOA4", "FTH1", "FTL", "IREB2", "HMOX1", "NFS1",
  # 铁死亡抑制基因
  "SLC3A2", "GSS", "GCLC", "GCLM", "NFE2L2", "KEAP1",
  # 线粒体相关
  "VDAC2", "VDAC3", "CISD1", "CS", "ACO2"
)

#======================== 创建目录 ========================
if (!dir.exists(workDir)) {
  dir.create(workDir, recursive = TRUE)
}
setwd(workDir)

#======================== 加载数据 ========================
cat("正在加载数据...\n")
load(file.path(dataDir, paste0(geoID, "_rawdata.RData")))

#======================== 数据预处理 ========================
# 1. 探针到基因的映射
cat("正在进行探针注释...\n")

# 根据不同平台调整（这里以常见格式为例）
if ("Gene Symbol" %in% colnames(annot)) {
  probe2gene <- annot[, c("ID", "Gene Symbol")]
  colnames(probe2gene) <- c("probe", "gene")
} else if ("GENE_SYMBOL" %in% colnames(annot)) {
  probe2gene <- annot[, c("ID", "GENE_SYMBOL")]
  colnames(probe2gene) <- c("probe", "gene")
} else {
  # 手动查找基因符号列
  cat("请手动指定基因符号列！\n")
  cat("注释列名：\n")
  print(colnames(annot))
  stop("无法自动识别基因符号列")
}

# 去除空基因和NA
probe2gene <- probe2gene[!is.na(probe2gene$gene) & probe2gene$gene != "", ]
probe2gene <- probe2gene[!grepl("///", probe2gene$gene), ]  # 去除多基因探针

# 2. 表达矩阵转换为基因矩阵
expMatrix <- as.data.frame(expMatrix)
expMatrix$probe <- rownames(expMatrix)
expMatrix <- merge(probe2gene, expMatrix, by = "probe")
expMatrix <- expMatrix[, -1]  # 删除probe列

# 同基因取平均
geneMatrix <- aggregate(. ~ gene, data = expMatrix, FUN = mean)
rownames(geneMatrix) <- geneMatrix$gene
geneMatrix <- geneMatrix[, -1]
geneMatrix <- as.matrix(geneMatrix)

cat("基因表达矩阵维度：", nrow(geneMatrix), "个基因 x", ncol(geneMatrix), "个样本\n")

# 3. 数据标准化
qx <- as.numeric(quantile(geneMatrix, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0))
if (LogC) {
  geneMatrix[geneMatrix < 0] <- 0
  geneMatrix <- log2(geneMatrix + 1)
  cat("已进行log2转换\n")
}
geneMatrix <- normalizeBetweenArrays(geneMatrix)
cat("已进行分位数标准化\n")

#======================== 设置分组（需要根据实际数据修改）========================
cat("\n请检查样本分组信息...\n")
# 查看phenoData中的分组信息
print(head(phenoData[, c("title", "source_name_ch1")]))

# 根据实际情况修改分组
# 示例：假设有"control"和"heart failure"
# 你需要根据phenoData的实际列来修改这里

# 创建分组向量
# 这里需要根据你的数据修改！
group <- ifelse(grepl("control|normal|healthy", phenoData$source_name_ch1, ignore.case = TRUE),
                "Control", "Disease")
names(group) <- rownames(phenoData)
table(group)

#======================== 差异分析 ========================
cat("\n正在进行差异分析...\n")

# 确保样本顺序一致
geneMatrix <- geneMatrix[, names(group)]

# 构建设计矩阵
group <- factor(group, levels = c("Control", "Disease"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 拟合线性模型
fit <- lmFit(geneMatrix, design)
contrast.matrix <- makeContrasts(Disease - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 获取所有基因的结果
allDiff <- topTable(fit2, adjust = "fdr", number = Inf)
allDiff$gene <- rownames(allDiff)

# 筛选差异基因
diffSig <- allDiff[abs(allDiff$logFC) > logFCfilter & allDiff$adj.P.Val < pvalFilter, ]
cat("差异基因数量：", nrow(diffSig), "\n")
cat("  上调基因：", sum(diffSig$logFC > 0), "\n")
cat("  下调基因：", sum(diffSig$logFC < 0), "\n")

#======================== 保存差异分析结果 ========================
write.table(allDiff, file = "all_genes_diff.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diffSig, file = "sig_genes_diff.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#======================== 目标基因表达情况 ========================
cat("\n========== 目标基因表达情况 ==========\n")
allGenes <- c(targetGenes, ferroptosisGenes)
targetResult <- allDiff[allDiff$gene %in% allGenes, ]
targetResult <- targetResult[order(targetResult$adj.P.Val), ]
print(targetResult[, c("gene", "logFC", "P.Value", "adj.P.Val")])

# 保存目标基因结果
write.table(targetResult, file = "target_genes_result.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#======================== MDH2和SIRT5详细信息 ========================
cat("\n========== MDH2 和 SIRT5 详细信息 ==========\n")
for (gene in c("MDH2", "SIRT5")) {
  if (gene %in% rownames(allDiff)) {
    info <- allDiff[gene, ]
    cat(paste0("\n", gene, ":\n"))
    cat(paste0("  logFC = ", round(info$logFC, 3), "\n"))
    cat(paste0("  P.Value = ", format(info$P.Value, digits = 3), "\n"))
    cat(paste0("  adj.P.Val = ", format(info$adj.P.Val, digits = 3), "\n"))
    if (info$logFC > 0) {
      cat("  → 在心衰中上调\n")
    } else {
      cat("  → 在心衰中下调\n")
    }
  } else {
    cat(paste0("\n", gene, ": 未在数据中找到\n"))
  }
}

#======================== 绑制火山图 ========================
cat("\n正在绑制火山图...\n")

allDiff$Significance <- "Not Sig"
allDiff$Significance[allDiff$logFC > logFCfilter & allDiff$adj.P.Val < pvalFilter] <- "Up"
allDiff$Significance[allDiff$logFC < -logFCfilter & allDiff$adj.P.Val < pvalFilter] <- "Down"

# 标注目标基因
allDiff$label <- ""
allDiff$label[allDiff$gene %in% c("MDH2", "SIRT5")] <- allDiff$gene[allDiff$gene %in% c("MDH2", "SIRT5")]

p1 <- ggplot(allDiff, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Down" = "blue", "Not Sig" = "grey", "Up" = "red")) +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(pvalFilter), linetype = "dashed", color = "grey40") +
  geom_text(aes(label = label), hjust = -0.1, vjust = 0, size = 3, color = "black") +
  labs(title = paste0(geoID, " Differential Expression"),
       x = "log2(Fold Change)", y = "-log10(adj.P.Value)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf("volcano_plot.pdf", width = 8, height = 6)
print(p1)
dev.off()

#======================== 目标基因热图 ========================
cat("正在绑制目标基因热图...\n")

# 提取目标基因表达
targetExp <- geneMatrix[rownames(geneMatrix) %in% allGenes, ]
targetExp <- targetExp[rowSums(is.na(targetExp)) == 0, ]

if (nrow(targetExp) > 0) {
  # 准备注释
  annotation_col <- data.frame(
    Group = group,
    row.names = names(group)
  )

  pdf("target_genes_heatmap.pdf", width = 12, height = 8)
  pheatmap(targetExp,
           scale = "row",
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           annotation_col = annotation_col,
           show_colnames = FALSE,
           main = "Target Genes Expression (MDH2, SIRT5, Ferroptosis)")
  dev.off()
}

#======================== 保存处理后的数据 ========================
save(geneMatrix, group, allDiff, diffSig, targetResult,
     file = "diff_analysis_results.RData")

cat("\n========== 差异分析完成 ==========\n")
cat("输出文件：\n")
cat("  - all_genes_diff.txt（所有基因差异结果）\n")
cat("  - sig_genes_diff.txt（显著差异基因）\n")
cat("  - target_genes_result.txt（目标基因结果）\n")
cat("  - volcano_plot.pdf（火山图）\n")
cat("  - target_genes_heatmap.pdf（目标基因热图）\n")
cat("  - diff_analysis_results.RData（R数据）\n")
