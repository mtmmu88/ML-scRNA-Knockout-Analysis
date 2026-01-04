#===============================================================================
# MDH2心衰研究 - 测试版本
# 步骤3：差异分析
#
# 参考原代码：09.diff/geo09.diff.R
# 重点关注：MDH2、SIRT5、铁死亡基因
#===============================================================================

#加载包
library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)

#======================== 参数设置（可修改）========================
geoID <- "GSE26887"           #GEO数据集ID
logFCfilter <- 0.585          #logFC阈值（0.585=1.5倍，1=2倍）
adj.P.Val.Filter <- 0.05      #校正P值阈值
diseaseName <- "HeartFailure" #疾病名称（用于图例）

#目标基因（你关注的基因）
targetGenes <- c("MDH2", "SIRT5")

#铁死亡相关基因
ferroptosisGenes <- c(
  #铁死亡促进基因
  "ACSL4", "LPCAT3", "ALOX15", "TFRC", "NCOA4",
  #铁死亡抑制基因
  "SLC7A11", "GPX4", "FTH1", "FTL", "NFE2L2",
  #线粒体相关
  "VDAC2", "VDAC3", "CS", "ACO2", "HMOX1"
)

#设置工作目录（修改为你的路径）
setwd("D:/test_MDH2")         #Windows用户修改这里
#setwd("~/test_MDH2")         #Mac/Linux用户用这个

#======================== 加载数据 ========================
cat("正在加载标准化数据...\n")
inputFile <- paste0(geoID, ".normalize.txt")

#读取表达矩阵
rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

cat("基因数：", nrow(data), "\n")
cat("样本数：", ncol(data), "\n")

#======================== 获取分组信息 ========================
#从列名提取分组（格式：样本名_分组）
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, order(Type)]  #按分组排序
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data) <- gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

cat("\n分组情况：\n")
print(table(Type))

#======================== 差异分析 ========================
cat("\n正在进行差异分析...\n")

#构建设计矩阵
design <- model.matrix(~ 0 + factor(Type))
colnames(design) <- c("Control", "Treat")

#拟合线性模型
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(Treat - Control, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#获取所有差异结果
allDiff <- topTable(fit2, adjust = 'fdr', number = 200000)
allDiff$gene <- rownames(allDiff)
allDiffOut <- rbind(id = colnames(allDiff), allDiff)
write.table(allDiffOut, file = "diff.bindao_all.txt", sep = "\t", quote = FALSE, col.names = FALSE)

#筛选显著差异基因
diffSig <- allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
diffSigOut <- rbind(id = colnames(diffSig), diffSig)
write.table(diffSigOut, file = "diff.bindao_sig.txt", sep = "\t", quote = FALSE, col.names = FALSE)

cat("\n差异分析结果：\n")
cat("  显著差异基因：", nrow(diffSig), "个\n")
cat("  上调基因：", sum(diffSig$logFC > 0), "个\n")
cat("  下调基因：", sum(diffSig$logFC < 0), "个\n")

#======================== 目标基因结果 ========================
cat("\n========== MDH2和SIRT5表达情况 ==========\n")
for (gene in c("MDH2", "SIRT5")) {
  if (gene %in% rownames(allDiff)) {
    info <- allDiff[gene, ]
    cat(paste0("\n", gene, ":\n"))
    cat(paste0("  logFC = ", round(info$logFC, 4), "\n"))
    cat(paste0("  P.Value = ", format(info$P.Value, digits = 4), "\n"))
    cat(paste0("  adj.P.Val = ", format(info$adj.P.Val, digits = 4), "\n"))
    if (info$logFC > 0) {
      cat("  结论：在心衰中【上调】\n")
    } else {
      cat("  结论：在心衰中【下调】\n")
    }
    if (info$adj.P.Val < 0.05) {
      cat("  显著性：【显著】\n")
    } else {
      cat("  显著性：不显著\n")
    }
  } else {
    cat(paste0("\n", gene, ": 未在数据中找到\n"))
  }
}

#======================== 铁死亡基因结果 ========================
cat("\n========== 铁死亡相关基因 ==========\n")
allGenes <- c(targetGenes, ferroptosisGenes)
targetResult <- allDiff[rownames(allDiff) %in% allGenes, ]
targetResult <- targetResult[order(targetResult$P.Value), ]
print(targetResult[, c("logFC", "P.Value", "adj.P.Val")])

#保存目标基因结果
write.table(cbind(gene = rownames(targetResult), targetResult),
            file = "diff.bindao_targetGenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#======================== 绑制火山图 ========================
cat("\n正在绑制火山图...\n")

Sig <- ifelse((allDiff$adj.P.Val < adj.P.Val.Filter) & (abs(allDiff$logFC) > logFCfilter),
              ifelse(allDiff$logFC > logFCfilter, "Up", "Down"), "Not")
allDiff$Sig <- Sig

#标记目标基因
allDiff$label <- ""
allDiff$label[rownames(allDiff) %in% c("MDH2", "SIRT5")] <-
  rownames(allDiff)[rownames(allDiff) %in% c("MDH2", "SIRT5")]

p <- ggplot(allDiff, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(col = Sig), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down" = "green3", "Not" = "grey", "Up" = "red3")) +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(adj.P.Val.Filter), linetype = "dashed", color = "grey40") +
  geom_text(aes(label = label), hjust = -0.2, vjust = 0, size = 3, fontface = "bold") +
  labs(title = paste0(geoID, " - Differential Expression"),
       x = "log2(Fold Change)", y = "-log10(adj.P.Value)") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

pdf(file = "bindao_bindao_volcano.pdf", width = 7, height = 6)
print(p)
dev.off()

cat("火山图已保存：bindao_bindao_volcano.pdf\n")

#======================== 绑制热图 ========================
cat("正在绘制目标基因热图...\n")

#提取目标基因表达
hmGenes <- allGenes[allGenes %in% rownames(data)]
if (length(hmGenes) > 0) {
  hmExp <- data[hmGenes, ]

  #准备注释
  Type_plot <- ifelse(Type == "Control", "Control", diseaseName)
  Type_plot <- factor(Type_plot, levels = c("Control", diseaseName))
  annotation_col <- data.frame(Group = Type_plot, row.names = colnames(data))

  pdf(file = "bindao_bindao_heatmap_target.pdf", width = 10, height = 8)
  pheatmap(hmExp,
           annotation_col = annotation_col,
           color = colorRampPalette(c("blue2", "white", "red2"))(50),
           cluster_cols = FALSE,
           show_colnames = FALSE,
           scale = "row",
           fontsize = 10,
           fontsize_row = 8,
           main = "MDH2 / SIRT5 / Ferroptosis Genes")
  dev.off()

  cat("热图已保存：bindao_bindao_heatmap_target.pdf\n")
}

#======================== 差异基因热图（Top50）========================
geneNum <- 50
diffUp <- diffSig[diffSig$logFC > 0, ]
diffDown <- diffSig[diffSig$logFC < 0, ]
geneUp <- rownames(diffUp)
geneDown <- rownames(diffDown)
if (nrow(diffUp) > geneNum) geneUp <- rownames(diffUp)[1:geneNum]
if (nrow(diffDown) > geneNum) geneDown <- rownames(diffDown)[1:geneNum]

if (length(c(geneUp, geneDown)) > 0) {
  hmExp <- data[c(geneUp, geneDown), ]

  pdf(file = "bindao_bindao_heatmap_diff.pdf", width = 10, height = 10)
  pheatmap(hmExp,
           annotation_col = annotation_col,
           color = colorRampPalette(c("blue2", "white", "red2"))(50),
           cluster_cols = FALSE,
           show_colnames = FALSE,
           scale = "row",
           fontsize = 8,
           fontsize_row = 5,
           main = paste0("Top ", geneNum, " Up/Down Genes"))
  dev.off()

  cat("差异基因热图已保存：bindao_bindao_heatmap_diff.pdf\n")
}

#======================== 保存结果 ========================
save(allDiff, diffSig, targetResult, data, Type,
     file = paste0(geoID, ".bindao_diff.RData"))

cat("\n========== 差异分析完成 ==========\n")
cat("保存文件：\n")
cat("  - diff.bindao_all.txt（所有基因结果）\n")
cat("  - diff.bindao_sig.txt（显著差异基因）\n")
cat("  - diff.bindao_targetGenes.txt（目标基因结果）\n")
cat("  - bindao_bindao_volcano.pdf（火山图）\n")
cat("  - bindao_bindao_heatmap_target.pdf（目标基因热图）\n")
cat("  - bindao_bindao_heatmap_diff.pdf（差异基因热图）\n")
cat("\n请运行下一步：step04_bindao_bindao_ML.R\n")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
