#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤4：差异分析
#
# 参考原代码：09.diff/geo09.diff.R
#===============================================================================

#加载包
library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)

#======================== 参数设置 ========================
logFCfilter <- 0.585          #logFC阈值
adj.P.Val.Filter <- 0.05      #校正P值阈值
diseaseName <- "HeartFailure"

#目标基因
targetGenes <- c("MDH2", "SIRT5")

#铁死亡基因
ferroptosisGenes <- c(
  "ACSL4", "LPCAT3", "ALOX15", "TFRC", "NCOA4",
  "SLC7A11", "GPX4", "FTH1", "FTL", "NFE2L2",
  "VDAC2", "VDAC3", "CS", "ACO2", "HMOX1",
  "IREB2", "NFS1", "SLC3A2", "GSS", "GCLC"
)

#设置工作目录
setwd("D:/MDH2_HF_full")         #Windows用户
#setwd("~/MDH2_HF_full")         #Mac/Linux用户

#======================== 加载数据 ========================
cat("正在加载批次校正后的数据...\n")
inputFile <- "merge.normalize.txt"

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
#从列名提取（格式：GSE_Sample_Group）
Type <- gsub(".*_", "", colnames(data))
Project <- gsub("_.*", "", colnames(data))

#按分组排序
data <- data[, order(Type)]
Type <- gsub(".*_", "", colnames(data))
Project <- gsub("_.*", "", colnames(data))

#简化列名（去掉分组后缀，保留GSE_Sample）
colnames(data) <- gsub("_Control$|_Treat$", "", colnames(data))

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
write.table(cbind(gene = rownames(allDiff), allDiff),
            file = "diff.bindao_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#筛选显著差异基因
diffSig <- allDiff[abs(allDiff$logFC) > logFCfilter & allDiff$adj.P.Val < adj.P.Val.Filter, ]
write.table(cbind(gene = rownames(diffSig), diffSig),
            file = "diff.bindao_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n差异分析结果：\n")
cat("  显著差异基因：", nrow(diffSig), "个\n")
cat("  上调基因：", sum(diffSig$logFC > 0), "个\n")
cat("  下调基因：", sum(diffSig$logFC < 0), "个\n")

#======================== MDH2和SIRT5结果 ========================
cat("\n========== MDH2 和 SIRT5 表达情况 ==========\n")
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
  } else {
    cat(paste0("\n", gene, ": 未找到\n"))
  }
}

#======================== 铁死亡基因结果 ========================
cat("\n========== 铁死亡相关基因 ==========\n")
allGenes <- c(targetGenes, ferroptosisGenes)
targetResult <- allDiff[rownames(allDiff) %in% allGenes, ]
targetResult <- targetResult[order(targetResult$P.Value), ]
print(targetResult[, c("logFC", "P.Value", "adj.P.Val")])

write.table(cbind(gene = rownames(targetResult), targetResult),
            file = "diff.bindao_targetGenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#======================== 火山图 ========================
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
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed") +
  geom_hline(yintercept = -log10(adj.P.Val.Filter), linetype = "dashed") +
  geom_text_repel(aes(label = label), size = 4, fontface = "bold", color = "black") +
  labs(title = "Differential Expression (Heart Failure vs Control)",
       x = "log2(Fold Change)", y = "-log10(adj.P.Value)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf("bindao_volcano.pdf", width = 8, height = 6)
print(p)
dev.off()
cat("火山图已保存\n")

#======================== 热图 ========================
cat("正在绘制热图...\n")

#目标基因热图
hmGenes <- allGenes[allGenes %in% rownames(data)]
if (length(hmGenes) > 0) {
  hmExp <- data[hmGenes, ]

  annotation_col <- data.frame(
    Group = ifelse(Type == "Control", "Control", diseaseName),
    Dataset = Project,
    row.names = colnames(data)
  )

  pdf("bindao_heatmap_target.pdf", width = 14, height = 8)
  pheatmap(hmExp,
           annotation_col = annotation_col,
           color = colorRampPalette(c("blue2", "white", "red2"))(50),
           cluster_cols = FALSE,
           show_colnames = FALSE,
           scale = "row",
           main = "MDH2 / SIRT5 / Ferroptosis Genes")
  dev.off()
  cat("目标基因热图已保存\n")
}

#======================== 保存结果 ========================
save(allDiff, diffSig, targetResult, data, Type, Project,
     file = "diff.RData")

cat("\n========== 差异分析完成 ==========\n")
cat("输出文件：\n")
cat("  - diff.bindao_all.txt\n")
cat("  - diff.bindao_sig.txt\n")
cat("  - diff.bindao_targetGenes.txt\n")
cat("  - bindao_volcano.pdf\n")
cat("  - bindao_heatmap_target.pdf\n")
cat("\n请运行下一步：step05_bindao_ML.R\n")


######生信自学网: https://www.biowolf.cn/
######MDH2心衰研究完整版
