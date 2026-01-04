#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤6：单细胞虚拟敲除分析（scTenifoldKnk）
#
# 参考原代码：26.scTenifoldKnk/toxicityKnk26.scTenifoldKnk.R
# 数据集：GSE135310（小鼠心梗后单细胞）
#===============================================================================

#加载包
library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
set.seed(123)

#======================== 参数设置 ========================
target_gene <- "Mdh2"         #小鼠基因用首字母大写
#target_gene <- "Sirt5"       #也可以敲除Sirt5

#设置工作目录
setwd("D:/MDH2_HF_full")       #Windows用户
#setwd("~/MDH2_HF_full")       #Mac/Linux用户

#======================== 数据准备说明 ========================
cat("========== 单细胞数据准备说明 ==========\n")
cat("\n你需要：\n")
cat("1. 从GEO下载 GSE135310（小鼠心梗后单细胞）\n")
cat("   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135310\n")
cat("\n2. 用Seurat处理成Seurat对象，保存为 seurat.Rdata\n")
cat("   对象名称需要是 'pbmc'\n")
cat("\n3. 确保对象包含 'Type' 列（Control/Disease）\n")

#======================== 检查数据文件 ========================
if (!file.exists("seurat.Rdata")) {
  cat("\n[警告] 未找到 seurat.Rdata 文件\n")
  cat("请先准备单细胞数据，或跳过此步骤\n")
  cat("\n如果你想跳过单细胞分析，直接运行：step07_bindao_bindao_enrichment.R\n")
  stop("需要准备单细胞数据")
}

#======================== 加载Seurat对象 ========================
cat("\n正在加载Seurat对象...\n")
load("seurat.Rdata")

#检查目标基因
if (!target_gene %in% rownames(pbmc)) {
  #尝试不同大小写
  possible_names <- c(target_gene, toupper(target_gene), tolower(target_gene),
                      paste0(toupper(substr(target_gene, 1, 1)),
                             tolower(substr(target_gene, 2, nchar(target_gene)))))
  for (name in possible_names) {
    if (name %in% rownames(pbmc)) {
      target_gene <- name
      break
    }
  }
}

if (!target_gene %in% rownames(pbmc)) {
  stop(paste0("目标基因 ", target_gene, " 不在数据中"))
}
cat("目标基因：", target_gene, "\n")

#只保留疾病组样本
if ("Type" %in% colnames(pbmc@meta.data)) {
  pbmc <- subset(pbmc, subset = Type == "Disease")
}

cat("细胞数：", ncol(pbmc), "\n")
cat("基因数：", nrow(pbmc), "\n")

#======================== 提取表达矩阵 ========================
countMat <- GetAssayData(pbmc, layer = "counts")

#提取高可变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
hvgs <- VariableFeatures(pbmc)

#确保目标基因包含在内
data <- as.data.frame(countMat[unique(c(target_gene, hvgs)), ])
cat("分析基因数：", nrow(data), "\n")

#清理内存
rm(pbmc, countMat)
gc()

#======================== 执行虚拟敲除 ========================
cat("\n========== 开始scTenifoldKnk分析 ==========\n")
cat("这可能需要30分钟到数小时...\n\n")

result <- scTenifoldKnk(
  countMatrix = data,
  gKO = target_gene,
  qc_mtThreshold = 0.1,
  qc_minLSize = 1000,
  nc_nNet = 10,
  nc_nCells = 500
)

#======================== 提取结果 ========================
df <- result$diffRegulation
df <- df[df$gene != target_gene, ]

#显著差异基因
sigDiff <- df[df$p.adj < 0.05, ]
cat("\n敲除", target_gene, "后显著变化的基因：", nrow(sigDiff), "个\n")

#保存结果
write.table(df, file = paste0("knockout_", target_gene, "_all.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigDiff, file = paste0("knockout_", target_gene, "_sig.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#======================== 检查铁死亡基因 ========================
cat("\n========== 铁死亡相关基因变化 ==========\n")

#小鼠铁死亡基因（首字母大写格式）
ferroptosisGenes_mouse <- c(
  "Acsl4", "Lpcat3", "Alox15", "Tfrc", "Ncoa4",
  "Slc7a11", "Gpx4", "Fth1", "Ftl1", "Nfe2l2",
  "Vdac2", "Vdac3", "Cs", "Aco2", "Hmox1"
)

ferroptosis_result <- df[df$gene %in% ferroptosisGenes_mouse, ]
if (nrow(ferroptosis_result) > 0) {
  ferroptosis_result <- ferroptosis_result[order(ferroptosis_result$p.adj), ]
  print(ferroptosis_result)
  write.table(ferroptosis_result,
              file = paste0("knockout_", target_gene, "_ferroptosis.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

#======================== 柱状图 ========================
cat("\n正在绑制图形...\n")

top_genes <- head(df[order(-df$FC), ], 20)
p1 <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC)) +
  geom_bar(stat = 'identity', fill = '#5A9BD4') +
  coord_flip() +
  labs(title = paste0("Top 20 Genes After ", target_gene, " Knockout"),
       x = "Gene", y = "Fold Change") +
  theme_minimal()

pdf(paste0("knockout_", target_gene, "_barplot.pdf"), width = 6, height = 5)
print(p1)
dev.off()

#======================== 火山图 ========================
df$log_p.adj <- -log10(df$p.adj)
df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")
label_genes <- subset(df, p.adj < 0.05 & abs(Z) > 2)

y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)
p2 <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = 30) +
  labs(title = paste0(target_gene, " Knockout - Volcano Plot"),
       x = "Z-score", y = "-log10(p.adj)") +
  theme_classic() +
  coord_cartesian(ylim = c(0, y_upper)) +
  theme(legend.position = "none")

pdf(paste0("knockout_", target_gene, "_volcano.pdf"), width = 8, height = 6)
print(p2)
dev.off()

#======================== 保存结果 ========================
save(result, df, sigDiff, file = paste0("knockout_", target_gene, ".RData"))

cat("\n========== 敲除分析完成 ==========\n")
cat("输出文件：\n")
cat(paste0("  - knockout_", target_gene, "_all.txt\n"))
cat(paste0("  - knockout_", target_gene, "_sig.txt\n"))
cat(paste0("  - knockout_", target_gene, "_ferroptosis.txt\n"))
cat(paste0("  - knockout_", target_gene, "_barplot.pdf\n"))
cat(paste0("  - knockout_", target_gene, "_volcano.pdf\n"))
cat("\n请运行下一步：step07_bindao_bindao_enrichment.R\n")


######生信自学网: https://www.biowolf.cn/
######MDH2心衰研究完整版
