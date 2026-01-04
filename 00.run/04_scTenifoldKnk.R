#===============================================================================
# 步骤4：scTenifoldKnk 虚拟基因敲除分析
# 模拟敲除MDH2后对其他基因的影响
# 需要单细胞数据！（如GSE132146）
#===============================================================================

library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
set.seed(123)

#======================== 参数设置（需要修改）========================
target_gene <- "MDH2"       # 要敲除的基因（可以改成SIRT5测试）
geoID <- "GSE132146"        # 单细胞数据集ID
workDir <- "./04_knockout"
dataDir <- "./01_data"

#======================== 创建目录 ========================
if (!dir.exists(workDir)) {
  dir.create(workDir, recursive = TRUE)
}
setwd(workDir)

#======================== 下载单细胞数据（如果没有）========================
# 注意：单细胞数据通常很大，建议手动下载
# 这里提供两种方式

# 方式1：从GEO下载（可能很慢）
if (!file.exists("seurat.Rdata")) {
  cat("正在下载单细胞数据，这可能需要较长时间...\n")
  cat("建议从GEO网站手动下载：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geoID, "\n")

  # 尝试自动下载
  library(GEOquery)
  tryCatch({
    gse <- getGEO(geoID, GSEMatrix = FALSE)
    # 获取supplementary files
    getGEOSuppFiles(geoID, baseDir = dataDir)
    cat("下载完成，请解压文件并处理成Seurat对象\n")
  }, error = function(e) {
    cat("自动下载失败，请手动下载数据\n")
  })

  stop("请先准备好单细胞数据的Seurat对象（seurat.Rdata）")
}

#======================== 加载Seurat对象 ========================
cat("正在加载Seurat对象...\n")
load("seurat.Rdata")  # 需要包含名为 "pbmc" 的Seurat对象

# 检查目标基因是否存在
if (!target_gene %in% rownames(pbmc)) {
  # 尝试小写/大写转换（小鼠基因）
  if (tolower(target_gene) %in% tolower(rownames(pbmc))) {
    target_gene_new <- rownames(pbmc)[tolower(rownames(pbmc)) == tolower(target_gene)]
    cat("找到对应基因：", target_gene_new, "\n")
    target_gene <- target_gene_new
  } else {
    stop(paste0("目标基因 ", target_gene, " 不在数据中！"))
  }
}

# 如果有分组信息，只保留疾病组
if ("Type" %in% colnames(pbmc@meta.data)) {
  pbmc <- subset(pbmc, subset = Type == "Disease")
  cat("已筛选疾病组样本\n")
}

cat("细胞数量：", ncol(pbmc), "\n")
cat("基因数量：", nrow(pbmc), "\n")

#======================== 提取表达矩阵 ========================
countMat <- GetAssayData(pbmc, layer = "counts")

# 提取高可变基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 5000)
hvgs <- VariableFeatures(pbmc)

# 确保目标基因在分析中
data <- as.data.frame(countMat[unique(c(target_gene, hvgs)), ])
cat("用于分析的基因数：", nrow(data), "\n")
cat("用于分析的细胞数：", ncol(data), "\n")

# 清理内存
rm(pbmc, countMat)
gc()

#======================== 执行虚拟敲除 ========================
cat("\n========== 开始scTenifoldKnk分析 ==========\n")
cat("这可能需要30分钟到数小时，取决于数据大小...\n\n")

result <- scTenifoldKnk(
  countMatrix = data,
  gKO = target_gene,           # 需要敲除的基因
  qc_mtThreshold = 0.1,        # 线粒体基因阈值
  qc_minLSize = 1000,          # 最小文库大小
  nc_nNet = 10,                # 子网络数量
  nc_nCells = 500              # 每个网络的细胞数
)

#======================== 提取差异调控结果 ========================
df <- result$diffRegulation
df <- df[df$gene != target_gene, ]  # 移除目标基因本身

# 显著差异基因
sigDiff <- df[df$p.adj < 0.05, ]
cat("\n敲除", target_gene, "后显著变化的基因数：", nrow(sigDiff), "\n")

# 保存结果
write.table(df, file = paste0(target_gene, "_knockout_all.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigDiff, file = paste0(target_gene, "_knockout_significant.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#======================== 检查铁死亡基因 ========================
cat("\n========== 铁死亡相关基因变化 ==========\n")
ferroptosisGenes <- c(
  "ACSL4", "LPCAT3", "ALOX15", "TFRC", "SLC7A11", "GPX4",
  "NCOA4", "FTH1", "FTL", "IREB2", "HMOX1", "NFS1",
  "SLC3A2", "GSS", "GCLC", "GCLM", "NFE2L2", "KEAP1",
  "VDAC2", "VDAC3", "CISD1", "CS", "ACO2"
)

# 小鼠基因可能是小写
ferroptosis_mouse <- tolower(ferroptosisGenes)
ferroptosis_check <- c(ferroptosisGenes, ferroptosis_mouse,
                       paste0(substr(ferroptosisGenes, 1, 1),
                              tolower(substr(ferroptosisGenes, 2, nchar(ferroptosisGenes)))))

ferroptosis_result <- df[df$gene %in% ferroptosis_check, ]
if (nrow(ferroptosis_result) > 0) {
  ferroptosis_result <- ferroptosis_result[order(ferroptosis_result$p.adj), ]
  print(ferroptosis_result)
  write.table(ferroptosis_result, file = paste0(target_gene, "_ferroptosis_genes.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("未找到铁死亡相关基因的显著变化\n")
}

#======================== 柱状图（Top20基因）========================
cat("\n正在绑制柱状图...\n")
top_genes <- head(df[order(-df$FC), ], 20)
p1 <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC)) +
  geom_bar(stat = 'identity', fill = '#5A9BD4') +
  coord_flip() +
  labs(title = paste0("Top 20 Genes After ", target_gene, " Knockout"),
       x = "Gene", y = "Fold Change") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(paste0(target_gene, "_barplot.pdf"), width = 6, height = 5)
print(p1)
dev.off()

#======================== 火山图 ========================
cat("正在绘制火山图...\n")
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
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

pdf(paste0(target_gene, "_volcano.pdf"), width = 8, height = 6)
print(p2)
dev.off()

#======================== 保存结果 ========================
save(result, df, sigDiff, file = paste0(target_gene, "_knockout_results.RData"))

cat("\n========== scTenifoldKnk分析完成 ==========\n")
cat("输出文件：\n")
cat(paste0("  - ", target_gene, "_knockout_all.txt（所有基因变化）\n"))
cat(paste0("  - ", target_gene, "_knockout_significant.txt（显著变化基因）\n"))
cat(paste0("  - ", target_gene, "_ferroptosis_genes.txt（铁死亡基因）\n"))
cat(paste0("  - ", target_gene, "_barplot.pdf（柱状图）\n"))
cat(paste0("  - ", target_gene, "_volcano.pdf（火山图）\n"))
cat(paste0("  - ", target_gene, "_knockout_results.RData（R数据）\n"))

cat("\n下一步：用显著变化的基因做富集分析（GO/KEGG）\n")
