#===============================================================================
# 步骤1：下载GEO数据
# 数据集：GSE59867（人类心梗后心衰，111例，Bulk数据）
#===============================================================================

library(GEOquery)
library(limma)

#======================== 参数设置（可修改）========================
geoID <- "GSE59867"     # GEO数据集ID
workDir <- "./01_data"  # 数据保存目录

#======================== 创建目录 ========================
if (!dir.exists(workDir)) {
  dir.create(workDir, recursive = TRUE)
}
setwd(workDir)

#======================== 下载数据 ========================
cat("正在下载", geoID, "数据集...\n")
cat("这可能需要几分钟，请耐心等待...\n\n")

# 下载GEO数据
gset <- getGEO(geoID, GSEMatrix = TRUE, getGPL = TRUE)
gset <- gset[[1]]

# 提取表达矩阵
expMatrix <- exprs(gset)
cat("表达矩阵维度：", nrow(expMatrix), "个探针 x", ncol(expMatrix), "个样本\n")

# 提取样本信息
phenoData <- pData(gset)
cat("样本信息列：\n")
print(colnames(phenoData))

#======================== 查看分组信息 ========================
cat("\n\n========== 样本分组信息 ==========\n")
# 常见的分组列名
groupCols <- c("characteristics_ch1", "source_name_ch1", "title",
               "characteristics_ch1.1", "disease state:ch1")
for (col in groupCols) {
  if (col %in% colnames(phenoData)) {
    cat(paste0("\n", col, ":\n"))
    print(table(phenoData[[col]]))
  }
}

#======================== 探针注释 ========================
cat("\n正在进行探针注释...\n")
# 获取注释信息
gpl <- getGEO(annotation(gset))
annot <- Table(gpl)

# 查看注释列
cat("注释信息列：\n")
print(head(colnames(annot), 20))

# 通常基因符号在这些列中
geneCols <- c("Gene Symbol", "GENE_SYMBOL", "gene_assignment", "Symbol")
geneCol <- NULL
for (col in geneCols) {
  if (col %in% colnames(annot)) {
    geneCol <- col
    break
  }
}

if (!is.null(geneCol)) {
  cat(paste0("使用 ", geneCol, " 列进行基因注释\n"))
}

#======================== 保存原始数据 ========================
# 保存表达矩阵
write.table(expMatrix, file = paste0(geoID, "_expression_raw.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# 保存样本信息
write.table(phenoData, file = paste0(geoID, "_phenotype.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# 保存注释信息
write.table(annot, file = paste0(geoID, "_annotation.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 保存R对象
save(gset, expMatrix, phenoData, annot, file = paste0(geoID, "_rawdata.RData"))

cat("\n========== 下载完成 ==========\n")
cat("保存文件：\n")
cat(paste0("  - ", geoID, "_expression_raw.txt（表达矩阵）\n"))
cat(paste0("  - ", geoID, "_phenotype.txt（样本信息）\n"))
cat(paste0("  - ", geoID, "_annotation.txt（探针注释）\n"))
cat(paste0("  - ", geoID, "_rawdata.RData（R对象）\n"))
cat("\n请查看样本分组信息，确定对照组和疾病组！\n")
