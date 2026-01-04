#===============================================================================
# MDH2心衰研究 - 测试版本
# 步骤1：下载GEO数据
#
# 数据集：GSE26887（心衰，24例，用于测试）
# 正式分析时可换成：GSE59867（心梗后心衰，111例）
#===============================================================================

#加载包
library(GEOquery)
library(limma)

#======================== 参数设置（可修改）========================
geoID <- "GSE26887"           #GEO数据集ID（测试用小数据集）
#geoID <- "GSE59867"          #正式分析用这个（心梗后心衰，111例）

#设置工作目录（修改为你的路径）
setwd("D:/test_MDH2")         #Windows用户修改这里
#setwd("~/test_MDH2")         #Mac/Linux用户用这个

#======================== 下载数据 ========================
cat("正在下载", geoID, "数据集...\n")
cat("这可能需要几分钟，请耐心等待...\n\n")

#下载GEO数据
gset <- getGEO(geoID, GSEMatrix = TRUE, getGPL = TRUE, destdir = ".")

#如果返回的是列表，取第一个元素
if (length(gset) > 1) {
  idx <- grep(geoID, attr(gset, "names"))
  gset <- gset[[idx]]
} else {
  gset <- gset[[1]]
}

#======================== 提取表达矩阵 ========================
#提取表达矩阵
expMatrix <- exprs(gset)
cat("表达矩阵维度：", nrow(expMatrix), "个探针 x", ncol(expMatrix), "个样本\n")

#提取样本信息
phenoData <- pData(gset)
cat("\n样本信息列：\n")
print(colnames(phenoData))

#======================== 查看分组信息 ========================
cat("\n\n========== 样本分组信息 ==========\n")
#查看常见的分组列
if ("source_name_ch1" %in% colnames(phenoData)) {
  cat("\nsource_name_ch1:\n")
  print(table(phenoData$source_name_ch1))
}
if ("characteristics_ch1" %in% colnames(phenoData)) {
  cat("\ncharacteristics_ch1:\n")
  print(table(phenoData$characteristics_ch1))
}

#======================== 探针注释 ========================
cat("\n正在获取探针注释...\n")
#获取平台注释
gpl <- annotation(gset)
cat("平台ID：", gpl, "\n")

#从GEO获取注释表
platInfo <- getGEO(gpl, destdir = ".")
annot <- Table(platInfo)
cat("注释表维度：", nrow(annot), "行 x", ncol(annot), "列\n")
cat("注释列名：\n")
print(colnames(annot))

#======================== 探针转基因 ========================
#根据平台选择基因符号列（不同平台列名不同）
if ("Gene Symbol" %in% colnames(annot)) {
  geneCol <- "Gene Symbol"
} else if ("gene_assignment" %in% colnames(annot)) {
  geneCol <- "gene_assignment"
} else if ("GENE_SYMBOL" %in% colnames(annot)) {
  geneCol <- "GENE_SYMBOL"
} else {
  cat("警告：未找到基因符号列，请手动检查annot对象\n")
  geneCol <- NULL
}

if (!is.null(geneCol)) {
  cat("使用", geneCol, "列进行基因注释\n")

  #创建探针到基因的映射
  probe2gene <- annot[, c("ID", geneCol)]
  colnames(probe2gene) <- c("probe_id", "gene_symbol")

  #清理基因符号
  probe2gene$gene_symbol <- gsub(" ///.*", "", probe2gene$gene_symbol)  #去除多基因
  probe2gene <- probe2gene[probe2gene$gene_symbol != "", ]              #去除空值
  probe2gene <- probe2gene[!is.na(probe2gene$gene_symbol), ]            #去除NA

  cat("有效探针-基因映射：", nrow(probe2gene), "个\n")
}

#======================== 构建基因表达矩阵 ========================
cat("\n正在构建基因表达矩阵...\n")

#合并探针和基因符号
expDF <- as.data.frame(expMatrix)
expDF$probe_id <- rownames(expDF)
expDF <- merge(probe2gene, expDF, by = "probe_id")
expDF <- expDF[, -1]  #删除probe_id列

#同基因多探针取平均
geneMatrix <- aggregate(. ~ gene_symbol, data = expDF, FUN = mean)
rownames(geneMatrix) <- geneMatrix$gene_symbol
geneMatrix <- geneMatrix[, -1]
geneMatrix <- as.matrix(geneMatrix)

cat("基因表达矩阵：", nrow(geneMatrix), "个基因 x", ncol(geneMatrix), "个样本\n")

#======================== 检查目标基因 ========================
cat("\n========== 检查目标基因 ==========\n")
targetGenes <- c("MDH2", "SIRT5", "GPX4", "ACSL4", "TFRC", "SLC7A11")
for (gene in targetGenes) {
  if (gene %in% rownames(geneMatrix)) {
    cat("[OK]", gene, "存在\n")
  } else {
    cat("[MISS]", gene, "不存在\n")
  }
}

#======================== 保存数据 ========================
cat("\n正在保存数据...\n")

#保存基因表达矩阵
outMatrix <- rbind(id = colnames(geneMatrix), geneMatrix)
write.table(outMatrix, file = paste0(geoID, ".geneMatrix.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE)

#保存样本信息
write.table(phenoData, file = paste0(geoID, ".phenotype.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

#保存R对象（方便后续使用）
save(geneMatrix, phenoData, probe2gene, file = paste0(geoID, ".RData"))

cat("\n========== 下载完成 ==========\n")
cat("保存文件：\n")
cat("  -", paste0(geoID, ".geneMatrix.txt"), "（基因表达矩阵）\n")
cat("  -", paste0(geoID, ".phenotype.txt"), "（样本信息）\n")
cat("  -", paste0(geoID, ".RData"), "（R数据对象）\n")
cat("\n请查看样本分组信息，确定Control和Disease分组！\n")
cat("然后运行下一步：step02_bindao_normalize.R\n")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
