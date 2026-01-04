#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤3：批次校正（ComBat）
#
# 参考原代码：06.sva/geo06.sva.R
# 合并多个数据集，去除批次效应
#===============================================================================

#加载包
library(limma)
library(sva)

#======================== 参数设置 ========================
#设置工作目录
setwd("D:/MDH2_HF_full")         #Windows用户
#setwd("~/MDH2_HF_full")         #Mac/Linux用户

#======================== 读取所有normalize文件 ========================
cat("正在读取标准化数据...\n")

files <- dir(pattern = "normalize.txt$")
files <- files[!grepl("merge", files)]  #排除已合并的文件
cat("找到文件：\n")
print(files)

#======================== 获取共同基因 ========================
geneList <- list()

for (file in files) {
  rt <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE)
  geneNames <- as.vector(rt[, 1])
  uniqGene <- unique(geneNames)
  header <- unlist(strsplit(file, "\\.|\\-"))[1]
  geneList[[header]] <- uniqGene
  cat(header, ":", length(uniqGene), "基因\n")
}

#取交集
interGenes <- Reduce(intersect, geneList)
cat("\n共同基因数：", length(interGenes), "\n")

#======================== 数据合并 ========================
cat("\n正在合并数据...\n")

allTab <- data.frame()
batchType <- c()

for (i in 1:length(files)) {
  inputFile <- files[i]
  header <- unlist(strsplit(inputFile, "\\.|\\-"))[1]

  #读取数据
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, 2:ncol(rt)]
  dimnames <- list(rownames(exp), colnames(exp))
  data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
  data <- avereps(data)

  #只保留共同基因
  data <- data[interGenes, ]

  #合并
  if (i == 1) {
    allTab <- data
  } else {
    allTab <- cbind(allTab, data)
  }
  batchType <- c(batchType, rep(i, ncol(data)))

  cat("合并：", header, ",", ncol(data), "样本\n")
}

cat("\n合并后：", nrow(allTab), "基因 x", ncol(allTab), "样本\n")
cat("批次分布：\n")
print(table(batchType))

#======================== 保存合并前的数据 ========================
outTab <- rbind(geneNames = colnames(allTab), allTab)
write.table(outTab, file = "merge.preNorm.txt", sep = "\t", quote = FALSE, col.names = FALSE)
cat("已保存：merge.preNorm.txt（批次校正前）\n")

#======================== ComBat批次校正 ========================
cat("\n正在进行ComBat批次校正...\n")

#执行ComBat
outTab <- ComBat(allTab, batchType, par.prior = TRUE)

#保存校正后的数据
outTab_save <- rbind(geneNames = colnames(outTab), outTab)
write.table(outTab_save, file = "merge.normalize.txt", sep = "\t", quote = FALSE, col.names = FALSE)
cat("已保存：merge.normalize.txt（批次校正后）\n")

#======================== 提取分组信息 ========================
#从列名提取分组（格式：GSE_Sample_Group）
Type <- gsub(".*_", "", colnames(outTab))
Project <- gsub("_.*", "", colnames(outTab))

cat("\n分组统计：\n")
print(table(Type))
cat("\n数据集统计：\n")
print(table(Project))

#保存分组信息
groupInfo <- data.frame(
  Sample = colnames(outTab),
  Project = Project,
  Type = Type
)
write.table(groupInfo, file = "merge.groupInfo.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#保存R对象
mergedData <- outTab
save(mergedData, Type, Project, batchType, file = "merge.RData")

cat("\n========== 批次校正完成 ==========\n")
cat("输出文件：\n")
cat("  - merge.preNorm.txt（校正前）\n")
cat("  - merge.normalize.txt（校正后）\n")
cat("  - merge.groupInfo.txt（分组信息）\n")
cat("  - merge.RData（R对象）\n")
cat("\n请运行下一步：step04_bindao_diff.R\n")


######生信自学网: https://www.biowolf.cn/
######MDH2心衰研究完整版
