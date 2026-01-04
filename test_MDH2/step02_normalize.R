#===============================================================================
# MDH2心衰研究 - 测试版本
# 步骤2：数据标准化
#
# 参考原代码：05.normalize/geo05.normalize.R
#===============================================================================

#加载包
library(limma)

#======================== 参数设置（可修改）========================
geoID <- "GSE26887"           #GEO数据集ID

#设置工作目录（修改为你的路径）
setwd("D:/test_MDH2")         #Windows用户修改这里
#setwd("~/test_MDH2")         #Mac/Linux用户用这个

#======================== 加载数据 ========================
cat("正在加载数据...\n")
load(paste0(geoID, ".RData"))
cat("基因数：", nrow(geneMatrix), "\n")
cat("样本数：", ncol(geneMatrix), "\n")

#======================== 查看分组信息 ========================
cat("\n========== 样本分组信息 ==========\n")
print(phenoData[, c("title", "source_name_ch1")])

#======================== 设置分组（需要根据实际数据修改）========================
#GSE26887的分组信息在source_name_ch1列
#  - "Healthy control" = Control
#  - 其他 = Disease（心衰）

#根据phenoData设置分组
group <- ifelse(grepl("control|healthy|normal", phenoData$source_name_ch1, ignore.case = TRUE),
                "Control", "Treat")
names(group) <- rownames(phenoData)

cat("\n分组结果：\n")
print(table(group))

#获取Control和Treat的样本名
sampleName_Control <- names(group)[group == "Control"]
sampleName_Treat <- names(group)[group == "Treat"]

cat("\nControl样本数：", length(sampleName_Control), "\n")
cat("Treat样本数：", length(sampleName_Treat), "\n")

#======================== 数据标准化 ========================
cat("\n正在进行数据标准化...\n")

#检查是否需要log2转换
data <- geneMatrix
qx <- as.numeric(quantile(data, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0))

if (LogC) {
  data[data < 0] <- 0
  data <- log2(data + 1)
  cat("已进行log2转换\n")
} else {
  cat("数据已是log2格式，无需转换\n")
}

#分位数标准化
data <- normalizeBetweenArrays(data)
cat("已进行分位数标准化\n")

#======================== 整理输出数据 ========================
#按分组排序
conData <- data[, sampleName_Control, drop = FALSE]
treatData <- data[, sampleName_Treat, drop = FALSE]
data <- cbind(conData, treatData)

#添加分组标签到列名
Type <- c(rep("Control", length(sampleName_Control)),
          rep("Treat", length(sampleName_Treat)))
colnames(data) <- paste0(colnames(data), "_", Type)

#======================== 保存标准化数据 ========================
outData <- rbind(id = colnames(data), data)
write.table(outData, file = paste0(geoID, ".normalize.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE)

#保存R对象
normalizedMatrix <- data
save(normalizedMatrix, group, file = paste0(geoID, ".normalized.RData"))

cat("\n========== 标准化完成 ==========\n")
cat("保存文件：\n")
cat("  -", paste0(geoID, ".normalize.txt"), "\n")
cat("  -", paste0(geoID, ".normalized.RData"), "\n")
cat("\n请运行下一步：step03_bindao_diff.R\n")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
