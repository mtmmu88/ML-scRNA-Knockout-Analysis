#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤2：数据标准化（每个数据集单独标准化）
#
# 参考原代码：05.normalize/geo05.normalize.R
#===============================================================================

#加载包
library(limma)

#======================== 参数设置 ========================
datasets <- c("GSE59867", "GSE57338", "GSE66360")

#设置工作目录
setwd("D:/MDH2_HF_full")         #Windows用户
#setwd("~/MDH2_HF_full")         #Mac/Linux用户

#======================== 标准化函数 ========================
normalize_dataset <- function(geoID) {
  cat("\n========================================\n")
  cat("正在标准化：", geoID, "\n")
  cat("========================================\n")

  #加载数据
  load(paste0(geoID, ".RData"))
  cat("基因数：", nrow(geneMatrix), "\n")
  cat("样本数：", ncol(geneMatrix), "\n")

  #======================== 查看分组信息 ========================
  cat("\n样本信息列：\n")
  print(colnames(phenoData)[1:min(10, ncol(phenoData))])

  #查找分组信息（不同数据集列名不同）
  groupCol <- NULL
  possibleCols <- c("source_name_ch1", "characteristics_ch1",
                    "disease state:ch1", "group:ch1", "condition:ch1")
  for (col in possibleCols) {
    if (col %in% colnames(phenoData)) {
      groupCol <- col
      break
    }
  }

  if (!is.null(groupCol)) {
    cat("\n分组列：", groupCol, "\n")
    print(table(phenoData[[groupCol]]))
  }

  #======================== 设置分组 ========================
  #根据不同数据集设置分组规则
  if (geoID == "GSE59867") {
    #GSE59867: STEMI患者 vs 对照
    group <- ifelse(grepl("STEMI|MI|infarction", phenoData$source_name_ch1, ignore.case = TRUE),
                    "Treat", "Control")
  } else if (geoID == "GSE57338") {
    #GSE57338: 心衰 vs 正常
    group <- ifelse(grepl("heart failure|HF|failing", phenoData$source_name_ch1, ignore.case = TRUE),
                    "Treat", "Control")
  } else if (geoID == "GSE66360") {
    #GSE66360: 急性心梗 vs 对照
    group <- ifelse(grepl("AMI|acute|infarction", phenoData$source_name_ch1, ignore.case = TRUE),
                    "Treat", "Control")
  } else {
    #默认规则
    group <- ifelse(grepl("control|normal|healthy", phenoData$source_name_ch1, ignore.case = TRUE),
                    "Control", "Treat")
  }

  names(group) <- rownames(phenoData)
  cat("\n分组结果：\n")
  print(table(group))

  #======================== 数据标准化 ========================
  data <- geneMatrix

  #检查是否需要log2转换
  qx <- as.numeric(quantile(data, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0))

  if (LogC) {
    data[data < 0] <- 0
    data <- log2(data + 1)
    cat("已进行log2转换\n")
  }

  #分位数标准化
  data <- normalizeBetweenArrays(data)
  cat("已进行分位数标准化\n")

  #======================== 整理输出 ========================
  #确保样本顺序与分组一致
  data <- data[, names(group)]

  #添加分组标签到列名
  colnames(data) <- paste0(geoID, "_", colnames(data), "_", group)

  #======================== 保存 ========================
  outData <- rbind(id = colnames(data), data)
  write.table(outData, file = paste0(geoID, ".normalize.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE)

  save(data, group, file = paste0(geoID, ".normalized.RData"))

  cat("已保存：", geoID, ".normalize.txt\n")

  return(list(data = data, group = group))
}

#======================== 批量标准化 ========================
cat("开始标准化所有数据集...\n")

normalized_data <- list()
for (geoID in datasets) {
  tryCatch({
    normalized_data[[geoID]] <- normalize_dataset(geoID)
  }, error = function(e) {
    cat("标准化失败：", geoID, "\n")
    cat("错误：", conditionMessage(e), "\n")
    cat("请检查phenoData的分组列，可能需要手动调整分组代码\n")
  })
}

#======================== 汇总 ========================
cat("\n\n========== 标准化完成 ==========\n")
for (geoID in names(normalized_data)) {
  if (!is.null(normalized_data[[geoID]])) {
    cat("[OK] ", geoID, "\n")
    print(table(normalized_data[[geoID]]$group))
  }
}

cat("\n请运行下一步：step03_bindao_bindao_bindao_sva.R（批次校正）\n")


######生信自学网: https://www.biowolf.cn/
######MDH2心衰研究完整版
