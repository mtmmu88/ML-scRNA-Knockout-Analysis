#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤1：下载GEO数据（3个Bulk数据集）
#
# 训练集：GSE59867（心梗后心衰，111例）
# 验证集1：GSE57338（缺血性心衰，231例）
# 验证集2：GSE66360（急性心梗，99例）
#===============================================================================

#加载包
library(GEOquery)
library(limma)

#======================== 参数设置 ========================
#数据集列表
datasets <- c("GSE59867", "GSE57338", "GSE66360")

#设置工作目录（修改为你的路径）
setwd("D:/MDH2_HF_full")         #Windows用户
#setwd("~/MDH2_HF_full")         #Mac/Linux用户

#目标基因
targetGenes <- c("MDH2", "SIRT5", "GPX4", "ACSL4", "TFRC", "SLC7A11",
                 "LPCAT3", "NCOA4", "FTH1", "FTL", "NFE2L2", "VDAC2")

#======================== 下载函数 ========================
download_and_process <- function(geoID) {
  cat("\n========================================\n")
  cat("正在处理：", geoID, "\n")
  cat("========================================\n")

  #下载数据
  cat("下载中...\n")
  gset <- getGEO(geoID, GSEMatrix = TRUE, getGPL = TRUE, destdir = ".")

  if (length(gset) > 1) {
    idx <- grep(geoID, attr(gset, "names"))
    gset <- gset[[idx]]
  } else {
    gset <- gset[[1]]
  }

  #提取表达矩阵
  expMatrix <- exprs(gset)
  cat("探针数：", nrow(expMatrix), "\n")
  cat("样本数：", ncol(expMatrix), "\n")

  #提取样本信息
  phenoData <- pData(gset)

  #获取平台注释
  gpl <- annotation(gset)
  cat("平台：", gpl, "\n")

  platInfo <- getGEO(gpl, destdir = ".")
  annot <- Table(platInfo)

  #======================== 探针转基因 ========================
  #根据平台选择基因符号列
  if ("Gene Symbol" %in% colnames(annot)) {
    geneCol <- "Gene Symbol"
  } else if ("gene_assignment" %in% colnames(annot)) {
    geneCol <- "gene_assignment"
  } else if ("GENE_SYMBOL" %in% colnames(annot)) {
    geneCol <- "GENE_SYMBOL"
  } else if ("Symbol" %in% colnames(annot)) {
    geneCol <- "Symbol"
  } else {
    cat("警告：未找到基因符号列\n")
    print(colnames(annot))
    return(NULL)
  }

  cat("使用列：", geneCol, "\n")

  #创建探针到基因的映射
  probe2gene <- annot[, c("ID", geneCol)]
  colnames(probe2gene) <- c("probe_id", "gene_symbol")

  #清理基因符号（处理不同格式）
  #格式1: "NM_xxx // GENE // description // ..."
  if (grepl("//", probe2gene$gene_symbol[1])) {
    probe2gene$gene_symbol <- gsub(".*// (.*?) //.*", "\\1", probe2gene$gene_symbol)
  }
  #格式2: "GENE /// GENE2"
  probe2gene$gene_symbol <- gsub(" ///.*", "", probe2gene$gene_symbol)
  probe2gene$gene_symbol <- trimws(probe2gene$gene_symbol)

  #去除空值和NA
  probe2gene <- probe2gene[probe2gene$gene_symbol != "", ]
  probe2gene <- probe2gene[probe2gene$gene_symbol != "---", ]
  probe2gene <- probe2gene[!is.na(probe2gene$gene_symbol), ]

  cat("有效探针-基因映射：", nrow(probe2gene), "\n")

  #======================== 构建基因表达矩阵 ========================
  expDF <- as.data.frame(expMatrix)
  expDF$probe_id <- rownames(expDF)
  expDF <- merge(probe2gene, expDF, by = "probe_id")
  expDF <- expDF[, -1]  #删除probe_id列

  #同基因多探针取平均
  geneMatrix <- aggregate(. ~ gene_symbol, data = expDF, FUN = mean)
  rownames(geneMatrix) <- geneMatrix$gene_symbol
  geneMatrix <- geneMatrix[, -1]
  geneMatrix <- as.matrix(geneMatrix)

  cat("基因表达矩阵：", nrow(geneMatrix), "基因 x", ncol(geneMatrix), "样本\n")

  #======================== 检查目标基因 ========================
  cat("\n目标基因检查：\n")
  found <- sum(targetGenes %in% rownames(geneMatrix))
  cat("  找到", found, "/", length(targetGenes), "个目标基因\n")

  if ("MDH2" %in% rownames(geneMatrix)) {
    cat("  [OK] MDH2 存在\n")
  } else {
    cat("  [MISS] MDH2 不存在\n")
  }

  if ("SIRT5" %in% rownames(geneMatrix)) {
    cat("  [OK] SIRT5 存在\n")
  } else {
    cat("  [MISS] SIRT5 不存在\n")
  }

  #======================== 保存数据 ========================
  #保存基因表达矩阵
  outMatrix <- rbind(id = colnames(geneMatrix), geneMatrix)
  write.table(outMatrix, file = paste0(geoID, ".geneMatrix.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE)

  #保存样本信息
  write.table(phenoData, file = paste0(geoID, ".phenotype.txt"),
              sep = "\t", quote = FALSE, row.names = TRUE)

  #保存R对象
  save(geneMatrix, phenoData, file = paste0(geoID, ".RData"))

  cat("已保存：", geoID, ".geneMatrix.txt\n")

  return(list(geneMatrix = geneMatrix, phenoData = phenoData))
}

#======================== 批量下载 ========================
cat("开始下载3个数据集...\n")
cat("这可能需要10-20分钟，请耐心等待...\n")

results <- list()
for (geoID in datasets) {
  tryCatch({
    results[[geoID]] <- download_and_process(geoID)
  }, error = function(e) {
    cat("下载失败：", geoID, "\n")
    cat("错误：", conditionMessage(e), "\n")
  })
}

#======================== 汇总 ========================
cat("\n\n========== 下载完成 ==========\n")
for (geoID in names(results)) {
  if (!is.null(results[[geoID]])) {
    cat("[OK] ", geoID, ":", nrow(results[[geoID]]$geneMatrix), "基因,",
        ncol(results[[geoID]]$geneMatrix), "样本\n")
  }
}

cat("\n请运行下一步：step02_normalize.R\n")


######生信自学网: https://www.biowolf.cn/
######MDH2心衰研究完整版
