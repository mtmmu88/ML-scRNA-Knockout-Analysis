#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤1：下载GEO数据（3个Bulk数据集）
#===============================================================================

library(GEOquery)
library(limma)

options(timeout = 600)

#设置工作目录
setwd("D:/MDH2_HF_full")

#目标基因
targetGenes <- c("MDH2", "SIRT5", "GPX4", "ACSL4", "TFRC", "SLC7A11",
                 "LPCAT3", "NCOA4", "FTH1", "FTL", "NFE2L2", "VDAC2")

#======================== 获取探针映射函数 ========================
get_probe2gene <- function(gpl, probeIDs) {

  if (gpl == "GPL6244") {
    library(hugene10sttranscriptcluster.db)
    mapping <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                      keys = probeIDs,
                                      columns = "SYMBOL",
                                      keytype = "PROBEID")
  } else if (gpl == "GPL11532") {
    library(hugene11sttranscriptcluster.db)
    mapping <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
                                      keys = probeIDs,
                                      columns = "SYMBOL",
                                      keytype = "PROBEID")
  } else if (gpl == "GPL570") {
    library(hgu133plus2.db)
    mapping <- AnnotationDbi::select(hgu133plus2.db,
                                      keys = probeIDs,
                                      columns = "SYMBOL",
                                      keytype = "PROBEID")
  } else {
    stop("未知平台：", gpl)
  }

  colnames(mapping) <- c("probe_id", "gene_symbol")
  mapping <- mapping[!is.na(mapping$gene_symbol), ]
  mapping <- mapping[mapping$gene_symbol != "", ]
  return(mapping)
}

#======================== 处理函数 ========================
process_dataset <- function(geoID) {
  cat("\n========================================\n")
  cat("正在处理：", geoID, "\n")
  cat("========================================\n")

  #下载数据（不下载GPL）
  cat("下载中...\n")
  gset <- getGEO(geoID, GSEMatrix = TRUE, getGPL = FALSE, destdir = ".")

  if (length(gset) > 1) {
    gset <- gset[[1]]
  } else {
    gset <- gset[[1]]
  }

  #提取表达矩阵
  expMatrix <- exprs(gset)
  cat("探针数：", nrow(expMatrix), "\n")
  cat("样本数：", ncol(expMatrix), "\n")

  #提取样本信息
  phenoData <- pData(gset)

  #获取平台
  gpl <- annotation(gset)
  cat("平台：", gpl, "\n")

  #获取探针-基因映射
  cat("获取基因注释...\n")
  probe2gene <- get_probe2gene(gpl, rownames(expMatrix))
  cat("有效探针-基因映射：", nrow(probe2gene), "\n")

  #构建基因表达矩阵
  expDF <- as.data.frame(expMatrix)
  expDF$probe_id <- rownames(expDF)
  expDF <- merge(probe2gene, expDF, by = "probe_id")
  expDF <- expDF[, -1]

  #同基因多探针取平均
  geneMatrix <- aggregate(. ~ gene_symbol, data = expDF, FUN = mean)
  rownames(geneMatrix) <- geneMatrix$gene_symbol
  geneMatrix <- geneMatrix[, -1]
  geneMatrix <- as.matrix(geneMatrix)

  cat("基因表达矩阵：", nrow(geneMatrix), "基因 x", ncol(geneMatrix), "样本\n")

  #检查目标基因
  cat("\n目标基因检查：\n")
  found <- sum(targetGenes %in% rownames(geneMatrix))
  cat("  找到", found, "/", length(targetGenes), "个目标基因\n")
  cat("  MDH2:", ifelse("MDH2" %in% rownames(geneMatrix), "[OK]", "[MISS]"), "\n")
  cat("  SIRT5:", ifelse("SIRT5" %in% rownames(geneMatrix), "[OK]", "[MISS]"), "\n")

  #保存数据
  outMatrix <- rbind(id = colnames(geneMatrix), geneMatrix)
  write.table(outMatrix, file = paste0(geoID, ".geneMatrix.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE)
  write.table(phenoData, file = paste0(geoID, ".phenotype.txt"),
              sep = "\t", quote = FALSE, row.names = TRUE)
  save(geneMatrix, phenoData, file = paste0(geoID, ".RData"))

  cat("已保存：", geoID, ".geneMatrix.txt\n")
  return(list(geneMatrix = geneMatrix, phenoData = phenoData))
}

#======================== 批量处理 ========================
datasets <- c("GSE59867", "GSE57338", "GSE66360")

cat("提示：请先安装注释包：\n")
cat("BiocManager::install(c('hugene10sttranscriptcluster.db',\n")
cat("                       'hugene11sttranscriptcluster.db',\n")
cat("                       'hgu133plus2.db'))\n\n")

results <- list()
for (geoID in datasets) {
  tryCatch({
    results[[geoID]] <- process_dataset(geoID)
  }, error = function(e) {
    cat("处理失败：", geoID, "\n")
    cat("错误：", conditionMessage(e), "\n")
  })
}

#汇总
cat("\n\n========== 完成 ==========\n")
for (geoID in names(results)) {
  if (!is.null(results[[geoID]])) {
    cat("[OK]", geoID, ":", nrow(results[[geoID]]$geneMatrix), "基因,",
        ncol(results[[geoID]]$geneMatrix), "样本\n")
  }
}

cat("\n请运行下一步：step02_normalize.R\n")
