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
library(AnnotationDbi)

#======================== 网络设置 ========================
options(timeout = 600)
maxRetry <- 3

#======================== 参数设置 ========================
datasets <- c("GSE59867", "GSE57338", "GSE66360")

#设置工作目录（修改为你的路径）
setwd("D:/MDH2_HF_full")         #Windows用户
#setwd("~/MDH2_HF_full")         #Mac/Linux用户

#目标基因
targetGenes <- c("MDH2", "SIRT5", "GPX4", "ACSL4", "TFRC", "SLC7A11",
                 "LPCAT3", "NCOA4", "FTH1", "FTL", "NFE2L2", "VDAC2")

#平台对应的注释包
platformDB <- list(
  "GPL6244" = "hugene10sttranscriptcluster.db",
  "GPL11532" = "hugene11sttranscriptcluster.db",
  "GPL570" = "hgu133plus2.db"
)

#======================== 下载函数（带重试） ========================
download_with_retry <- function(geoID, getGPL = FALSE, maxRetry = 3) {
  for (attempt in 1:maxRetry) {
    cat("  尝试第", attempt, "次下载...\n")
    result <- tryCatch({
      gset <- getGEO(geoID, GSEMatrix = TRUE, getGPL = getGPL, destdir = ".")
      return(gset)
    }, error = function(e) {
      cat("  下载失败：", conditionMessage(e), "\n")
      if (attempt < maxRetry) {
        waitTime <- 2^attempt * 5
        cat("  等待", waitTime, "秒后重试...\n")
        Sys.sleep(waitTime)
      }
      return(NULL)
    })
    if (!is.null(result)) return(result)
  }
  return(NULL)
}

#======================== 处理函数 ========================
download_and_process <- function(geoID) {
  cat("\n========================================\n")
  cat("正在处理：", geoID, "\n")
  cat("========================================\n")

  #下载数据（不下载GPL，使用本地注释包）
  cat("下载中（超时设置：10分钟）...\n")

  #先尝试不下载GPL
  gset <- download_with_retry(geoID, getGPL = FALSE, maxRetry)

  if (is.null(gset)) {
    cat("错误：", geoID, "下载失败\n")
    return(NULL)
  }

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

  #获取平台
  gpl <- annotation(gset)
  cat("平台：", gpl, "\n")

  #======================== 探针转基因（使用注释包） ========================
  probe2gene <- NULL

  #检查是否有对应的注释包
  if (gpl %in% names(platformDB)) {
    dbName <- platformDB[[gpl]]
    cat("使用注释包：", dbName, "\n")

    #尝试加载注释包
    if (require(dbName, character.only = TRUE, quietly = TRUE)) {
      db <- get(dbName)

      #获取探针到基因符号的映射
      mapping <- tryCatch({
        select(db, keys = rownames(expMatrix), columns = "SYMBOL", keytype = "PROBEID")
      }, error = function(e) {
        cat("注释包查询失败，尝试其他方法...\n")
        NULL
      })

      if (!is.null(mapping)) {
        probe2gene <- mapping
        colnames(probe2gene) <- c("probe_id", "gene_symbol")
      }
    } else {
      cat("注释包未安装，尝试下载GPL文件...\n")
    }
  }

  #如果注释包不可用，尝试下载GPL或使用本地缓存
  if (is.null(probe2gene)) {
    cat("尝试从GPL文件获取注释...\n")

    platInfo <- tryCatch({
      getGEO(gpl, destdir = ".")
    }, error = function(e) {
      cat("GPL下载失败：", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(platInfo)) {
      annot <- Table(platInfo)

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

      probe2gene <- annot[, c("ID", geneCol)]
      colnames(probe2gene) <- c("probe_id", "gene_symbol")

      #清理基因符号
      if (geneCol == "gene_assignment") {
        probe2gene$gene_symbol <- gsub(" ///.*", "", probe2gene$gene_symbol)
        probe2gene$gene_symbol <- sapply(strsplit(probe2gene$gene_symbol, " // "), function(x) {
          if (length(x) >= 2) return(x[2]) else return("")
        })
      } else {
        probe2gene$gene_symbol <- gsub(" ///.*", "", probe2gene$gene_symbol)
      }
    } else {
      cat("错误：无法获取探针注释\n")
      return(NULL)
    }
  }

  #清理
  probe2gene$gene_symbol <- trimws(probe2gene$gene_symbol)
  probe2gene <- probe2gene[probe2gene$gene_symbol != "", ]
  probe2gene <- probe2gene[probe2gene$gene_symbol != "---", ]
  probe2gene <- probe2gene[!is.na(probe2gene$gene_symbol), ]
  probe2gene <- unique(probe2gene)

  cat("有效探针-基因映射：", nrow(probe2gene), "\n")

  #======================== 构建基因表达矩阵 ========================
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
  outMatrix <- rbind(id = colnames(geneMatrix), geneMatrix)
  write.table(outMatrix, file = paste0(geoID, ".geneMatrix.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE)

  write.table(phenoData, file = paste0(geoID, ".phenotype.txt"),
              sep = "\t", quote = FALSE, row.names = TRUE)

  save(geneMatrix, phenoData, file = paste0(geoID, ".RData"))

  cat("已保存：", geoID, ".geneMatrix.txt\n")

  return(list(geneMatrix = geneMatrix, phenoData = phenoData))
}

#======================== 批量下载 ========================
cat("开始下载3个数据集...\n")
cat("这可能需要10-20分钟，请耐心等待...\n")
cat("\n提示：如果注释包未安装，请先运行：\n")
cat("BiocManager::install(c('hugene10sttranscriptcluster.db',\n")
cat("                       'hugene11sttranscriptcluster.db',\n")
cat("                       'hgu133plus2.db'))\n\n")

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
