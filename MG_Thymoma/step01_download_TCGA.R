#===============================================================================
# 步骤1：下载TCGA-THYM胸腺瘤数据
# 用于胸腺瘤相关重症肌无力(MG)研究
# 数据：116例胸腺瘤，其中34例有MG并发症
#===============================================================================

# ======================== 安装必要的包 ========================
# 首次运行需要安装，之后可以注释掉

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# 安装TCGA相关包
install_if_missing("BiocManager")
install_if_missing("TCGAbiolinks", bioc = TRUE)
install_if_missing("SummarizedExperiment", bioc = TRUE)
install_if_missing("dplyr")
install_if_missing("tidyr")

# ======================== 加载包 ========================
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)

# ======================== 参数设置 ========================
project <- "TCGA-THYM"  # 胸腺瘤项目
workDir <- "./01_TCGA_data"

# ======================== 创建目录 ========================
if (!dir.exists(workDir)) {
  dir.create(workDir, recursive = TRUE)
}
setwd(workDir)

cat("========================================\n")
cat("  TCGA-THYM 胸腺瘤数据下载\n")
cat("  目标：116例样本，34例MG+，82例MG-\n")
cat("========================================\n\n")

# ======================== 方法1：使用TCGAbiolinks下载 ========================
cat("正在查询TCGA-THYM数据...\n")

# 查询RNA-seq表达数据
query_exp <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"  # 使用STAR比对的counts数据
)

# 显示查询结果
cat("\n查询到的样本数量：", nrow(getResults(query_exp)), "\n")

# 下载数据（可能需要较长时间）
cat("\n正在下载表达数据，请稍候...\n")
cat("（首次下载可能需要10-30分钟，取决于网络速度）\n\n")

tryCatch({
  GDCdownload(query_exp, method = "api", files.per.chunk = 10)
  cat("下载完成！\n")
}, error = function(e) {
  cat("API下载失败，尝试使用client方式...\n")
  GDCdownload(query_exp, method = "client")
})

# 准备数据
cat("\n正在整理表达数据...\n")
data_exp <- GDCprepare(query_exp)

# ======================== 提取表达矩阵 ========================
cat("正在提取表达矩阵...\n")

# 获取counts矩阵
counts_matrix <- assay(data_exp, "unstranded")  # raw counts

# 获取TPM矩阵（如果可用）
if ("tpm_unstrand" %in% assayNames(data_exp)) {
  tpm_matrix <- assay(data_exp, "tpm_unstrand")
} else {
  tpm_matrix <- NULL
  cat("注意：TPM数据不可用，将使用counts数据\n")
}

# 获取基因信息
gene_info <- rowData(data_exp)
gene_info <- as.data.frame(gene_info)

# 将Ensembl ID转换为Gene Symbol
cat("正在转换基因ID...\n")
gene_symbol <- gene_info$gene_name
names(gene_symbol) <- rownames(counts_matrix)

# 去除重复基因（取表达量最高的）
counts_df <- as.data.frame(counts_matrix)
counts_df$gene_symbol <- gene_symbol[rownames(counts_df)]
counts_df <- counts_df[!is.na(counts_df$gene_symbol) & counts_df$gene_symbol != "", ]

# 按基因symbol聚合（取均值）
counts_agg <- counts_df %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), mean)) %>%
  as.data.frame()

rownames(counts_agg) <- counts_agg$gene_symbol
counts_agg$gene_symbol <- NULL
expMatrix <- as.matrix(counts_agg)

cat("表达矩阵维度：", nrow(expMatrix), "个基因 x", ncol(expMatrix), "个样本\n")

# ======================== 下载临床数据 ========================
cat("\n正在下载临床数据...\n")

clinical <- GDCquery_clinic(project = project, type = "clinical")

# 查看临床数据列名
cat("\n临床数据包含以下信息：\n")
print(colnames(clinical))

# ======================== 提取MG分组信息 ========================
cat("\n正在提取MG分组信息...\n")

# TCGA-THYM的MG信息通常在特定字段中
# 需要检查实际的列名

# 尝试查找MG相关列
mg_cols <- grep("myasthenia|MG|autoimmune", colnames(clinical), ignore.case = TRUE, value = TRUE)
cat("找到的MG相关列：", paste(mg_cols, collapse = ", "), "\n")

# 如果没有直接的MG列，尝试从其他来源获取
if (length(mg_cols) == 0) {
  cat("\n临床数据中没有直接的MG列，尝试从biospecimen数据获取...\n")

  # 从样本信息中获取
  sample_info <- colData(data_exp)
  sample_info <- as.data.frame(sample_info)

  # 查看样本信息列
  mg_cols_sample <- grep("myasthenia|MG|autoimmune", colnames(sample_info), ignore.case = TRUE, value = TRUE)

  if (length(mg_cols_sample) > 0) {
    cat("在样本信息中找到MG列：", paste(mg_cols_sample, collapse = ", "), "\n")
  }
}

# ======================== 创建分组变量 ========================
# 注意：TCGA的MG信息可能需要从补充文件获取
# 这里先创建样本信息框架

sample_info <- colData(data_exp)
sample_info <- as.data.frame(sample_info)

# 简化样本ID（去除多余部分）
sample_info$sample_short <- substr(sample_info$barcode, 1, 15)

# 保存样本信息
phenoData <- sample_info[, c("barcode", "sample_short", "sample_type",
                              "age_at_index", "gender", "race",
                              "ajcc_pathologic_stage", "primary_diagnosis")]

# ======================== 处理MG分组 ========================
# TCGA-THYM的MG信息通常在 "myasthenia_gravis_history" 或类似字段
# 如果数据中没有，需要从文献或supplementary获取

# 检查是否有MG相关列
if ("myasthenia_gravis_associated_autoimmune_disease" %in% colnames(sample_info)) {
  phenoData$MG_status <- sample_info$myasthenia_gravis_associated_autoimmune_disease
} else if ("history_other_malignancy" %in% colnames(sample_info)) {
  # 尝试其他可能的列名
  cat("\n注意：需要手动添加MG分组信息\n")
  cat("请参考TCGA论文或GDC补充数据获取MG状态\n")
  phenoData$MG_status <- NA
}

# ======================== 备选方案：从cBioPortal获取MG信息 ========================
cat("\n========================================\n")
cat("备选方案：从cBioPortal下载带MG信息的临床数据\n")
cat("========================================\n")

cat("
如果上述方法无法获取MG分组，请手动下载：

1. 访问 https://www.cbioportal.org/
2. 搜索 'Thymoma (TCGA, PanCancer Atlas)'
3. 点击 'Clinical Data' 标签
4. 下载临床数据CSV文件
5. 查找 'Myasthenia Gravis' 相关列

或者使用以下R代码从cBioPortal获取：
")

cat("
# 安装cBioPortal包
# install.packages('cBioPortalData')
# library(cBioPortalData)
#
# cbio <- cBioPortal()
# studies <- getStudies(cbio)
# thym_study <- studies[grep('thym', studies$studyId, ignore.case = TRUE), ]
#
# # 获取临床数据
# clinical_cbio <- clinicalData(cbio, studyId = 'thym_tcga_pan_can_atlas_2018')
")

# ======================== 保存数据 ========================
cat("\n正在保存数据...\n")

# 保存表达矩阵
save(expMatrix, file = "TCGA_THYM_expression.RData")
cat("  - 表达矩阵已保存：TCGA_THYM_expression.RData\n")

# 保存临床数据
save(phenoData, clinical, file = "TCGA_THYM_clinical.RData")
cat("  - 临床数据已保存：TCGA_THYM_clinical.RData\n")

# 保存为CSV便于查看
write.csv(phenoData, file = "TCGA_THYM_phenoData.csv", row.names = FALSE)
write.csv(clinical, file = "TCGA_THYM_clinical_full.csv", row.names = FALSE)
cat("  - CSV文件已保存\n")

# 保存基因信息
write.csv(gene_info, file = "TCGA_THYM_gene_info.csv", row.names = TRUE)

# ======================== 数据概览 ========================
cat("\n========================================\n")
cat("           数据下载完成！\n")
cat("========================================\n")
cat("\n数据概览：\n")
cat("  - 表达矩阵：", nrow(expMatrix), "个基因 x", ncol(expMatrix), "个样本\n")
cat("  - 临床样本：", nrow(phenoData), "例\n")

# 显示分组信息（如果有）
if ("MG_status" %in% colnames(phenoData) && !all(is.na(phenoData$MG_status))) {
  cat("\n  MG分组：\n")
  print(table(phenoData$MG_status))
} else {
  cat("\n  ⚠️ 注意：MG分组信息需要手动补充\n")
  cat("  请参考下一步的说明获取MG状态\n")
}

cat("\n输出文件：\n")
cat("  - TCGA_THYM_expression.RData（表达矩阵）\n")
cat("  - TCGA_THYM_clinical.RData（临床数据）\n")
cat("  - TCGA_THYM_phenoData.csv（样本信息CSV）\n")
cat("  - TCGA_THYM_clinical_full.csv（完整临床数据CSV）\n")
cat("  - TCGA_THYM_gene_info.csv（基因信息）\n")

cat("\n下一步：\n")
cat("  1. 检查phenoData.csv中的MG分组信息\n")
cat("  2. 如需补充MG信息，运行 step01b_add_MG_status.R\n")
cat("  3. 确认分组后，运行 step02_diff_analysis.R\n")
