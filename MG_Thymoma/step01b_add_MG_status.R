#===============================================================================
# 步骤1b：补充MG分组信息
# TCGA-THYM的MG状态可能需要从补充来源获取
#===============================================================================

library(dplyr)

# ======================== 设置目录 ========================
workDir <- "./01_TCGA_data"
setwd(workDir)

# 加载之前下载的数据
load("TCGA_THYM_clinical.RData")
load("TCGA_THYM_expression.RData")

cat("========================================\n")
cat("  补充TCGA-THYM的MG分组信息\n")
cat("========================================\n\n")

# ======================== 方法1：检查已有临床数据 ========================
cat("方法1：检查已有临床数据中的MG信息...\n\n")

# 查看所有可能包含MG信息的列
possible_mg_cols <- c(
  "myasthenia_gravis",
  "myasthenia_gravis_associated",
  "myasthenia_gravis_associated_autoimmune_disease",
  "history_myasthenia_gravis",
  "mg_status",
  "autoimmune_disease"
)

for (col in possible_mg_cols) {
  if (col %in% colnames(clinical)) {
    cat("找到列：", col, "\n")
    print(table(clinical[[col]], useNA = "ifany"))
    cat("\n")
  }
}

# 在phenoData中查找
cat("\nphenoData中的列：\n")
print(colnames(phenoData))

# ======================== 方法2：从GDC获取补充临床数据 ========================
cat("\n方法2：从GDC获取更详细的临床数据...\n")

library(TCGAbiolinks)

# 获取所有类型的临床数据
clinical_supplement <- tryCatch({
  GDCquery_clinic(project = "TCGA-THYM", type = "clinical")
}, error = function(e) {
  cat("获取临床数据失败：", e$message, "\n")
  NULL
})

if (!is.null(clinical_supplement)) {
  # 查找MG相关列
  mg_cols <- grep("myasthenia|autoimmune|MG", colnames(clinical_supplement),
                  ignore.case = TRUE, value = TRUE)
  if (length(mg_cols) > 0) {
    cat("找到MG相关列：\n")
    for (col in mg_cols) {
      cat("\n", col, ":\n")
      print(table(clinical_supplement[[col]], useNA = "ifany"))
    }
  }
}

# ======================== 方法3：使用已知的MG患者列表 ========================
cat("\n========================================\n")
cat("方法3：根据文献提供MG分组\n")
cat("========================================\n")

cat("
根据TCGA-THYM发表的研究（Nature Communications 2022等）：
- 总样本数：116例
- MG阳性：34例（约29%）
- MG阴性：82例

如果自动获取失败，需要手动添加MG信息。
以下是从cBioPortal或文献中获取MG状态的方法：
")

# ======================== 方法4：手动输入（基于cBioPortal数据）========================
cat("\n方法4：从cBioPortal下载临床数据\n")

cat("
请按以下步骤操作：

1. 访问 https://www.cbioportal.org/study/clinicalData?id=thym_tcga_pan_can_atlas_2018
2. 点击 'Download' 下载临床数据
3. 在Excel中打开，查找 'Myasthenia Gravis' 列
4. 将数据保存为 'cbioportal_clinical.csv' 放入当前目录
")

# 检查是否有cBioPortal下载的数据
if (file.exists("cbioportal_clinical.csv")) {
  cat("\n发现cBioPortal临床数据，正在读取...\n")

  cbio_clinical <- read.csv("cbioportal_clinical.csv", stringsAsFactors = FALSE)

  # 查找MG列
  mg_col <- grep("myasthenia|MG", colnames(cbio_clinical), ignore.case = TRUE, value = TRUE)

  if (length(mg_col) > 0) {
    cat("找到MG列：", mg_col[1], "\n")
    print(table(cbio_clinical[[mg_col[1]]], useNA = "ifany"))

    # 合并到phenoData
    # 需要匹配样本ID
  }
}

# ======================== 创建模拟分组（用于测试）========================
cat("\n========================================\n")
cat("创建分组数据框架\n")
cat("========================================\n")

# 获取样本ID
sample_ids <- colnames(expMatrix)
cat("表达矩阵中的样本数：", length(sample_ids), "\n")

# 创建分组数据框
group_df <- data.frame(
  sample_id = sample_ids,
  sample_short = substr(sample_ids, 1, 15),
  patient_id = substr(sample_ids, 1, 12),
  MG_status = NA,  # 待填充
  stringsAsFactors = FALSE
)

# ======================== 如果有MG信息，填充分组 ========================
# 这里需要根据实际数据源填充

# 示例：如果clinical中有MG信息
if ("submitter_id" %in% colnames(clinical)) {
  # 尝试匹配
  clinical$patient_short <- substr(clinical$submitter_id, 1, 12)

  # 查找可能的MG列
  for (col in colnames(clinical)) {
    if (grepl("myasthenia", col, ignore.case = TRUE)) {
      cat("\n使用列", col, "作为MG状态\n")

      # 创建映射
      mg_map <- clinical[, c("patient_short", col)]
      colnames(mg_map)[2] <- "MG_raw"

      # 合并
      group_df <- merge(group_df, mg_map, by.x = "patient_id", by.y = "patient_short", all.x = TRUE)

      # 标准化MG状态
      group_df$MG_status <- case_when(
        grepl("yes|true|positive|1", group_df$MG_raw, ignore.case = TRUE) ~ "MG_positive",
        grepl("no|false|negative|0", group_df$MG_raw, ignore.case = TRUE) ~ "MG_negative",
        TRUE ~ NA_character_
      )

      break
    }
  }
}

# ======================== 保存分组信息 ========================
cat("\n当前分组状态：\n")
print(table(group_df$MG_status, useNA = "ifany"))

# 保存
write.csv(group_df, file = "sample_group_info.csv", row.names = FALSE)
save(group_df, file = "sample_group_info.RData")

cat("\n分组信息已保存至：\n")
cat("  - sample_group_info.csv\n")
cat("  - sample_group_info.RData\n")

# ======================== 使用说明 ========================
cat("\n========================================\n")
cat("           下一步操作说明\n")
cat("========================================\n")

cat("
如果MG_status列全是NA，请按以下方式补充：

【选项1】从cBioPortal网页下载
1. 访问：https://www.cbioportal.org/study/clinicalData?id=thym_tcga_pan_can_atlas_2018
2. 下载临床数据表格
3. 找到MG相关列，整理后添加到 sample_group_info.csv

【选项2】从文献Supplementary Table获取
参考文献：
- Radovich M, et al. The Integrated Genomic Landscape of Thymic Epithelial Tumors. Cancer Cell. 2018
- Marx A, et al. Nature Communications. 2022

【选项3】手动编辑CSV文件
1. 打开 sample_group_info.csv
2. 在MG_status列填入：MG_positive 或 MG_negative
3. 保存文件
4. 重新运行此脚本的最后部分加载数据

填充完成后，运行 step02_diff_analysis.R
")

# ======================== 验证函数 ========================
cat("\n验证分组数据的函数：\n")

validate_group <- function() {
  if (file.exists("sample_group_info.csv")) {
    df <- read.csv("sample_group_info.csv", stringsAsFactors = FALSE)
    cat("\n当前分组状态：\n")
    print(table(df$MG_status, useNA = "ifany"))

    n_mg_pos <- sum(df$MG_status == "MG_positive", na.rm = TRUE)
    n_mg_neg <- sum(df$MG_status == "MG_negative", na.rm = TRUE)
    n_na <- sum(is.na(df$MG_status))

    cat("\nMG阳性：", n_mg_pos, "例\n")
    cat("MG阴性：", n_mg_neg, "例\n")
    cat("未知：", n_na, "例\n")

    if (n_mg_pos >= 30 && n_mg_neg >= 30) {
      cat("\n✅ 分组数据充足，可以进行后续分析！\n")
      return(TRUE)
    } else {
      cat("\n⚠️ 分组数据不完整，请补充MG状态信息\n")
      return(FALSE)
    }
  } else {
    cat("未找到分组文件\n")
    return(FALSE)
  }
}

# 运行验证
validate_group()
