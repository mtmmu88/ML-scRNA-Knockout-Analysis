#===============================================================================
# 步骤3：机器学习建模
# 使用114种机器学习方法评估MDH2和铁死亡基因的诊断价值
#===============================================================================

library(glmnet)
library(randomForestSRC)
library(gbm)
library(xgboost)
library(e1071)
library(caret)
library(mboost)
library(plsRglm)
library(MASS)
library(pROC)
library(ComplexHeatmap)
library(RColorBrewer)

#======================== 参数设置 ========================
workDir <- "./03_ML"
diffDir <- "./02_diff"

# 分类变量名
classVar <- "Type"

# 最小特征数阈值
min.selected.var <- 3

#======================== 加载工具函数 ========================
# 复制原代码的ML函数
source("../21.ML/refer.ML.R")

#======================== 创建目录 ========================
if (!dir.exists(workDir)) {
  dir.create(workDir, recursive = TRUE)
}
setwd(workDir)

#======================== 加载数据 ========================
cat("正在加载差异分析结果...\n")
load(file.path(diffDir, "diff_analysis_results.RData"))

#======================== 准备机器学习数据 ========================
cat("正在准备机器学习数据...\n")

# 使用目标基因 + 差异基因作为特征
# 铁死亡基因
ferroptosisGenes <- c(
  "ACSL4", "LPCAT3", "ALOX15", "TFRC", "SLC7A11", "GPX4",
  "NCOA4", "FTH1", "FTL", "IREB2", "HMOX1", "NFS1",
  "SLC3A2", "GSS", "GCLC", "GCLM", "NFE2L2", "KEAP1",
  "VDAC2", "VDAC3", "CISD1", "CS", "ACO2"
)

# 核心目标基因
coreGenes <- c("MDH2", "SIRT5")

# 合并特征基因
featureGenes <- unique(c(coreGenes, ferroptosisGenes))
featureGenes <- featureGenes[featureGenes %in% rownames(geneMatrix)]
cat("使用特征基因数量：", length(featureGenes), "\n")

# 准备表达矩阵
expr_data <- t(geneMatrix[featureGenes, ])
expr_data <- as.data.frame(expr_data)

# 添加分组信息
expr_data$Type <- ifelse(group == "Control", 0, 1)

# 划分训练集和测试集（70% vs 30%）
set.seed(123456)
trainIndex <- createDataPartition(expr_data$Type, p = 0.7, list = FALSE)
Train_data <- expr_data[trainIndex, ]
Test_data <- expr_data[-trainIndex, ]

cat("训练集样本：", nrow(Train_data), "\n")
cat("测试集样本：", nrow(Test_data), "\n")

# 分离特征和标签
Train_expr <- as.matrix(Train_data[, 1:(ncol(Train_data)-1)])
Train_class <- data.frame(Type = Train_data$Type)
rownames(Train_class) <- rownames(Train_data)

Test_expr <- as.matrix(Test_data[, 1:(ncol(Test_data)-1)])
Test_class <- data.frame(Type = Test_data$Type, Cohort = "Test")
rownames(Test_class) <- rownames(Test_data)

# 数据标准化
Train_set <- scale(Train_expr, center = TRUE, scale = TRUE)
Test_set <- scale(Test_expr,
                  center = attr(Train_set, "scaled:center"),
                  scale = attr(Train_set, "scaled:scale"))

#======================== 机器学习方法列表 ========================
# 简化版：使用常用的20种方法（完整版有114种）
methods <- c(
  "Lasso",
  "Ridge",
  "Enet[alpha=0.5]",
  "SVM",
  "RF",
  "GBM",
  "XGBoost",
  "glmBoost",
  "Stepglm[forward]",
  "Stepglm[backward]",
  "Stepglm[both]",
  "LDA",
  "NaiveBayes",
  "Lasso+SVM",
  "Lasso+RF",
  "RF+SVM",
  "glmBoost+Lasso",
  "RF+XGBoost",
  "Lasso+XGBoost",
  "glmBoost+RF"
)

#======================== 训练模型 ========================
cat("\n========== 开始训练机器学习模型 ==========\n")

# 第一阶段：特征选择
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))

preTrain.var <- list()
set.seed(123456)
for (method in preTrain.method) {
  cat("特征选择：", method, "\n")
  tryCatch({
    preTrain.var[[method]] <- RunML(
      method = method,
      Train_set = Train_set,
      Train_label = Train_class,
      mode = "Variable",
      classVar = classVar
    )
  }, error = function(e) {
    cat("  跳过（错误）\n")
  })
}
preTrain.var[["simple"]] <- colnames(Train_set)

# 第二阶段：模型训练
model <- list()
set.seed(123456)
Train_set_bk <- Train_set

for (method in methods) {
  cat(match(method, methods), "/", length(methods), ":", method, "\n")
  method_name <- method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)

  Variable <- preTrain.var[[method[1]]]
  if (is.null(Variable) || length(Variable) == 0) {
    Variable <- colnames(Train_set_bk)
  }

  Train_set <- Train_set_bk[, Variable, drop = FALSE]
  Train_label <- Train_class

  tryCatch({
    model[[method_name]] <- RunML(
      method = method[2],
      Train_set = Train_set,
      Train_label = Train_label,
      mode = "Model",
      classVar = classVar
    )

    # 检查特征数量
    if (length(ExtractVar(model[[method_name]])) <= min.selected.var) {
      model[[method_name]] <- NULL
    }
  }, error = function(e) {
    cat("  跳过（错误）:", conditionMessage(e), "\n")
  })
}
Train_set <- Train_set_bk

# 保存模型
saveRDS(model, "model.MLmodel.rds")
cat("\n模型已保存\n")

#======================== 模型评估 ========================
cat("\n========== 模型评估 ==========\n")

methodsValid <- names(model)
cat("有效模型数量：", length(methodsValid), "\n")

# 计算AUC
AUC_list <- list()
for (method in methodsValid) {
  tryCatch({
    AUC_list[[method]] <- RunEval(
      fit = model[[method]],
      Test_set = Test_set,
      Test_label = Test_class,
      Train_set = Train_set,
      Train_label = Train_class,
      Train_name = "Train",
      cohortVar = "Cohort",
      classVar = classVar
    )
  }, error = function(e) {
    cat("评估失败：", method, "\n")
  })
}

AUC_mat <- do.call(rbind, AUC_list)
AUC_mat <- as.data.frame(AUC_mat)

# 按平均AUC排序
avg_AUC <- rowMeans(AUC_mat, na.rm = TRUE)
AUC_mat <- AUC_mat[order(avg_AUC, decreasing = TRUE), ]

# 保存AUC结果
write.table(cbind(Method = rownames(AUC_mat), AUC_mat),
            file = "model.AUCmatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#======================== 输出最佳模型信息 ========================
cat("\n========== 最佳模型 ==========\n")
bestModel <- rownames(AUC_mat)[1]
cat("最佳模型：", bestModel, "\n")
cat("训练集AUC：", round(AUC_mat[1, "Train"], 3), "\n")
cat("测试集AUC：", round(AUC_mat[1, "Test"], 3), "\n")

# 最佳模型使用的基因
bestFeatures <- ExtractVar(model[[bestModel]])
cat("\n最佳模型使用的基因：\n")
print(bestFeatures)
write.table(bestFeatures, file = "best_model_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#======================== MDH2单独的诊断价值 ========================
cat("\n========== MDH2 单独诊断价值 ==========\n")
if ("MDH2" %in% colnames(Train_set_bk)) {
  MDH2_train <- Train_set_bk[, "MDH2"]
  MDH2_test <- Test_set[, "MDH2"]

  roc_train <- roc(Train_class$Type, MDH2_train, quiet = TRUE)
  roc_test <- roc(Test_class$Type, MDH2_test, quiet = TRUE)

  cat("MDH2 训练集 AUC：", round(auc(roc_train), 3), "\n")
  cat("MDH2 测试集 AUC：", round(auc(roc_test), 3), "\n")

  # 绑制ROC曲线
  pdf("MDH2_ROC.pdf", width = 6, height = 6)
  plot(roc_test, main = "MDH2 ROC Curve", col = "red", lwd = 2)
  legend("bottomright",
         legend = paste0("AUC = ", round(auc(roc_test), 3)),
         col = "red", lwd = 2)
  dev.off()
}

#======================== 绘制AUC热图 ========================
cat("\n正在绑制AUC热图...\n")

# 准备热图数据
AUC_matrix <- as.matrix(AUC_mat)
avg_AUC_sorted <- rowMeans(AUC_matrix, na.rm = TRUE)

# 颜色
CohortCol <- c("Train" = "red", "Test" = "blue")

# 绘制热图
pdf("model.AUCheatmap.pdf", width = 8, height = max(6, nrow(AUC_matrix) * 0.3))
hm <- Heatmap(
  AUC_matrix,
  name = "AUC",
  col = c("#4195C1", "#FFFFFF", "#CB5746"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.text(sprintf("%.3f", AUC_matrix[i, j]), x, y, gp = gpar(fontsize = 8))
  }
)
draw(hm)
dev.off()

#======================== 保存结果 ========================
save(model, AUC_mat, bestModel, bestFeatures,
     file = "ML_results.RData")

cat("\n========== 机器学习分析完成 ==========\n")
cat("输出文件：\n")
cat("  - model.MLmodel.rds（所有模型）\n")
cat("  - model.AUCmatrix.txt（AUC结果）\n")
cat("  - model.AUCheatmap.pdf（AUC热图）\n")
cat("  - best_model_genes.txt（最佳模型基因）\n")
cat("  - MDH2_ROC.pdf（MDH2 ROC曲线）\n")
cat("  - ML_results.RData（R数据）\n")
