#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤5：机器学习建模（114种方法）
#
# 参考原代码：21.ML/toxicityKnk21.model.R
#===============================================================================

#加载包
library(openxlsx)
library(glmnet)
library(randomForestSRC)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(MASS)
library(plsRglm)
library(xgboost)
library(pROC)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#======================== 参数设置 ========================
classVar <- "Type"
min.selected.var <- 3

#设置工作目录
setwd("D:/MDH2_HF_full")         #Windows用户
#setwd("~/MDH2_HF_full")         #Mac/Linux用户

#加载ML函数（修改为你的路径）
source("../21.ML/refer.ML.R")

#======================== 加载数据 ========================
cat("正在加载数据...\n")
load("diff.RData")

#======================== 准备特征基因 ========================
#目标基因 + 铁死亡基因
targetGenes <- c("MDH2", "SIRT5")
ferroptosisGenes <- c(
  "ACSL4", "LPCAT3", "ALOX15", "TFRC", "NCOA4",
  "SLC7A11", "GPX4", "FTH1", "FTL", "NFE2L2",
  "VDAC2", "VDAC3", "CS", "ACO2", "HMOX1",
  "IREB2", "NFS1", "SLC3A2", "GSS", "GCLC"
)

featureGenes <- unique(c(targetGenes, ferroptosisGenes))
featureGenes <- featureGenes[featureGenes %in% rownames(data)]
cat("特征基因数：", length(featureGenes), "\n")

#======================== 准备训练/测试数据 ========================
#提取特征基因表达
expr_data <- t(data[featureGenes, ])
expr_data <- as.data.frame(expr_data)
expr_data$Type <- ifelse(Type == "Control", 0, 1)

#按数据集划分：GSE59867作为训练集，其他作为测试集
trainIdx <- grepl("GSE59867", rownames(expr_data))
Train_data <- expr_data[trainIdx, ]
Test_data <- expr_data[!trainIdx, ]

cat("训练集（GSE59867）：", nrow(Train_data), "样本\n")
cat("测试集（其他）：", nrow(Test_data), "样本\n")

#分离特征和标签
Train_expr <- as.matrix(Train_data[, 1:(ncol(Train_data) - 1)])
Train_class <- data.frame(Type = Train_data$Type)
rownames(Train_class) <- rownames(Train_data)

Test_expr <- as.matrix(Test_data[, 1:(ncol(Test_data) - 1)])
Test_class <- data.frame(
  Type = Test_data$Type,
  Cohort = gsub("_.*", "", rownames(Test_data))
)
rownames(Test_class) <- rownames(Test_data)

#数据标准化
Train_set <- scale(Train_expr, center = TRUE, scale = TRUE)
Test_set <- scale(Test_expr,
                  center = attr(Train_set, "scaled:center"),
                  scale = attr(Train_set, "scaled:scale"))

#======================== 114种机器学习方法 ========================
methods <- c(
  "Lasso", "Ridge",
  "Enet[alpha=0.1]", "Enet[alpha=0.2]", "Enet[alpha=0.3]",
  "Enet[alpha=0.4]", "Enet[alpha=0.5]", "Enet[alpha=0.6]",
  "Enet[alpha=0.7]", "Enet[alpha=0.8]", "Enet[alpha=0.9]",
  "SVM", "RF", "GBM", "XGBoost", "LDA", "NaiveBayes",
  "glmBoost", "plsRglm",
  "Stepglm[forward]", "Stepglm[backward]", "Stepglm[both]",
  #组合方法
  "Lasso+SVM", "Lasso+RF", "Lasso+Ridge", "Lasso+GBM", "Lasso+XGBoost",
  "Lasso+glmBoost", "Lasso+plsRglm", "Lasso+LDA", "Lasso+NaiveBayes",
  "Lasso+Stepglm[forward]", "Lasso+Stepglm[backward]", "Lasso+Stepglm[both]",
  "Ridge+SVM", "Ridge+RF",
  "RF+SVM", "RF+Ridge", "RF+GBM", "RF+XGBoost", "RF+LDA",
  "RF+Enet[alpha=0.1]", "RF+Enet[alpha=0.5]", "RF+Enet[alpha=0.9]",
  "RF+Stepglm[forward]", "RF+Stepglm[backward]", "RF+Stepglm[both]",
  "glmBoost+SVM", "glmBoost+RF", "glmBoost+Ridge", "glmBoost+Lasso",
  "glmBoost+GBM", "glmBoost+XGBoost", "glmBoost+LDA", "glmBoost+NaiveBayes",
  "glmBoost+Enet[alpha=0.1]", "glmBoost+Enet[alpha=0.5]", "glmBoost+Enet[alpha=0.9]",
  "glmBoost+Stepglm[forward]", "glmBoost+Stepglm[backward]", "glmBoost+Stepglm[both]",
  "glmBoost+plsRglm",
  "Stepglm[forward]+SVM", "Stepglm[forward]+RF", "Stepglm[forward]+Ridge",
  "Stepglm[forward]+GBM", "Stepglm[forward]+LDA",
  "Stepglm[backward]+SVM", "Stepglm[backward]+RF", "Stepglm[backward]+Ridge",
  "Stepglm[backward]+GBM", "Stepglm[backward]+LDA",
  "Stepglm[both]+SVM", "Stepglm[both]+RF", "Stepglm[both]+Ridge",
  "Stepglm[both]+GBM", "Stepglm[both]+LDA", "Stepglm[both]+glmBoost",
  "Stepglm[both]+Enet[alpha=0.1]", "Stepglm[both]+Enet[alpha=0.5]",
  "Stepglm[both]+Lasso", "Stepglm[both]+plsRglm"
)

cat("\n使用", length(methods), "种机器学习方法\n")

#======================== 第一阶段：特征选择 ========================
cat("\n========== 第一阶段：特征选择 ==========\n")

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
    cat("  [跳过]\n")
  })
}
preTrain.var[["simple"]] <- colnames(Train_set)

#======================== 第二阶段：模型训练 ========================
cat("\n========== 第二阶段：模型训练 ==========\n")

model <- list()
set.seed(123456)
Train_set_bk <- Train_set

for (i in 1:length(methods)) {
  method <- methods[i]
  cat(i, "/", length(methods), ":", method, "\n")
  method_name <- method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)

  Variable <- preTrain.var[[method[1]]]
  if (is.null(Variable) || length(Variable) == 0) {
    Variable <- colnames(Train_set_bk)
  }

  Train_set <- Train_set_bk[, Variable, drop = FALSE]

  tryCatch({
    model[[method_name]] <- RunML(
      method = method[2],
      Train_set = Train_set,
      Train_label = Train_class,
      mode = "Model",
      classVar = classVar
    )

    if (length(ExtractVar(model[[method_name]])) <= min.selected.var) {
      model[[method_name]] <- NULL
    }
  }, error = function(e) {
    #静默跳过
  })
}
Train_set <- Train_set_bk

#保存模型
saveRDS(model, "ML.model.rds")
cat("\n有效模型：", length(model), "个\n")

#======================== 模型评估 ========================
cat("\n========== 模型评估 ==========\n")

methodsValid <- names(model)
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
    #静默跳过
  })
}

AUC_mat <- do.call(rbind, AUC_list)
AUC_mat <- as.data.frame(AUC_mat)

#按平均AUC排序
avg_AUC <- rowMeans(AUC_mat, na.rm = TRUE)
AUC_mat <- AUC_mat[order(avg_AUC, decreasing = TRUE), ]

write.table(cbind(Method = rownames(AUC_mat), AUC_mat),
            file = "ML.AUCmatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#======================== 最佳模型 ========================
cat("\n========== 最佳模型 ==========\n")
bestModel <- rownames(AUC_mat)[1]
cat("最佳模型：", bestModel, "\n")
cat("平均AUC：", round(mean(as.numeric(AUC_mat[1, ]), na.rm = TRUE), 3), "\n")

bestFeatures <- ExtractVar(model[[bestModel]])
cat("使用基因：", paste(bestFeatures, collapse = ", "), "\n")
write.table(bestFeatures, file = "ML.bestModelGenes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#======================== MDH2单独诊断价值 ========================
cat("\n========== MDH2单独诊断价值 ==========\n")
if ("MDH2" %in% colnames(Train_set_bk)) {
  roc_train <- roc(Train_class$Type, Train_set_bk[, "MDH2"], quiet = TRUE)
  roc_test <- roc(Test_class$Type, Test_set[, "MDH2"], quiet = TRUE)

  cat("MDH2 训练集 AUC：", round(auc(roc_train), 3), "\n")
  cat("MDH2 测试集 AUC：", round(auc(roc_test), 3), "\n")

  pdf("ML.MDH2_ROC.pdf", width = 6, height = 6)
  plot(roc_test, main = "MDH2 Diagnostic Value", col = "red", lwd = 2,
       print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.2)
  dev.off()
}

#======================== AUC热图 ========================
cat("\n正在绘制AUC热图...\n")
AUC_matrix <- as.matrix(AUC_mat)

pdf("ML.AUCheatmap.pdf", width = 8, height = max(10, nrow(AUC_matrix) * 0.25))
hm <- Heatmap(
  AUC_matrix,
  name = "AUC",
  col = colorRamp2(c(0.5, 0.75, 1), c("#4195C1", "#FFFFFF", "#CB5746")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 6),
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.text(sprintf("%.2f", AUC_matrix[i, j]), x, y, gp = gpar(fontsize = 5))
  }
)
draw(hm)
dev.off()

#======================== 保存结果 ========================
save(model, AUC_mat, bestModel, bestFeatures,
     Train_set, Test_set, Train_class, Test_class,
     file = "ML.results.RData")

cat("\n========== 机器学习完成 ==========\n")
cat("输出文件：\n")
cat("  - ML.model.rds\n")
cat("  - ML.AUCmatrix.txt\n")
cat("  - ML.AUCheatmap.pdf\n")
cat("  - ML.MDH2_ROC.pdf\n")
cat("  - ML.bestModelGenes.txt\n")
cat("\n请运行下一步：step06_bindao_scKnockout.R\n")


######生信自学网: https://www.biowolf.cn/
######MDH2心衰研究完整版
