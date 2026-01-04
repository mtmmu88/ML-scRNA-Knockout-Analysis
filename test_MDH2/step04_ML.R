#===============================================================================
# MDH2心衰研究 - 测试版本
# 步骤4：机器学习建模
#
# 参考原代码：21.ML/toxicityKnk21.model.R
# 测试版本：先用10种方法，跑通后可扩展到114种
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
library(pROC)
library(ComplexHeatmap)
library(RColorBrewer)

#======================== 参数设置（可修改）========================
geoID <- "GSE26887"
classVar <- "Type"              #分类变量名
min.selected.var <- 3           #最小特征数阈值

#设置工作目录（修改为你的路径）
setwd("D:/test_MDH2")           #Windows用户修改这里
#setwd("~/test_MDH2")           #Mac/Linux用户用这个

#======================== 加载ML函数 ========================
#从原代码目录加载ML函数（修改为你的路径）
source("../21.ML/refer.ML.R")

#======================== 加载数据 ========================
cat("正在加载差异分析结果...\n")
load(paste0(geoID, ".bindao_diff.RData"))

#======================== 准备特征基因 ========================
#目标基因 + 铁死亡基因作为特征
targetGenes <- c("MDH2", "SIRT5")
ferroptosisGenes <- c(
  "ACSL4", "LPCAT3", "ALOX15", "TFRC", "NCOA4",
  "SLC7A11", "GPX4", "FTH1", "FTL", "NFE2L2",
  "VDAC2", "VDAC3", "CS", "ACO2", "HMOX1"
)

#合并特征基因
featureGenes <- unique(c(targetGenes, ferroptosisGenes))
featureGenes <- featureGenes[featureGenes %in% rownames(data)]
cat("使用特征基因数量：", length(featureGenes), "\n")
cat("特征基因：", paste(featureGenes, collapse = ", "), "\n")

#======================== 准备训练数据 ========================
#提取特征基因表达
expr_data <- t(data[featureGenes, ])
expr_data <- as.data.frame(expr_data)

#添加分组标签（0=Control, 1=Disease）
expr_data$Type <- ifelse(Type == "Control", 0, 1)

#划分训练集和测试集（70% vs 30%）
set.seed(123456)
trainIndex <- createDataPartition(expr_data$Type, p = 0.7, list = FALSE)
Train_data <- expr_data[trainIndex, ]
Test_data <- expr_data[-trainIndex, ]

cat("\n训练集样本：", nrow(Train_data), "\n")
cat("测试集样本：", nrow(Test_data), "\n")

#分离特征和标签
Train_expr <- as.matrix(Train_data[, 1:(ncol(Train_data) - 1)])
Train_class <- data.frame(Type = Train_data$Type)
rownames(Train_class) <- rownames(Train_data)

Test_expr <- as.matrix(Test_data[, 1:(ncol(Test_data) - 1)])
Test_class <- data.frame(Type = Test_data$Type, Cohort = "Test")
rownames(Test_class) <- rownames(Test_data)

#数据标准化
Train_set <- scale(Train_expr, center = TRUE, scale = TRUE)
Test_set <- scale(Test_expr,
                  center = attr(Train_set, "scaled:center"),
                  scale = attr(Train_set, "scaled:scale"))

#======================== 机器学习方法（测试版：10种）========================
#测试版本用10种方法，正式分析可扩展到114种
methods <- c(
  "Lasso",
  "Ridge",
  "Enet[alpha=0.5]",
  "SVM",
  "RF",
  "glmBoost",
  "Stepglm[forward]",
  "Stepglm[backward]",
  "LDA",
  "Lasso+RF"
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
    cat("  [跳过] 错误：", conditionMessage(e), "\n")
  })
}
preTrain.var[["simple"]] <- colnames(Train_set)

#======================== 第二阶段：模型训练 ========================
cat("\n========== 第二阶段：模型训练 ==========\n")

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

    #检查特征数量
    if (length(ExtractVar(model[[method_name]])) <= min.selected.var) {
      model[[method_name]] <- NULL
      cat("  [移除] 特征数不足\n")
    }
  }, error = function(e) {
    cat("  [跳过] 错误：", conditionMessage(e), "\n")
  })
}
Train_set <- Train_set_bk

#保存模型
saveRDS(model, "ML.model.rds")
cat("\n模型已保存：ML.model.rds\n")

#======================== 模型评估 ========================
cat("\n========== 模型评估 ==========\n")

methodsValid <- names(model)
cat("有效模型数量：", length(methodsValid), "\n")

#计算AUC
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

#按平均AUC排序
avg_AUC <- rowMeans(AUC_mat, na.rm = TRUE)
AUC_mat <- AUC_mat[order(avg_AUC, decreasing = TRUE), ]

#保存AUC结果
write.table(cbind(Method = rownames(AUC_mat), AUC_mat),
            file = "ML.AUCmatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#======================== 输出最佳模型 ========================
cat("\n========== 最佳模型 ==========\n")
bestModel <- rownames(AUC_mat)[1]
cat("最佳模型：", bestModel, "\n")
cat("训练集AUC：", round(AUC_mat[1, "Train"], 3), "\n")
cat("测试集AUC：", round(AUC_mat[1, "Test"], 3), "\n")

#最佳模型的特征
bestFeatures <- ExtractVar(model[[bestModel]])
cat("\n最佳模型使用的基因：\n")
print(bestFeatures)
write.table(bestFeatures, file = "ML.bestModelGenes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#======================== MDH2单独诊断价值 ========================
cat("\n========== MDH2单独诊断价值 ==========\n")
if ("MDH2" %in% colnames(Train_set_bk)) {
  MDH2_train <- Train_set_bk[, "MDH2"]
  MDH2_test <- Test_set[, "MDH2"]

  roc_train <- roc(Train_class$Type, MDH2_train, quiet = TRUE)
  roc_test <- roc(Test_class$Type, MDH2_test, quiet = TRUE)

  cat("MDH2 训练集 AUC：", round(auc(roc_train), 3), "\n")
  cat("MDH2 测试集 AUC：", round(auc(roc_test), 3), "\n")

  #绑制ROC曲线
  pdf("ML.MDH2_ROC.pdf", width = 6, height = 6)
  plot(roc_test, main = "MDH2 Diagnostic Value", col = "red", lwd = 2,
       print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.2)
  dev.off()
  cat("MDH2 ROC曲线已保存：ML.MDH2_ROC.pdf\n")
}

#======================== 所有模型AUC展示 ========================
cat("\n========== 所有模型AUC ==========\n")
print(round(AUC_mat, 3))

#======================== 绘制AUC热图 ========================
cat("\n正在绘制AUC热图...\n")

AUC_matrix <- as.matrix(AUC_mat)
avg_AUC_sorted <- rowMeans(AUC_matrix, na.rm = TRUE)

pdf("ML.AUCheatmap.pdf", width = 6, height = max(4, nrow(AUC_matrix) * 0.4))
hm <- Heatmap(
  AUC_matrix,
  name = "AUC",
  col = circlize::colorRamp2(c(0.5, 0.75, 1), c("#4195C1", "#FFFFFF", "#CB5746")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.text(sprintf("%.3f", AUC_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)
draw(hm)
dev.off()
cat("AUC热图已保存：ML.AUCheatmap.pdf\n")

#======================== 保存结果 ========================
save(model, AUC_mat, bestModel, bestFeatures, Train_set, Test_set,
     Train_class, Test_class, file = "ML.results.RData")

cat("\n========== 机器学习分析完成 ==========\n")
cat("保存文件：\n")
cat("  - ML.model.rds（所有模型）\n")
cat("  - ML.AUCmatrix.txt（AUC结果）\n")
cat("  - ML.AUCheatmap.pdf（AUC热图）\n")
cat("  - ML.MDH2_ROC.pdf（MDH2 ROC曲线）\n")
cat("  - ML.bestModelGenes.txt（最佳模型基因）\n")
cat("  - ML.results.RData（R数据）\n")

cat("\n========== 测试流程完成！==========\n")
cat("如果一切正常，你可以：\n")
cat("  1. 换成大数据集（GSE59867）重新跑\n")
cat("  2. 扩展到114种机器学习方法\n")
cat("  3. 继续做单细胞敲除分析\n")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
