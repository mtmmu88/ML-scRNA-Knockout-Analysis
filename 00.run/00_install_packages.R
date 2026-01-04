#===============================================================================
# MDH2-心衰研究：R包安装脚本
# 研究方向：MDH2乳酸化修饰 - Sirt5 - 线粒体代谢 - 铁死亡 - 心衰/心梗
#===============================================================================

# 设置CRAN镜像（中国用户用清华镜像更快）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# ======================== 基础包 ========================
install.packages(c(
  "dplyr",        # 数据处理
  "ggplot2",      # 绑图
  "pheatmap",     # 热图
  "openxlsx",     # 读写Excel
  "VennDiagram",  # 韦恩图
  "ggrepel",      # 标签避免重叠
  "RColorBrewer"  # 颜色
))

# ======================== GEO数据下载 ========================
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",     # 下载GEO数据
  "limma",        # 差异分析
  "sva",          # 批次校正
  "org.Hs.eg.db", # 人类基因注释
  "clusterProfiler", # 富集分析
  "enrichplot"    # 富集分析可视化
))

# ======================== 机器学习包 ========================
install.packages(c(
  "glmnet",       # Lasso/Ridge/Elastic Net
  "randomForestSRC", # 随机森林
  "gbm",          # 梯度提升
  "xgboost",      # XGBoost
  "e1071",        # SVM
  "caret",        # 机器学习框架
  "mboost",       # glmBoost
  "plsRglm",      # PLS回归
  "MASS",         # LDA
  "pROC"          # ROC曲线
))

# ======================== 单细胞分析 ========================
install.packages("Seurat")
install.packages("remotes")
remotes::install_github('cailab-tamu/scTenifoldKnk')

# ======================== 复杂热图 ========================
BiocManager::install("ComplexHeatmap")

# ======================== 检查安装 ========================
cat("\n========== 检查安装结果 ==========\n")
packages <- c("GEOquery", "limma", "glmnet", "Seurat", "scTenifoldKnk",
              "clusterProfiler", "ComplexHeatmap", "pROC")
for (pkg in packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("✓ ", pkg, " 安装成功\n"))
  } else {
    cat(paste0("✗ ", pkg, " 安装失败\n"))
  }
}

cat("\n所有包安装完成！\n")
