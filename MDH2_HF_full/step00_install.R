#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤0：安装R包
#
# 研究假设：Sirt5 → MDH2乳酸化 → 线粒体代谢紊乱 → 铁死亡 → 心梗后心衰
# 数据集：GSE59867（心梗后心衰，111例）
#===============================================================================

#设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

#======================== 基础包 ========================
install.packages("dplyr")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("openxlsx")
install.packages("VennDiagram")
install.packages("ggrepel")
install.packages("RColorBrewer")
install.packages("circlize")

#======================== Bioconductor包 ========================
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("sva")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("ComplexHeatmap")
BiocManager::install("DOSE")

#======================== 机器学习包 ========================
install.packages("glmnet")
install.packages("randomForestSRC")
install.packages("gbm")
install.packages("xgboost")
install.packages("e1071")
install.packages("caret")
install.packages("mboost")
install.packages("plsRglm")
install.packages("MASS")
install.packages("pROC")
install.packages("rms")           #列线图

#======================== 单细胞分析包 ========================
install.packages("Seurat")
install.packages("remotes")
remotes::install_github('cailab-tamu/scTenifoldKnk')

#======================== 检查安装结果 ========================
cat("\n========== 检查安装结果 ==========\n")
packages <- c("GEOquery", "limma", "dplyr", "ggplot2", "pheatmap",
              "glmnet", "randomForestSRC", "pROC", "ComplexHeatmap",
              "clusterProfiler", "Seurat", "rms")

for (pkg in packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("[OK] ", pkg, "\n"))
  } else {
    cat(paste0("[FAIL] ", pkg, " - 请手动安装\n"))
  }
}

cat("\n安装完成！请运行下一步：step01_download.R\n")


######生信自学网: https://www.biowolf.cn/
######MDH2心衰研究完整版
