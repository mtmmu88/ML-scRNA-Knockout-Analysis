#===============================================================================
# MDH2心衰研究 - 测试版本
# 步骤0：安装R包
#
# 研究假设：Sirt5 → MDH2乳酸化 → 线粒体代谢紊乱 → 铁死亡 → 心梗后心衰
#===============================================================================

#设置CRAN镜像（国内用户用清华镜像更快）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

#======================== 基础包 ========================
install.packages("dplyr")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("openxlsx")
install.packages("VennDiagram")
install.packages("ggrepel")
install.packages("RColorBrewer")

#======================== Bioconductor包 ========================
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")      #下载GEO数据
BiocManager::install("limma")         #差异分析
BiocManager::install("sva")           #批次校正
BiocManager::install("org.Hs.eg.db")  #基因注释
BiocManager::install("clusterProfiler") #富集分析
BiocManager::install("ComplexHeatmap")  #复杂热图

#======================== 机器学习包 ========================
install.packages("glmnet")            #Lasso/Ridge/Elastic Net
install.packages("randomForestSRC")   #随机森林
install.packages("gbm")               #梯度提升
install.packages("xgboost")           #XGBoost
install.packages("e1071")             #SVM
install.packages("caret")             #机器学习框架
install.packages("mboost")            #glmBoost
install.packages("plsRglm")           #PLS回归
install.packages("MASS")              #LDA
install.packages("pROC")              #ROC曲线

#======================== 检查安装结果 ========================
cat("\n========== 检查安装结果 ==========\n")
packages <- c("GEOquery", "limma", "dplyr", "ggplot2", "pheatmap",
              "glmnet", "randomForestSRC", "pROC", "ComplexHeatmap")

for (pkg in packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("[OK] ", pkg, "\n"))
  } else {
    cat(paste0("[FAIL] ", pkg, " - 请手动安装\n"))
  }
}

cat("\n安装完成！请运行下一步：step01_download.R\n")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
