#install.packages("clustree")
#install.packages("Seurat")
#install.packages("harmony")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle")

#install.packages("devtools")
#devtools::install_github('immunogenomics/presto')

#install.packages("assertthat")
#install.packages("SCpubr")



#######################01.数据前期处理和矫正#######################
#引用包
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(clustree)
library(harmony)
library(assertthat)
library(SCpubr)

#设置工作目录
workDir="C:\\Users\\Administrator\\Desktop\\toxicityKnk\\25.geneView"
setwd(workDir)

#导入单细胞的对象
load("Seurat.Rdata")

#读取模型基因的文件
geneRT=read.table("optimalModelGenes.txt", header=F, sep="\t", check.names=F)
hubGenes=as.vector(geneRT[,1])
hubGenes=hubGenes[hubGenes %in% rownames(pbmc)]

#绘制散点图
pdf(file="Scatter.pdf", width=12, height=9)
FeaturePlot(object = pbmc, features = hubGenes, cols = c("green", "red"), ncol=4)
dev.off()

#绘制小提琴图
pdf(file="Violin.pdf", width=15, height=12)
VlnPlot(object = pbmc, features = hubGenes, ncol=4)
dev.off()

#绘制气泡图
pdf(file="Bubble.pdf", width=8, height=5)
DotPlot(object = pbmc, features = hubGenes) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


