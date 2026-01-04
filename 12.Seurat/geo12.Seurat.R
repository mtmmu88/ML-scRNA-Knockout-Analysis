#install.packages("clustree")
#install.packages("Seurat")
#install.packages("harmony")
#install.packages("gridExtra")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("scrapper")
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
library(gridExtra)

logFCfilter=1           #logFC的过滤条件
adjPvalFilter=0.05     #矫正后的pvalue的过滤条件

#设置工作目录
workDir="C:\\Users\\Administrator\\Desktop\\GEO\\12.Seurat"
setwd(workDir)

#读取数据
dirs=list.dirs(workDir)
dirs_sample=dirs[-1]
names(dirs_sample)=gsub(".+\\/(.+)", "\\1", dirs_sample)
counts <- Read10X(data.dir = dirs_sample)
pbmc = CreateSeuratObject(counts, min.cells=5, min.features=100)

#使用PercentageFeatureSet函数计算线粒体基因的百分比
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#绘制基因特征的小提琴图
pdf(file="01.featureViolin.pdf", width=10, height=6.5)
# 分别绘制3个小提琴图，移除每个图的图例
p1 <- VlnPlot(pbmc, features = "nFeature_RNA", raster = FALSE) + theme(legend.position = "none")
p2 <- VlnPlot(pbmc, features = "nCount_RNA", raster = FALSE) + theme(legend.position = "none")
p3 <- VlnPlot(pbmc, features = "percent.mt", raster = FALSE) + theme(legend.position = "none")
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()
pbmc=subset(x = pbmc, subset = nFeature_RNA > 100 & percent.mt < 10)    #对数据进行过滤

#绘制测序深度相关性的图形
pdf(file="01.featureCor.pdf", width=13, height=7)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#对数据进行标准化
pbmc <- NormalizeData(object=pbmc, normalization.method="LogNormalize", scale.factor=10000)
#提取细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=1500)
#输出特征方差图
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf", width=10, height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()



#######################02.PCA主成分分析#######################
##PCA分析
pbmc=ScaleData(pbmc)        #PCA降维之前的标准预处理步骤
pbmc=RunPCA(object= pbmc, npcs=20, pc.genes=VariableFeatures(object=pbmc))     #PCA分析
pbmc=RunHarmony(pbmc, "orig.ident")

#绘制每个PCA成分的特征基因
pdf(file="02.pcaGene.pdf", width=10, height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures=20)
dev.off()

#绘制PCA图形
pdf(file="02.PCA.pdf", width=7.5, height=5)
DimPlot(object=pbmc, reduction="pca")
dev.off()

#PCA的热图
pdf(file="02.pcaHeatmap.pdf", width=10, height=8)
DimHeatmap(object=pbmc, dims=1:4, cells=500, balanced=TRUE, nfeatures=30, ncol=2)
dev.off()

#得到每个PC的p值分布
pbmc <- JackStraw(object=pbmc, num.replicate=100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object=pbmc, dims=1:20)
dev.off()



#######################03.细胞聚类分析和marker基因#######################
##聚类分析
pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)     #计算邻接距离

#对细胞分组,对细胞标准模块化
pbmc <- FindClusters(pbmc, resolution=seq(0.5, 1.2, by=0.1))
pbmc <- FindClusters(object = pbmc, resolution=0.6)
#输出聚类图形
pdf(file="03.cluster.pdf", width=7.5, height=6)
#pbmc <-RunUMAP(object = pbmc, dims = 1:pcSelect)        #UMAP聚类
#DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)      #UMAP可视化
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)             #TSNE聚类
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)     #TSNE可视化
dev.off()
write.table(pbmc$seurat_clusters,file="03.Cluster.txt",quote=F,sep="\t",col.names=F)

##查找每个聚类的差异基因
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#绘制marker在每个聚类的热图
pdf(file="03.clusterHeatmap.pdf",width=15, height=15)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()



#######################04.SingleR R包注释细胞类型#######################
pbmc_for_SingleR <- GetAssayData(pbmc, layer="data")
clusters<-pbmc@meta.data$seurat_clusters
ref1=get(load("ref_Human_all.RData"))
ref2=get(load("ref_Hematopoietic.RData"))
ref3=get(load("DatabaseImmuneCellExpressionData.Rdata"))
ref4=get(load("BlueprintEncode_bpe.se_human.RData"))
ref5=get(load("HumanPrimaryCellAtlas_hpca.se_human.RData"))
ref6=get(load("MonacoImmuneData.Rdata"))
ref7=get(load("NovershternHematopoieticData.Rdata"))
singler=SingleR(test=pbmc_for_SingleR, ref =list(ref1, ref2, ref3, ref4, ref5, ref6, ref7),
              labels=list(ref1$label.main,ref2$label.main,ref3$label.main,ref4$label.main,ref5$label.main,ref6$label.main,ref7$label.main), clusters = clusters)

singler$labels=gsub("_|-", " ", singler$labels)
singler$labels=gsub("T cells, CD4\\+", "CD4+ T cells", singler$labels)
singler$labels=gsub("T cells, CD8\\+", "CD8+ T cells", singler$labels)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)

#输出细胞注释的结果
cellAnn <- clusterAnn[match(pbmc$seurat_clusters, clusterAnn[,1]), 2]
cellAnnOut=cbind(names(pbmc$seurat_clusters), cellAnn)
colnames(cellAnnOut)=c("id", "labels")
write.table(cellAnnOut, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#细胞注释后的可视化
newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pdf(file="04.cellAnn.pdf", width=8, height=6)
#DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)     #UMAP可视化
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)           #TSNE可视化
dev.off()
#分组可视化
Type=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
names(Type)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=Type, col.name="Type")
pdf(file="04.group.cellAnn.pdf", width=12, height=6)
#DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE, split.by="Type")       #UMAP可视化
TSNEPlot(object = pbmc, pt.size = 1, label = TRUE, split.by="Type")     #TSNE可视化
dev.off()

#细胞类型的差异分析
pbmc.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)

#细胞组间差异分析
groups=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
groups=paste0(groups, "_", cellAnn)
names(groups)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
for(cellName in unique(cellAnn)){
	conName=paste0("Control_", cellName)
	treatName=paste0("Disease_", cellName)
	if( (length(groups[groups==conName])>5) & (length(groups[groups==treatName])>5) ){
		pbmc.markers=FindMarkers(pbmc, ident.1=treatName, ident.2=conName, group.by='group', logfc.threshold=0.1)
		sig.markersGroup=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
		sig.markersGroup=cbind(Gene=row.names(sig.markersGroup), sig.markersGroup)
		write.table(sig.markersGroup,file=paste0("05.", cellName, ".diffGene.txt"),sep="\t",row.names=F,quote=F)
	}
}

#保存单细胞数据的对象
save(pbmc, cellAnn, clusterAnn, sig.markers, file="Seurat.Rdata")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


