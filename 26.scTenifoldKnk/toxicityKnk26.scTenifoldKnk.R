# install.packages("remotes")
# remotes::install_github('cailab-tamu/scTenifoldKnk')
# install.packages("Seurat")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("ggrepel")


# 加载包
library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggrepel)
set.seed(123)


target_gene="HSPA1A"     #目标基因的名称(需要修改)
setwd("C:\\Users\\Administrator\\Desktop\\toxicityKnk\\26.scTenifoldKnk")      #设置工作目录


# 加载Seurat的对象
load("seurat.Rdata")
# 只保留疾病组的样品
pbmc <- subset(pbmc, subset = Type == "Disease")

#提取表达矩阵
countMat = GetAssayData(pbmc, layer = "counts")
#提取高可变基因
pbmc <- FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=5000)
hvgs <- VariableFeatures(pbmc)
data=as.data.frame(countMat[unique(c(target_gene,hvgs)),])

#选取表达均值最高的前20000个细胞进行分析
#cell_mean <- colMeans(data, na.rm = TRUE)
#top_cell_idx <- order(cell_mean, decreasing = TRUE)[1:min(20000, length(cell_mean))]
#data <- data[, top_cell_idx]
rm(pbmc, countMat);gc()

#执行虚拟敲除
result <- scTenifoldKnk(countMatrix = data, 
                        gKO = target_gene,          #需要敲除的基因
                        qc_mtThreshold = 0.1,       #mt的阈值
                        qc_minLSize = 1000,         #文库阈值(细胞测到的基因总数)
                        nc_nNet = 10,                #子网络数量
                        nc_nCells = 500              #每个网络中随机抽取的细胞数
                        )

#输出差异分析的结果
df=result$diffRegulation
df <- df[df$gene != target_gene, ]
outTab=df[df$p.adj<0.05,]
write.table(outTab, file="sigDiff.txt", sep="\t", quote=F, row.names=F)


###########################柱状图###########################
#绘制柱状图
top_genes <- head(df[order(-df$FC), ], 20)
p1=ggplot(top_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='#5A9BD4') +
  coord_flip() + 
  labs(title="Top 20 Differentially Regulated Genes", x="Gene", y="FC") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
#输出图形
pdf(file="barplot.pdf", width=6, height=5)
print(p1)
dev.off()


###########################火山图###########################
#准备火山图的数据
df$log_p.adj <- -log10(df$p.adj)
df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")
label_genes <- subset(df, p.adj < 0.05)

# 绘制火山图
y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)
p2=ggplot(df, aes(x=Z, y=log_p.adj, color=significant)) +
  geom_point(alpha=0.7, size=1) +  # 设置点的透明度和大小
  # 定义图形的颜色
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  #在图形中标注显著差异基因的名称
  geom_text_repel(data=label_genes, aes(label=gene), size=3, max.overlaps=50) +
  labs(title="", x="Z-score", y="-log10(p.adj)") +
  theme_classic() + 
  # 设置纵坐标范围
  coord_cartesian(ylim = c(0, y_upper)) + theme(legend.position = "none")
#输出图形
pdf(file="vol.pdf", width=6, height=5)
print(p2)
dev.off()


###########################饼图###########################
#准备饼图数据
sig_count <- table(df$significant)
sig_df <- as.data.frame(sig_count)
colnames(sig_df) <- c("category", "count")
sig_df$percentage <- paste0(round(sig_df$count / sum(sig_df$count) * 100, 1), "%")
# 绘制饼图
p3 <- ggplot(sig_df, aes(x="", y=count, fill=category)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label=percentage), position=position_stack(vjust=0.5), size=4) + 
  coord_polar("y", start=0) + 
  scale_fill_manual(values=c("Significant"="red", "Not significant"="lightgray")) + 
  labs(title="Proportion of Significant Genes", fill="") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_blank(),
    panel.grid = element_blank(), legend.position = "right", plot.title = element_text(hjust = 0.5))
# 输出饼图
pdf(file="pie.pdf", width=6, height=5)
print(p3)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


