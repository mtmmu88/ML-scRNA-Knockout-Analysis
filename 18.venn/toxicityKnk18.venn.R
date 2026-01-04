#install.packages("ggvenn")


#引用包
library(ggvenn)

compoundName="DEP"            #化合物名称
diseaseName="Alzheimer"      #疾病名称
setwd("C:\\Users\\Administrator\\Desktop\\toxicityKnk\\18.venn")    #设置工作目录
geneList=list()

#读取化合物靶基因的文件
rt=read.table("Compound.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])      #提取化合物的靶基因
geneList[[compoundName]]=geneNames

#读取差异分析的结果文件
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])      #提取疾病的差异基因
geneList[[diseaseName]]=geneNames

#绘制venn图
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
	stroke_color = "white", stroke_size = 0.5,
	fill_color = c("#E41A1C","#1E90FF"),
	set_name_color =c("#E41A1C","#1E90FF"),
	set_name_size=6, text_size=4.5)
dev.off()

#输出交集基因
interGenes=Reduce(intersect, geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)

#输出网络关系文件
networkTab=rbind(cbind(compoundName,interGenes, "Compound"), cbind(diseaseName,interGenes, "Disease"))
colnames(networkTab)=c("Node1", "Node2", "Type")
write.table(file="net.network.txt", networkTab, sep="\t", quote=F, row.names=F)
#输出节点属性文件
nodeTab=rbind(cbind(compoundName,"Compound"), cbind(diseaseName,"Disease"), cbind(interGenes,"Gene"))
colnames(nodeTab)=c("Node", "Type")
write.table(file="net.node.txt", nodeTab, sep="\t", quote=F, row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: seqBio


