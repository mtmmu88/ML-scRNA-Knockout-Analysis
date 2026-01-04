#install.packages("rms")
#install.packages("rmda")


#引用包
library(rms)
library(rmda)

inputFile="merge.normalize.txt"       #表达数据文件
geneFile="optimalModelGenes.txt"      #基因列表文件
setwd("C:\\Users\\Administrator\\Desktop\\toxicityKnk\\24.Nomo")      #设置工作目录

#读取表达数据文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#row.names(data)=gsub("-", "_", row.names(data))

#读取基因列表文件, 提取模型基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
paste(rownames(data), collapse="+")     #获取模型基因的列表

#获取样品分组信息(对照组和疾病组)
data=t(data)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
rt=cbind(as.data.frame(data), Type=group)

#数据打包
ddist=datadist(rt)
options(datadist="ddist")

#构建模型，绘制列线图
lrmModel=lrm(Type~ ., data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.001,0.01,0.1,0.5,0.9,0.99),
	lp=F, funlabel="Risk of Disease")
#输出列线图
pdf(file="Nomo.pdf", width=8, height=6)
plot(nomo)
dev.off()

#绘制校准曲线
cali=calibrate(lrmModel, method="boot", B=1000)
pdf(file="Calibration.pdf", width=5.5, height=5.5)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()

#绘制决策曲线
rt$Type=ifelse(rt$Type=="Control", 0, 1)
dc=decision_curve(Type ~ GABRG2+FUCA1+OGFRL1+HSPA1A+SPPL2A+SENP8+CAPN2+PPP1CA+CAPNS1+GABRA4+PPARG+CTDSP1+FOS, data=rt, 
	family = binomial(link ='logit'),
	thresholds= seq(0,1,by = 0.01),
	confidence.intervals = 0.95)
#输出DCA图形
pdf(file="DCA.pdf", width=6, height=6)
plot_decision_curve(dc,
	curve.names="Model",
	xlab="Threshold probability",
	cost.benefit.axis=T,
	col="red",
	confidence.intervals=FALSE,
	standardize=FALSE)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


