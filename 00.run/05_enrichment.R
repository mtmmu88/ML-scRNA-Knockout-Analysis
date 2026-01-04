#===============================================================================
# 步骤5：富集分析
# GO/KEGG/Reactome分析验证铁死亡、线粒体通路
#===============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

#======================== 参数设置 ========================
target_gene <- "MDH2"
workDir <- "./05_enrichment"
knockoutDir <- "./04_knockout"

# 富集分析参数
pvalueCutoff <- 0.05
qvalueCutoff <- 0.2

#======================== 创建目录 ========================
if (!dir.exists(workDir)) {
  dir.create(workDir, recursive = TRUE)
}
setwd(workDir)

#======================== 加载敲除分析结果 ========================
cat("正在加载敲除分析结果...\n")
load(file.path(knockoutDir, paste0(target_gene, "_knockout_results.RData")))

# 获取显著变化的基因
genes <- sigDiff$gene
cat("显著变化基因数量：", length(genes), "\n")

# 如果是小鼠基因，转换为人类基因
# 检测是否为小鼠基因（首字母大写，其余小写）
if (any(grepl("^[A-Z][a-z]", genes))) {
  cat("检测到可能是小鼠基因，尝试转换为人类基因符号...\n")
  genes <- toupper(genes)
}

#======================== 基因ID转换 ========================
cat("\n正在进行基因ID转换...\n")
gene_ids <- bitr(genes,
                 fromType = "SYMBOL",
                 toType = c("ENTREZID", "ENSEMBL"),
                 OrgDb = org.Hs.eg.db)

cat("成功转换基因数：", nrow(gene_ids), "/", length(genes), "\n")
entrez_ids <- gene_ids$ENTREZID

#======================== GO富集分析 ========================
cat("\n========== GO富集分析 ==========\n")

# Biological Process
go_BP <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = pvalueCutoff,
  qvalueCutoff = qvalueCutoff,
  readable = TRUE
)

# Cellular Component
go_CC <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = pvalueCutoff,
  qvalueCutoff = qvalueCutoff,
  readable = TRUE
)

# Molecular Function
go_MF <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = pvalueCutoff,
  qvalueCutoff = qvalueCutoff,
  readable = TRUE
)

# 保存GO结果
if (!is.null(go_BP) && nrow(go_BP) > 0) {
  write.table(as.data.frame(go_BP), file = "GO_BP.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("GO BP显著通路：", nrow(go_BP), "\n")
}

if (!is.null(go_CC) && nrow(go_CC) > 0) {
  write.table(as.data.frame(go_CC), file = "GO_CC.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("GO CC显著通路：", nrow(go_CC), "\n")
}

if (!is.null(go_MF) && nrow(go_MF) > 0) {
  write.table(as.data.frame(go_MF), file = "GO_MF.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("GO MF显著通路：", nrow(go_MF), "\n")
}

#======================== KEGG富集分析 ========================
cat("\n========== KEGG富集分析 ==========\n")

kegg <- enrichKEGG(
  gene = entrez_ids,
  organism = 'hsa',
  pAdjustMethod = "BH",
  pvalueCutoff = pvalueCutoff,
  qvalueCutoff = qvalueCutoff
)

if (!is.null(kegg) && nrow(kegg) > 0) {
  # 转换基因ID为符号
  kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.table(as.data.frame(kegg_readable), file = "KEGG.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("KEGG显著通路：", nrow(kegg), "\n")
}

#======================== 检查关键通路 ========================
cat("\n========== 关键通路检查 ==========\n")

# 铁死亡相关关键词
ferroptosis_keywords <- c("ferroptosis", "iron", "lipid peroxid", "glutathione",
                          "oxidative stress", "reactive oxygen", "mitochondri",
                          "electron transport", "TCA cycle", "citrate")

# 检查GO BP
if (!is.null(go_BP) && nrow(go_BP) > 0) {
  go_BP_df <- as.data.frame(go_BP)
  for (kw in ferroptosis_keywords) {
    matches <- go_BP_df[grepl(kw, go_BP_df$Description, ignore.case = TRUE), ]
    if (nrow(matches) > 0) {
      cat("\n找到关键通路 (", kw, "):\n")
      print(matches[, c("Description", "pvalue", "Count")])
    }
  }
}

# 检查KEGG
if (!is.null(kegg) && nrow(kegg) > 0) {
  kegg_df <- as.data.frame(kegg)
  cat("\n\nKEGG通路列表：\n")
  print(kegg_df[, c("Description", "pvalue", "Count")])
}

#======================== 绑制GO气泡图 ========================
cat("\n正在绑制富集分析图...\n")

if (!is.null(go_BP) && nrow(go_BP) > 0) {
  p1 <- dotplot(go_BP, showCategory = 15, title = "GO Biological Process")
  pdf("GO_BP_dotplot.pdf", width = 10, height = 8)
  print(p1)
  dev.off()

  # 条形图
  p2 <- barplot(go_BP, showCategory = 15, title = "GO Biological Process")
  pdf("GO_BP_barplot.pdf", width = 10, height = 8)
  print(p2)
  dev.off()
}

if (!is.null(go_CC) && nrow(go_CC) > 0) {
  p3 <- dotplot(go_CC, showCategory = 15, title = "GO Cellular Component")
  pdf("GO_CC_dotplot.pdf", width = 10, height = 8)
  print(p3)
  dev.off()
}

if (!is.null(go_MF) && nrow(go_MF) > 0) {
  p4 <- dotplot(go_MF, showCategory = 15, title = "GO Molecular Function")
  pdf("GO_MF_dotplot.pdf", width = 10, height = 8)
  print(p4)
  dev.off()
}

#======================== 绘制KEGG图 ========================
if (!is.null(kegg) && nrow(kegg) > 0) {
  p5 <- dotplot(kegg, showCategory = 15, title = "KEGG Pathway")
  pdf("KEGG_dotplot.pdf", width = 10, height = 8)
  print(p5)
  dev.off()

  p6 <- barplot(kegg, showCategory = 15, title = "KEGG Pathway")
  pdf("KEGG_barplot.pdf", width = 10, height = 8)
  print(p6)
  dev.off()
}

#======================== 合并GO结果绘图 ========================
# 合并三个GO本体的结果
go_all <- list(BP = go_BP, CC = go_CC, MF = go_MF)
go_all <- go_all[sapply(go_all, function(x) !is.null(x) && nrow(x) > 0)]

if (length(go_all) > 0) {
  pdf("GO_combined.pdf", width = 12, height = 10)
  for (ont in names(go_all)) {
    p <- dotplot(go_all[[ont]], showCategory = 10, title = paste("GO", ont))
    print(p)
  }
  dev.off()
}

#======================== 保存结果 ========================
save(go_BP, go_CC, go_MF, kegg, gene_ids,
     file = "enrichment_results.RData")

cat("\n========== 富集分析完成 ==========\n")
cat("输出文件：\n")
cat("  - GO_BP.txt / GO_CC.txt / GO_MF.txt（GO结果）\n")
cat("  - KEGG.txt（KEGG结果）\n")
cat("  - GO_BP_dotplot.pdf（GO气泡图）\n")
cat("  - GO_BP_barplot.pdf（GO条形图）\n")
cat("  - KEGG_dotplot.pdf（KEGG气泡图）\n")
cat("  - enrichment_results.RData（R数据）\n")

cat("\n请检查结果中是否包含：\n")
cat("  - 铁死亡相关通路（ferroptosis, iron homeostasis）\n")
cat("  - 线粒体相关通路（mitochondria, electron transport chain）\n")
cat("  - 氧化应激通路（oxidative stress, ROS）\n")
cat("  - 代谢通路（TCA cycle, lipid metabolism）\n")
