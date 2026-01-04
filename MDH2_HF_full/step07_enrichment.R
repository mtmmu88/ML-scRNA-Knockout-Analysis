#===============================================================================
# MDH2心衰研究 - 完整版
# 步骤7：富集分析（GO/KEGG）
#
# 参考原代码：27.GO、28.KEGG、29.Reactome
#===============================================================================

#加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)    #小鼠
library(enrichplot)
library(ggplot2)
library(DOSE)

#======================== 参数设置 ========================
target_gene <- "Mdh2"
pvalueCutoff <- 0.05
qvalueCutoff <- 0.2

#设置工作目录
setwd("D:/MDH2_HF_full")       #Windows用户
#setwd("~/MDH2_HF_full")       #Mac/Linux用户

#======================== 加载敲除结果 ========================
cat("正在加载敲除分析结果...\n")

knockoutFile <- paste0("knockout_", target_gene, ".RData")
if (file.exists(knockoutFile)) {
  load(knockoutFile)
  genes <- sigDiff$gene
  species <- "mouse"
  orgDb <- org.Mm.eg.db
  keggOrg <- "mmu"
  cat("使用单细胞敲除结果（小鼠）\n")
} else {
  #如果没有敲除结果，使用差异分析结果
  cat("未找到敲除结果，使用Bulk差异分析结果\n")
  load("diff.RData")
  genes <- rownames(diffSig)
  species <- "human"
  orgDb <- org.Hs.eg.db
  keggOrg <- "hsa"
}

cat("输入基因数：", length(genes), "\n")

#======================== 基因ID转换 ========================
cat("\n正在进行基因ID转换...\n")

#小鼠基因转人类（如果是小鼠）
if (species == "mouse") {
  #尝试直接转换
  genes_upper <- toupper(genes)
  gene_ids <- bitr(genes_upper,
                   fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  orgDb <- org.Hs.eg.db
  keggOrg <- "hsa"
} else {
  gene_ids <- bitr(genes,
                   fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = orgDb)
}

cat("成功转换：", nrow(gene_ids), "/", length(genes), "\n")
entrez_ids <- gene_ids$ENTREZID

#======================== GO富集分析 ========================
cat("\n========== GO富集分析 ==========\n")

#Biological Process
cat("GO BP...\n")
go_BP <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = pvalueCutoff,
  qvalueCutoff = qvalueCutoff,
  readable = TRUE
)

#Cellular Component
cat("GO CC...\n")
go_CC <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = pvalueCutoff,
  qvalueCutoff = qvalueCutoff,
  readable = TRUE
)

#Molecular Function
cat("GO MF...\n")
go_MF <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = pvalueCutoff,
  qvalueCutoff = qvalueCutoff,
  readable = TRUE
)

#保存GO结果
if (!is.null(go_BP) && nrow(go_BP) > 0) {
  write.table(as.data.frame(go_BP), file = "bindao_enrichment_GO_BP.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("GO BP显著通路：", nrow(go_BP), "\n")
}

if (!is.null(go_CC) && nrow(go_CC) > 0) {
  write.table(as.data.frame(go_CC), file = "bindao_enrichment_GO_CC.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("GO CC显著通路：", nrow(go_CC), "\n")
}

if (!is.null(go_MF) && nrow(go_MF) > 0) {
  write.table(as.data.frame(go_MF), file = "bindao_enrichment_GO_MF.txt",
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
  kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.table(as.data.frame(kegg_readable), file = "bindao_enrichment_KEGG.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("KEGG显著通路：", nrow(kegg), "\n")
}

#======================== 检查关键通路 ========================
cat("\n========== 关键通路检查 ==========\n")

keywords <- c("ferroptosis", "iron", "mitochondri", "oxidative",
              "lipid peroxid", "glutathione", "TCA", "electron transport")

if (!is.null(go_BP) && nrow(go_BP) > 0) {
  go_BP_df <- as.data.frame(go_BP)
  cat("\nGO BP中的关键通路：\n")
  for (kw in keywords) {
    matches <- go_BP_df[grepl(kw, go_BP_df$Description, ignore.case = TRUE), ]
    if (nrow(matches) > 0) {
      cat(paste0("\n[", kw, "]\n"))
      print(matches[, c("Description", "pvalue", "Count")])
    }
  }
}

if (!is.null(kegg) && nrow(kegg) > 0) {
  kegg_df <- as.data.frame(kegg)
  cat("\n\nKEGG中的关键通路：\n")
  for (kw in keywords) {
    matches <- kegg_df[grepl(kw, kegg_df$Description, ignore.case = TRUE), ]
    if (nrow(matches) > 0) {
      cat(paste0("\n[", kw, "]\n"))
      print(matches[, c("Description", "pvalue", "Count")])
    }
  }
}

#======================== 绑制图形 ========================
cat("\n正在绘制富集分析图...\n")

#GO BP气泡图
if (!is.null(go_BP) && nrow(go_BP) > 0) {
  p1 <- dotplot(go_BP, showCategory = 15, title = "GO Biological Process")
  pdf("bindao_bindao_enrichment_GO_BP_bindao_bindao_bindao_dotplot.pdf", width = 10, height = 8)
  print(p1)
  dev.off()

  p2 <- barplot(go_BP, showCategory = 15, title = "GO Biological Process")
  pdf("bindao_enrichment_GO_BP_bindao_barplot.pdf", width = 10, height = 8)
  print(p2)
  dev.off()
}

#GO CC气泡图
if (!is.null(go_CC) && nrow(go_CC) > 0) {
  p3 <- dotplot(go_CC, showCategory = 15, title = "GO Cellular Component")
  pdf("bindao_enrichment_GO_CC_bindao_bindao_bindao_dotplot.pdf", width = 10, height = 8)
  print(p3)
  dev.off()
}

#GO MF气泡图
if (!is.null(go_MF) && nrow(go_MF) > 0) {
  p4 <- dotplot(go_MF, showCategory = 15, title = "GO Molecular Function")
  pdf("bindao_enrichment_GO_MF_bindao_bindao_bindao_dotplot.pdf", width = 10, height = 8)
  print(p4)
  dev.off()
}

#KEGG气泡图
if (!is.null(kegg) && nrow(kegg) > 0) {
  p5 <- dotplot(kegg, showCategory = 15, title = "KEGG Pathway")
  pdf("bindao_enrichment_KEGG_bindao_bindao_bindao_dotplot.pdf", width = 10, height = 8)
  print(p5)
  dev.off()

  p6 <- barplot(kegg, showCategory = 15, title = "KEGG Pathway")
  pdf("bindao_enrichment_KEGG_bindao_barplot.pdf", width = 10, height = 8)
  print(p6)
  dev.off()
}

#======================== 保存结果 ========================
save(go_BP, go_CC, go_MF, kegg, gene_ids,
     file = "bindao_enrichment_results.RData")

cat("\n========== 富集分析完成 ==========\n")
cat("输出文件：\n")
cat("  - bindao_enrichment_GO_BP/CC/MF.txt\n")
cat("  - bindao_enrichment_KEGG.txt\n")
cat("  - bindao_enrichment_*_bindao_bindao_bindao_dotplot.pdf\n")
cat("  - bindao_enrichment_*_bindao_barplot.pdf\n")

cat("\n========== 全部分析完成！ ==========\n")
cat("\n请检查结果中是否包含：\n")
cat("  ✓ 铁死亡通路（Ferroptosis）\n")
cat("  ✓ 线粒体代谢（mitochondria, TCA cycle）\n")
cat("  ✓ 氧化应激（oxidative stress, ROS）\n")
cat("  ✓ 脂质代谢（lipid metabolism）\n")


######生信自学网: https://www.biowolf.cn/
######MDH2心衰研究完整版
