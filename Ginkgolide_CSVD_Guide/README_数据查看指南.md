# 银杏内酯B-脑小血管病网络药理学分析 - 数据查看指南

## 项目结构

```
Ginkgolide_CSVD/
├── 01_data/          # 原始数据
├── 02_diff/          # 差异分析 ★★★
├── 03_compound/      # 化合物靶点
├── 04_venn/          # 交集基因 ★★★
├── 05_enrichment/    # 富集分析 ★★
├── 06_PPI/           # PPI网络 ★★★
└── 07_celltype/      # 细胞类型定位 ★★
```

---

## 核心结果文件（必看）

### 1. 差异表达基因（02_diff/）
| 文件 | 内容 | 用途 |
|------|------|------|
| `DEGs_VaD_vs_Control.csv` | 全部差异基因列表 | 了解VaD疾病特征 |
| `volcano_plot.pdf` | 火山图 | 论文Figure |

**关键列：**
- `logFC`: 表达变化倍数（正=上调，负=下调）
- `adj.P.Val`: 校正后P值（<0.05有意义）

---

### 2. 交集基因（04_venn/） ★最重要
| 文件 | 内容 |
|------|------|
| `intersection_genes.csv` | **9个核心靶点** |
| `venn_diagram.pdf` | Venn图（论文Figure） |

**9个交集基因：**
```
PTAFR, SELP, F3, MTOR, GSK3B, KEAP1, CAT, GPX1, BDNF
```

---

### 3. PPI网络（06_PPI/） ★关键
| 文件 | 内容 | 用途 |
|------|------|------|
| `ppi_network.csv` | 蛋白互作边列表 | Cytoscape可视化 |
| `hub_genes.csv` | **5个Hub基因** | 核心靶点 |
| `node_degree.csv` | 节点度数统计 | 筛选依据 |

**5个Hub基因（按度数排序）：**
| 排名 | 基因 | 度数 | 生物学功能 |
|------|------|------|-----------|
| 1 | CAT | 10 | 过氧化氢酶，清除ROS |
| 2 | KEAP1 | 8 | Nrf2抑制因子，氧化应激调控 |
| 3 | GSK3B | 8 | 糖原合成酶激酶，细胞存活 |
| 4 | MTOR | 8 | 自噬和蛋白合成调控 |
| 5 | BDNF | 6 | 神经营养因子 |

---

### 4. 富集分析（05_enrichment/）
| 文件 | 内容 |
|------|------|
| `GO_BP.csv` | 生物过程富集 |
| `GO_MF.csv` | 分子功能富集 |
| `GO_CC.csv` | 细胞组分富集 |
| `KEGG_pathway.csv` | KEGG通路富集 |
| `enrichment_dotplot.pdf` | 气泡图（论文Figure） |

**主要富集通路：**
- 氧化应激响应（Response to oxidative stress）
- PI3K-Akt信号通路
- 细胞凋亡调控
- 神经营养因子信号

---

### 5. 细胞类型定位（07_celltype/）
| 文件 | 内容 |
|------|------|
| `celltype_expression.csv` | Hub基因在各细胞类型的表达量 |
| `dotplot_celltype.pdf` | 点图（论文Figure） |
| `UMAP_celltype.pdf` | 细胞类型UMAP图 |

**表达特征总结：**
| 基因 | 主要表达细胞 | 意义 |
|------|-------------|------|
| CAT | 星形胶质细胞 | 胶质细胞介导抗氧化 |
| KEAP1 | 星形胶质细胞 | Nrf2通路调控 |
| GSK3B | OPC | 影响髓鞘再生 |
| MTOR | 寡突胶质细胞 | 白质自噬调控 |
| BDNF | 神经元 | 直接神经保护 |

---

## 论文写作用图

| 图号 | 文件 | 内容 |
|------|------|------|
| Fig 1 | `venn_diagram.pdf` | 韦恩图 |
| Fig 2 | `volcano_plot.pdf` | 火山图 |
| Fig 3 | `ppi_network.pdf` | PPI网络（或用Cytoscape重绘） |
| Fig 4 | `enrichment_dotplot.pdf` | 富集分析气泡图 |
| Fig 5 | `dotplot_celltype.pdf` | 细胞类型表达热图 |

---

## 快速结论

**银杏内酯B干预VaD的三条主要机制：**

1. **抗氧化应激**（CAT、KEAP1、GPX1）
   - 靶细胞：星形胶质细胞
   - 机制：清除ROS，激活Nrf2通路

2. **神经保护**（BDNF、MTOR）
   - 靶细胞：神经元、寡突胶质细胞
   - 机制：神经营养、自噬调控

3. **细胞存活/髓鞘保护**（GSK3B、MTOR）
   - 靶细胞：OPC
   - 机制：促进白质修复

---

## 下一步建议

1. **文献验证**：查阅5个Hub基因与VaD/CSVD的已有研究
2. **分子对接**：银杏内酯B与Hub基因蛋白的对接验证
3. **实验验证**：细胞实验验证银杏内酯B对Hub基因的调控
