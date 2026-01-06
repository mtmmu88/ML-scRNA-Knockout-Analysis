# 胸腺瘤相关重症肌无力(MG)研究项目

## 研究目标

通过机器学习筛选胸腺瘤患者发展为MG的预测标志物，并用单细胞虚拟敲除验证关键基因功能。

## 数据来源

| 数据集 | 类型 | 样本量 | 用途 |
|--------|------|--------|------|
| **TCGA-THYM** | Bulk RNA-seq | 116例 (34 MG+ / 82 MG-) | 差异分析 + 114种ML |
| **GSE233180** | 单细胞 scRNA-seq | 12例 EOMG | 基因定位 + 虚拟敲除 |

## 分析流程

```
Step 1: 下载TCGA数据
    ↓
Step 1b: 补充MG分组信息
    ↓
Step 2: limma差异分析 (MG+ vs MG-)
    ↓
Step 3: 114种机器学习筛选关键基因
    ↓
Step 4: 下载单细胞数据 (GSE233180)
    ↓
Step 5: Seurat单细胞分析 + 基因定位
    ↓
Step 6: scTenifoldKnk虚拟敲除
    ↓
Step 7: GO/KEGG富集分析
```

## 运行顺序

```r
# 第一阶段：TCGA Bulk分析
source("step01_download_TCGA.R")      # 下载TCGA-THYM数据
source("step01b_add_MG_status.R")     # 补充MG分组（可能需要手动）
source("step02_diff_analysis.R")      # 差异分析
source("step03_machine_learning.R")   # 114种ML筛选

# 第二阶段：单细胞验证
source("step04_download_scRNA.R")     # 下载GSE233180
source("step05_Seurat_analysis.R")    # 单细胞处理
source("step06_scTenifoldKnk.R")      # 虚拟敲除

# 第三阶段：富集分析
source("step07_enrichment.R")         # GO/KEGG富集
```

## 重要说明

### TCGA数据的MG分组

TCGA-THYM的MG状态可能不在标准下载中，需要从以下来源获取：

1. **cBioPortal**：https://www.cbioportal.org/study/clinicalData?id=thym_tcga_pan_can_atlas_2018
2. **GDC Portal**：https://portal.gdc.cancer.gov/projects/TCGA-THYM
3. **文献Supplementary**：参考TCGA胸腺瘤发表的论文

### 预期结果

| 阶段 | 产出 |
|------|------|
| 差异分析 | 火山图、热图、DEGs列表 |
| 机器学习 | AUC热图、ROC曲线、关键基因 |
| 单细胞 | UMAP、基因表达定位图 |
| 虚拟敲除 | 敲除火山图、受影响基因 |
| 富集分析 | GO/KEGG气泡图 |

## 目录结构

```
MG_Thymoma/
├── step01_download_TCGA.R        # TCGA数据下载
├── step01b_add_MG_status.R       # MG分组补充
├── step02_diff_analysis.R        # 差异分析
├── step03_machine_learning.R     # 机器学习
├── step04_download_scRNA.R       # 单细胞下载
├── step05_Seurat_analysis.R      # Seurat分析
├── step06_scTenifoldKnk.R        # 虚拟敲除
├── step07_enrichment.R           # 富集分析
├── MG_genes.txt                  # MG相关基因列表
├── README.md                     # 本文件
│
├── 01_TCGA_data/                 # TCGA数据
├── 02_diff/                      # 差异分析结果
├── 03_ML/                        # 机器学习结果
├── 04_scRNA/                     # 单细胞数据
├── 05_knockout/                  # 敲除分析结果
└── 06_enrichment/                # 富集分析结果
```

## 参考文献

1. Radovich M, et al. The Integrated Genomic Landscape of Thymic Epithelial Tumors. Cancer Cell. 2018
2. Marx A, et al. Myasthenia gravis-specific aberrant neuromuscular gene expression by medullary thymic epithelial cells in thymoma. Nature Communications. 2022
