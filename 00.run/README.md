# MDH2-心衰研究：机器学习+单细胞基因敲除分析流程

## 研究假设

```
Sirt5（去乳酸化酶）
    ↓ 调控
MDH2 乳酸化修饰
    ↓ 导致
线粒体代谢紊乱
    ↓ 激活
铁死亡通路
    ↓ 引发
心梗后心衰
```

---

## 运行顺序

```
00_install_packages.R  →  安装R包（只需运行一次）
        ↓
01_download_GEO.R      →  下载GEO数据（GSE59867心衰）
        ↓
02_diff_analysis.R     →  差异分析（MDH2/SIRT5表达变化）
        ↓
03_machine_learning.R  →  机器学习（114种方法，诊断价值）
        ↓
04_scTenifoldKnk.R     →  虚拟敲除MDH2（需要单细胞数据）
        ↓
05_enrichment.R        →  富集分析（验证铁死亡/线粒体通路）
```

---

## 数据集说明

### Bulk数据（用于差异分析+机器学习）
- **GSE59867**：人类心梗后心衰，111例，外周血PBMC

### 单细胞数据（用于敲除分析）
- **GSE132146**：小鼠心梗后心脏单细胞

---

## 输出结果

### 02_diff：差异分析
- `volcano_plot.pdf` - 火山图
- `target_genes_result.txt` - MDH2/SIRT5/铁死亡基因表达

### 03_ML：机器学习
- `model.AUCheatmap.pdf` - 114种模型AUC比较
- `MDH2_ROC.pdf` - MDH2单独诊断价值
- `best_model_genes.txt` - 最佳模型使用的基因

### 04_knockout：敲除分析
- `MDH2_knockout_significant.txt` - 敲除后显著变化的基因
- `MDH2_ferroptosis_genes.txt` - 铁死亡基因变化
- `MDH2_volcano.pdf` - 敲除效应火山图

### 05_enrichment：富集分析
- `GO_BP_dotplot.pdf` - GO富集气泡图
- `KEGG_dotplot.pdf` - KEGG通路图

---

## 快速开始

```r
# 1. 首先安装包
source("00_install_packages.R")

# 2. 下载数据
source("01_download_GEO.R")

# 3. 差异分析（需要根据实际数据修改分组）
source("02_diff_analysis.R")

# 4. 机器学习
source("03_machine_learning.R")

# 5. 敲除分析（需要先准备单细胞数据）
source("04_scTenifoldKnk.R")

# 6. 富集分析
source("05_enrichment.R")
```

---

## 注意事项

1. **分组信息**：运行02_diff_analysis.R时，需要根据实际的phenoData修改分组代码

2. **单细胞数据**：scTenifoldKnk需要单细胞数据，需要另外准备Seurat对象

3. **运行时间**：
   - 差异分析：几分钟
   - 机器学习：10-30分钟
   - scTenifoldKnk：30分钟-数小时

4. **内存需求**：
   - Bulk分析：8GB RAM够用
   - 单细胞分析：建议16GB以上

---

## 铁死亡基因集

详见 `ferroptosis_genes.txt`，包含约50个铁死亡相关基因

---

## 联系方式

如有问题，请检查各脚本中的参数设置
