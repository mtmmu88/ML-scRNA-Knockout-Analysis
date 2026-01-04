# MDH2心衰研究 - 测试版本

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

## 测试数据

- **GSE26887**：24例心衰样本（5健康 + 19心衰）
- 样本来源：左心室心肌活检
- 平台：Affymetrix Human Gene 1.0 ST

## 运行步骤

```
step00_install.R    →  安装R包（只需运行一次）
        ↓
step01_download.R   →  下载GEO数据
        ↓
step02_normalize.R  →  数据标准化
        ↓
step03_diff.R       →  差异分析（MDH2/SIRT5/铁死亡基因）
        ↓
step04_ML.R         →  机器学习（10种方法测试版）
```

## 运行前准备

1. 修改每个脚本中的 `setwd()` 为你的本地路径
2. 确保 `21.ML/refer.ML.R` 文件路径正确

## 预期输出

### step03_diff.R
- `diff.bindao_all.txt` - 所有基因差异结果
- `diff.bindao_sig.txt` - 显著差异基因
- `diff.bindao_targetGenes.txt` - MDH2/SIRT5/铁死亡基因结果
- `bindao_bindao_volcano.pdf` - 火山图
- `bindao_bindao_heatmap_target.pdf` - 目标基因热图

### step04_ML.R
- `ML.AUCmatrix.txt` - 各模型AUC值
- `ML.AUCheatmap.pdf` - AUC热图
- `ML.MDH2_ROC.pdf` - MDH2单独诊断ROC曲线
- `ML.bestModelGenes.txt` - 最佳模型使用的基因

## 测试通过后

1. 把数据集换成 **GSE59867**（111例心梗后心衰）
2. 把机器学习方法扩展到 **114种**
3. 添加单细胞敲除分析（scTenifoldKnk）

## 注意事项

- step02_normalize.R 中的分组代码可能需要根据实际数据调整
- step04_ML.R 需要 `21.ML/refer.ML.R` 中的函数

## 联系

参考：生信自学网 https://www.biowolf.cn/
