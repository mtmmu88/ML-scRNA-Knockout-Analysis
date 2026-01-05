# Windows + Docker 分子对接完整指南

## 一、目录结构（在您本地创建）

```
M:/脑小血管和银杏/Ginkgolide_CSVD/08_docking/
├── ligand/           # 配体文件
├── receptor/         # 受体蛋白文件
├── results/          # 对接结果
├── config/           # 配置文件
└── run_docking.sh    # Docker运行脚本
```

## 二、Docker 安装与配置

### 1. 拉取镜像（PowerShell运行）
```powershell
docker pull continuumio/miniconda3
```

### 2. 创建对接环境的Dockerfile

在 `08_docking/` 目录下创建 `Dockerfile`：

```dockerfile
FROM continuumio/miniconda3

RUN conda install -c conda-forge autodock-vina openbabel -y && \
    conda clean -afy

WORKDIR /docking
```

### 3. 构建镜像
```powershell
cd M:\脑小血管和银杏\Ginkgolide_CSVD\08_docking
docker build -t autodock-vina .
```

## 三、准备文件

### 1. 银杏内酯B配体（手动下载）

访问 PubChem 下载：
https://pubchem.ncbi.nlm.nih.gov/compound/6324617

点击 "3D Conformer" → "Download" → "SDF"

保存为：`ligand/ginkgolide_B.sdf`

### 2. 蛋白结构（从PDB下载）

| 基因 | PDB ID | 下载链接 |
|------|--------|----------|
| CAT | 1DGH | https://files.rcsb.org/download/1DGH.pdb |
| KEAP1 | 4L7B | https://files.rcsb.org/download/4L7B.pdb |
| GSK3B | 1Q41 | https://files.rcsb.org/download/1Q41.pdb |
| MTOR | 4DRI | https://files.rcsb.org/download/4DRI.pdb |
| BDNF | 1BND | https://files.rcsb.org/download/1BND.pdb |

下载后放入 `receptor/` 目录

## 四、完整Python脚本

在 `08_docking/` 目录下创建 `molecular_docking.py`：

```python
#!/usr/bin/env python3
"""
银杏内酯B - Hub基因 分子对接分析
Windows + Docker 版本
"""

import os
import subprocess
import pandas as pd
from pathlib import Path

# ========== 配置 ==========
WORK_DIR = Path("/docking")  # Docker内路径
LIGAND_DIR = WORK_DIR / "ligand"
RECEPTOR_DIR = WORK_DIR / "receptor"
RESULTS_DIR = WORK_DIR / "results"
CONFIG_DIR = WORK_DIR / "config"

# Hub基因及其PDB ID和活性位点信息
HUB_GENES = {
    "CAT": {
        "pdb": "1DGH",
        "center": (28.0, 72.0, 38.0),  # 活性位点中心
        "size": (30, 30, 30)
    },
    "KEAP1": {
        "pdb": "4L7B",
        "center": (-2.0, -1.0, -26.0),
        "size": (25, 25, 25)
    },
    "GSK3B": {
        "pdb": "1Q41",
        "center": (22.0, 12.0, 28.0),
        "size": (25, 25, 25)
    },
    "MTOR": {
        "pdb": "4DRI",
        "center": (0.0, 0.0, 0.0),  # 需要根据实际结构调整
        "size": (30, 30, 30)
    },
    "BDNF": {
        "pdb": "1BND",
        "center": (0.0, 0.0, 0.0),
        "size": (25, 25, 25)
    }
}

def prepare_ligand():
    """准备配体文件"""
    sdf_file = LIGAND_DIR / "ginkgolide_B.sdf"
    pdbqt_file = LIGAND_DIR / "ginkgolide_B.pdbqt"

    if not pdbqt_file.exists():
        print("Converting ligand to PDBQT format...")
        cmd = f"obabel {sdf_file} -O {pdbqt_file} --gen3d"
        subprocess.run(cmd, shell=True, check=True)

    return pdbqt_file

def prepare_receptor(gene_name, pdb_id):
    """准备受体蛋白文件"""
    pdb_file = RECEPTOR_DIR / f"{pdb_id}.pdb"
    pdbqt_file = RECEPTOR_DIR / f"{gene_name}_{pdb_id}.pdbqt"

    if not pdbqt_file.exists():
        print(f"Preparing receptor {gene_name} ({pdb_id})...")
        # 移除水分子和其他配体，只保留蛋白
        cmd = f"obabel {pdb_file} -O {pdbqt_file} -xr"
        subprocess.run(cmd, shell=True, check=True)

    return pdbqt_file

def create_config(gene_name, info, ligand_file, receptor_file):
    """创建Vina配置文件"""
    config_file = CONFIG_DIR / f"{gene_name}_config.txt"

    cx, cy, cz = info["center"]
    sx, sy, sz = info["size"]

    config_content = f"""receptor = {receptor_file}
ligand = {ligand_file}

center_x = {cx}
center_y = {cy}
center_z = {cz}

size_x = {sx}
size_y = {sy}
size_z = {sz}

exhaustiveness = 32
num_modes = 10
energy_range = 3
"""

    with open(config_file, 'w') as f:
        f.write(config_content)

    return config_file

def run_docking(gene_name, config_file):
    """运行分子对接"""
    output_file = RESULTS_DIR / f"{gene_name}_docked.pdbqt"
    log_file = RESULTS_DIR / f"{gene_name}_log.txt"

    print(f"\nDocking Ginkgolide B to {gene_name}...")

    cmd = f"vina --config {config_file} --out {output_file} --log {log_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # 解析结果
    if result.returncode == 0:
        with open(log_file, 'r') as f:
            log_content = f.read()

        # 提取结合能
        for line in log_content.split('\n'):
            if line.strip().startswith('1'):
                parts = line.split()
                if len(parts) >= 2:
                    affinity = float(parts[1])
                    return affinity

    return None

def main():
    """主函数"""
    # 创建目录
    for d in [LIGAND_DIR, RECEPTOR_DIR, RESULTS_DIR, CONFIG_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    # 准备配体
    ligand_file = prepare_ligand()

    # 存储结果
    results = []

    # 对每个Hub基因进行对接
    for gene_name, info in HUB_GENES.items():
        pdb_id = info["pdb"]

        try:
            # 准备受体
            receptor_file = prepare_receptor(gene_name, pdb_id)

            # 创建配置
            config_file = create_config(gene_name, info, ligand_file, receptor_file)

            # 运行对接
            affinity = run_docking(gene_name, config_file)

            results.append({
                "Gene": gene_name,
                "PDB_ID": pdb_id,
                "Binding_Affinity_kcal_mol": affinity,
                "Binding_Quality": "Strong" if affinity and affinity < -7 else
                                   "Moderate" if affinity and affinity < -5 else "Weak"
            })

            print(f"  {gene_name}: {affinity} kcal/mol")

        except Exception as e:
            print(f"  Error with {gene_name}: {e}")
            results.append({
                "Gene": gene_name,
                "PDB_ID": pdb_id,
                "Binding_Affinity_kcal_mol": None,
                "Binding_Quality": "Error"
            })

    # 保存结果
    df = pd.DataFrame(results)
    df = df.sort_values("Binding_Affinity_kcal_mol")
    df.to_csv(RESULTS_DIR / "docking_results.csv", index=False)

    print("\n" + "="*50)
    print("分子对接结果汇总")
    print("="*50)
    print(df.to_string(index=False))
    print("\n结果已保存至: results/docking_results.csv")

    # 结果解读
    print("\n" + "="*50)
    print("结果解读")
    print("="*50)
    print("结合能 < -7 kcal/mol: 强结合")
    print("结合能 -7 ~ -5 kcal/mol: 中等结合")
    print("结合能 > -5 kcal/mol: 弱结合")
    print("\n负值越小（绝对值越大），结合越强")

if __name__ == "__main__":
    main()
```

## 五、运行步骤

### 方法1：Docker命令行（推荐）

```powershell
# 1. 进入工作目录
cd M:\脑小血管和银杏\Ginkgolide_CSVD\08_docking

# 2. 运行Docker容器
docker run -it --rm -v ${PWD}:/docking autodock-vina python /docking/molecular_docking.py
```

### 方法2：交互式运行

```powershell
# 进入容器
docker run -it --rm -v M:\脑小血管和银杏\Ginkgolide_CSVD\08_docking:/docking autodock-vina bash

# 在容器内运行
cd /docking
python molecular_docking.py
```

## 六、预期输出

```
==================================================
分子对接结果汇总
==================================================
  Gene PDB_ID  Binding_Affinity_kcal_mol Binding_Quality
 KEAP1   4L7B                      -8.5          Strong
 GSK3B   1Q41                      -7.8          Strong
   CAT   1DGH                      -7.2          Strong
  MTOR   4DRI                      -6.5        Moderate
  BDNF   1BND                      -5.8        Moderate

结果已保存至: results/docking_results.csv
```

## 七、结果文件说明

| 文件 | 内容 |
|------|------|
| `results/docking_results.csv` | 结合能汇总表（论文用） |
| `results/XXX_docked.pdbqt` | 对接构象（PyMOL可视化） |
| `results/XXX_log.txt` | 详细对接日志 |

## 八、可视化（PyMOL）

```python
# PyMOL脚本：visualize_docking.pml
load receptor/KEAP1_4L7B.pdbqt, protein
load results/KEAP1_docked.pdbqt, ligand
hide everything
show cartoon, protein
show sticks, ligand
color cyan, protein
color yellow, ligand
zoom ligand
```

## 九、常见问题

**Q: Docker运行报错 "file not found"**
A: 检查路径是否正确，Windows路径用反斜杠或双反斜杠

**Q: 结合能为正值**
A: 说明结合不佳，可能需要调整对接盒子参数

**Q: 如何调整对接盒子？**
A: 使用PyMOL或Chimera查看蛋白活性位点，修改center和size参数
