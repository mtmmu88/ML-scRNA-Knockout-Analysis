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
        "center": (28.0, 72.0, 38.0),
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
        "center": (0.0, 0.0, 0.0),
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
        if sdf_file.exists():
            cmd = f"obabel {sdf_file} -O {pdbqt_file} --gen3d"
        else:
            # 使用SMILES生成
            smiles = "CC1C2(C(C3(C(C4(C(C5(C3(CC(=O)O5)O)C(=O)OC4O)O)OC2=O)C)O)OC1=O)C"
            mol2_file = LIGAND_DIR / "ginkgolide_B.mol2"
            cmd1 = f'echo "{smiles}" | obabel -ismi -O {mol2_file} --gen3d --best'
            subprocess.run(cmd1, shell=True, check=True)
            cmd = f"obabel {mol2_file} -O {pdbqt_file}"
        subprocess.run(cmd, shell=True, check=True)

    return pdbqt_file

def prepare_receptor(gene_name, pdb_id):
    """准备受体蛋白文件"""
    pdb_file = RECEPTOR_DIR / f"{pdb_id}.pdb"
    pdbqt_file = RECEPTOR_DIR / f"{gene_name}_{pdb_id}.pdbqt"

    if not pdb_file.exists():
        print(f"  Warning: {pdb_file} not found. Please download from PDB.")
        return None

    if not pdbqt_file.exists():
        print(f"Preparing receptor {gene_name} ({pdb_id})...")
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

    if result.returncode == 0 and log_file.exists():
        with open(log_file, 'r') as f:
            log_content = f.read()

        for line in log_content.split('\n'):
            if line.strip().startswith('1'):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        affinity = float(parts[1])
                        return affinity
                    except:
                        pass

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
            receptor_file = prepare_receptor(gene_name, pdb_id)

            if receptor_file is None:
                results.append({
                    "Gene": gene_name,
                    "PDB_ID": pdb_id,
                    "Binding_Affinity_kcal_mol": None,
                    "Binding_Quality": "Missing PDB"
                })
                continue

            config_file = create_config(gene_name, info, ligand_file, receptor_file)
            affinity = run_docking(gene_name, config_file)

            if affinity:
                quality = "Strong" if affinity < -7 else "Moderate" if affinity < -5 else "Weak"
            else:
                quality = "Error"

            results.append({
                "Gene": gene_name,
                "PDB_ID": pdb_id,
                "Binding_Affinity_kcal_mol": affinity,
                "Binding_Quality": quality
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

    print("\n" + "="*50)
    print("结果解读")
    print("="*50)
    print("结合能 < -7 kcal/mol: 强结合")
    print("结合能 -7 ~ -5 kcal/mol: 中等结合")
    print("结合能 > -5 kcal/mol: 弱结合")
    print("\n负值越小（绝对值越大），结合越强")

if __name__ == "__main__":
    main()
