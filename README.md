# MD Trajectory Analysis Toolkit

**AMBER / GROMACS · MDTraj · six-panel summary in one command**

一条命令完成分子动力学轨迹的常规分析：骨架 RMSD、逐残基 RMSF、回转半径、PCA 自由能景观、氢键统计，输出一张六合一图和一张逐帧指标表。支持 AMBER 与 GROMACS 轨迹。

## 功能

| 模块 | 内容 |
|---|---|
| RMSD | 蛋白骨架相对第一帧的均方根偏差 |
| RMSF | 逐残基均方根波动（先逐原子计算，再按残基取平均） |
| Rg | 回转半径 |
| PCA / FEL | 叠合后取前两个主成分，用核密度估计构造自由能景观 |
| H-bond | 蛋白内部氢键；指定 `--ligand` 时只统计蛋白–配体之间的氢键 |
| 指标表 | 导出的 CSV 含 time / RMSD / Rg / PC1 / PC2 |

## 安装

```bash
pip install -r requirements.txt
```

或安装为命令行工具：

```bash
pip install -e .
md-analyze --top 1.prmtop --traj md.nc --out results
```

## 使用

```bash
# AMBER
python src/md_analysis_tool.py --top 1.prmtop --traj md.nc --out results

# GROMACS
python src/md_analysis_tool.py --top topol.tpr --traj traj.xtc --out results

# 带配体，统计蛋白-配体氢键
python src/md_analysis_tool.py --top 1.prmtop --traj md.nc \
    --ligand "resname LIG" --temperature 310 --label "GLB-CYP2C9"
```

安装为命令后，把 `python src/md_analysis_tool.py` 换成 `md-analyze` 即可。

| 参数 | 默认 | 说明 |
|---|---|---|
| `--top` | 必填 | 拓扑文件：`.prmtop` / `.tpr` / `.gro` / `.pdb` |
| `--traj` | 必填 | 轨迹文件：`.nc` / `.xtc` / `.dcd` |
| `--out` | `../results` | 输出目录 |
| `--stride` | `50` | 载入时抽帧步长 |
| `--hb-stride` | `4` | 氢键统计的额外抽帧步长 |
| `--temperature` | `300` | 自由能换算温度（K） |
| `--ligand` | 无 | MDTraj 选择语法，如 `resname LIG` |
| `--label` | 空 | 图标题 |

## 输出

- `results/md_analysis_summary.png`：六合一图（RMSD / RMSF / Rg / FEL / PCA 投影 / 氢键分布）
- `results/md_analysis_metrics.csv`：逐帧指标（`time_ns, rmsd_nm, rg_nm, pc1, pc2`）

## 方法说明

- RMSD 与 RMSF 均基于蛋白骨架原子，参考帧为第一帧；RMSF **按残基取平均**，不是逐原子曲线。
- PCA 之前先对每一帧做最小二乘叠合，去掉整体平移与旋转，避免把全局转动混进主成分。
- 自由能景观取 `ΔG = −RT ln(P/Pmax)`，默认 300 K，上限裁剪到 25 kJ/mol，便于比较不同体系。
- 氢键用 Wernet–Nilsson 判据，一次遍历整条轨迹；未指定 `--ligand` 时统计蛋白内部氢键。
- 距离单位沿用 MDTraj 默认值 nm。

## 测试

```bash
python -m unittest discover -s tests -v
```

9 个用例覆盖逐残基 RMSF 聚合、PCA 叠合、自由能范围、氢键统计与输入校验；CI 会在每次 push 时跑一遍。

## English

MD trajectory analysis in one command: backbone RMSD, per-residue RMSF, radius of gyration, PCA-based free-energy landscape and hydrogen-bond statistics for AMBER or GROMACS trajectories. Writes a six-panel figure and a per-frame metrics CSV. Ligand–protein hydrogen bonds are counted when `--ligand` is given. Units are nanometres (MDTraj default).

## License

MIT — see [LICENSE](LICENSE).
