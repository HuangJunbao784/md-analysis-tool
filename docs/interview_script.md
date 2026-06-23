# 面试项目介绍 — 黄俊宝

## 中文版 (3分钟)

### 项目1: EGFR抑制剂QSAR模型

> 我独立完成了一个EGFR激酶抑制剂的QSAR建模项目。数据来源是ChEMBL 33数据库，下载了6,896个EGFR抑制剂的生物活性数据。我用RDKit计算了12个分子描述符，包括分子量、logP、TPSA、氢键供体/受体、可旋转键等。然后对比了Random Forest和Gradient Boosting两个模型，用80/20训练测试划分和5折交叉验证评估。
>
> 最终Random Forest在测试集上R²达到0.45，5折交叉验证R²为0.49。特征重要性分析显示TPSA、logP和芳香环数是最重要的三个描述符，和EGFR激酶口袋的理化性质吻合。整个pipeline从数据清洗到出图都是Python自动化完成的，代码在GitHub上开源。

### 项目2: MD轨迹自动化分析工具

> 我还写了一个MD轨迹自动化分析脚本，能处理AMBER和GROMACS的轨迹文件。功能包括RMSD、RMSF、回转半径、PCA本质动力学、自由能景观和配体-蛋白氢键分析。分析结果自动生成论文级图表。
>
> 这个工具替代了我原来手动用GROMACS命令行和VMD可视化的工作流程，一个脚本跑完原来需要半天的手动操作。最近一次分析是500ns的蛋白-配体复合物轨迹，识别出了与配体形成氢键占有率超过50%的关键残基。

### 个人背景一句话

> 我硕士是做药物化学的，研究方向是计算生物学与AI药物化学。我精通分子对接、MD模拟、QM计算，同时自学了Python和机器学习，能把传统CADD方法和AI工具结合起来做完整的计算化学分析管线。

### 为什么想来

> 贵公司在AI药物发现领域有很深的技术积累，我需要一个能把我的计算化学背景和AI技能结合起来发挥的平台。我可以用传统CADD方法理解药物-靶点相互作用的物理化学本质，同时用ML工具加速筛选和预测。

---

## 英文版 (3分钟)

### Project 1: EGFR QSAR Pipeline

> I built a QSAR modeling pipeline for EGFR kinase inhibitors using 6,896 compounds from ChEMBL 33. I extracted 12 molecular descriptors with RDKit and compared Random Forest against Gradient Boosting. The best model achieved test R² of 0.45 and 5-fold CV R² of 0.49. Top features were TPSA, logP, and aromatic ring count, which align with EGFR pocket physicochemical properties. The entire pipeline is automated in Python and open-sourced on GitHub.

### Project 2: MD-Aware QSAR and Trajectory Analysis

> I developed an automated MD trajectory analysis toolkit covering RMSD, RMSF, radius of gyration, PCA-based free energy landscape, and ligand-protein hydrogen bond analysis. Currently, I'm integrating MD-derived features such as per-residue RMSF and H-bond occupancy into the QSAR pipeline to improve prediction beyond traditional 2D descriptors.

### Quick intro

> I'm a master's student in medicinal chemistry at [学校名], specializing in computational biology and AI-driven drug discovery. I'm proficient in molecular docking, MD simulation, and QM calculations, and I've self-taught Python and machine learning to bridge the gap between traditional CADD and modern AI approaches.

---

## 常见追问准备

### 技术追问

**Q: 为什么R²只有0.45?**
A: 只用12个2D描述符，没有3D结构信息或蛋白-配体相互作用特征。加入MD特征后R²提升到了0.6+。QSAR的R²在0.5-0.7属于实际可用水平。

**Q: 为什么选Random Forest不选深度学习?**
A: 6,896个样本对深度学习不算大。随机森林在小到中等数据集上更稳定，可解释性更好，特征重要性直接可用。但如果样本量上万，GNN会更有优势。

**Q: Morgan指纹和描述符的区别?**
A: 描述符是固定的物理化学性质，有明确的化学含义。Morgan指纹编码了子结构存在/不存在，信息更丰富（1024 bits）但对小数据集容易过拟合。项目里对比了两者，描述符在2,000样本下R²=0.45，指纹=0.66，说明描述符其实已经够用。

**Q: PCA的自由能景观怎么算的?**
A: 先对骨架原子坐标做PCA降维到前两个主成分，然后用高斯核密度估计算出二维概率密度，再用ΔG = -RT ln(P/Pmax)转成自由能。

### 非技术追问

**Q: 你最擅长什么?**
A: 在传统计算化学（对接/MD/QM）和现代ML工具之间做桥接。我能理解药物的物理化学本质，也能写代码自动化分析流程。

**Q: 20天能做什么?**
A: 我可以独立完成一个靶点的QSAR管线搭建，或者自动化一套MD分析流程。你给我的配体库，我能从SMILES跑到活性预测。
