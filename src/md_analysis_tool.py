"""
MD Trajectory Analysis Tool
============================
自动化分析AMBER/GROMACS MD轨迹: RMSD/RMSF/Rg/PCA/FEL/H-bond

Author: HuangJunbao
Date: 2026-06
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import gaussian_kde
import warnings
warnings.filterwarnings("ignore")

# ============================================
# 1. 加载轨迹
# ============================================
print("=" * 60)
print("MD Trajectory Analysis")
print("=" * 60)

print("\n[1] Loading trajectory...")
traj = md.load("../../md.nc", top="../../1.prmtop")
print(f"  Atoms: {traj.n_atoms}, Frames: {traj.n_frames}")
print(f"  Length: {traj.time[-1]/1000:.1f} ns")

# 提取蛋白骨架
backbone = traj.topology.select("protein and backbone")
traj_bb = traj.atom_slice(backbone)
print(f"  Backbone atoms: {traj_bb.n_atoms}")

# 减采样 (每50帧取1帧)
traj_sub = traj_bb[::50]
time_ns = traj.time[::50] / 1000
print(f"  Subsampled frames: {traj_sub.n_frames} (stride=50)")

# ============================================
# 2. RMSD
# ============================================
print("\n[2] Computing RMSD...")
ref = traj_sub[0]
rmsd = md.rmsd(traj_sub, ref)
print(f"  RMSD: {rmsd.min():.2f} ~ {rmsd.max():.2f} nm")

# ============================================
# 3. RMSF (per-residue)
# ============================================
print("[3] Computing RMSF...")
rmsf = md.rmsf(traj_sub, ref)
residues = [atom.residue.resSeq for atom in traj_bb.topology.atoms]
print(f"  RMSF: {rmsf.min():.2f} ~ {rmsf.max():.2f} nm")

# ============================================
# 4. Radius of Gyration
# ============================================
print("[4] Computing Rg...")
rg = md.compute_rg(traj_sub)
print(f"  Rg: {rg.min():.2f} ~ {rg.max():.2f} nm")

# ============================================
# 5. PCA + Free Energy Landscape
# ============================================
print("[5] Computing PCA + FEL...")
xyz = traj_sub.xyz.reshape(traj_sub.n_frames, traj_sub.n_atoms * 3)
xyz_centered = xyz - xyz.mean(axis=0)

pca = PCA(n_components=2)
projected = pca.fit_transform(xyz_centered)

# 自由能
kT = 2.5  # kJ/mol at 300K
kde = gaussian_kde(np.vstack([projected[:, 0], projected[:, 1]]))
z = kde(np.vstack([projected[:, 0], projected[:, 1]]))
free_energy = -kT * np.log(z / z.max())
free_energy[free_energy > 20] = 20

print(f"  PC1: {pca.explained_variance_ratio_[0]*100:.1f}%")
print(f"  PC2: {pca.explained_variance_ratio_[1]*100:.1f}%")

# ============================================
# 6. H-bond count (protein internal)
# ============================================
print("[6] Computing H-bonds...")
traj_hb = traj[::200]
n_hbonds = np.zeros(traj_hb.n_frames)
for i in range(traj_hb.n_frames):
    n_hbonds[i] = len(md.baker_hubbard(traj_hb[i], periodic=False))
time_hb = traj_hb.time / 1000
print(f"  H-bonds: {n_hbonds.min():.0f} ~ {n_hbonds.max():.0f}, mean={n_hbonds.mean():.0f}")

# ============================================
# 7. 可视化 - 六合一图
# ============================================
print("[7] Generating figure...")
fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# A: RMSD
axes[0, 0].plot(time_ns, rmsd, linewidth=0.5, color="steelblue")
axes[0, 0].set_xlabel("Time (ns)")
axes[0, 0].set_ylabel("RMSD (nm)")
axes[0, 0].set_title(f"Backbone RMSD (max={rmsd.max():.2f} nm)")

# B: RMSF
axes[0, 1].plot(residues[:len(rmsf)], rmsf, linewidth=0.5, color="steelblue")
axes[0, 1].set_xlabel("Residue")
axes[0, 1].set_ylabel("RMSF (nm)")
axes[0, 1].set_title(f"Per-Residue Fluctuation")

# C: Rg + RMSD
ax_c = axes[0, 2]
ax_c.plot(time_ns, rmsd, linewidth=0.5, color="steelblue", alpha=0.7, label="RMSD")
ax_c.set_ylabel("RMSD (nm)", color="steelblue")
ax_c.set_xlabel("Time (ns)")
ax_c2 = ax_c.twinx()
ax_c2.plot(time_ns, rg, linewidth=0.5, color="darkorange", alpha=0.7, label="Rg")
ax_c2.set_ylabel("Rg (nm)", color="darkorange")
ax_c.set_title(f"RMSD vs Rg")

# D: Free Energy Landscape
sc_d = axes[1, 0].scatter(projected[:, 0], projected[:, 1], c=free_energy,
                           cmap="viridis", s=2, alpha=0.7)
axes[1, 0].set_xlabel("PC1")
axes[1, 0].set_ylabel("PC2")
axes[1, 0].set_title(f"FEL (PC1={pca.explained_variance_ratio_[0]*100:.0f}%, "
                     f"PC2={pca.explained_variance_ratio_[1]*100:.0f}%)")
plt.colorbar(sc_d, ax=axes[1, 0], label="Free Energy (kJ/mol)")

# E: PCA time projection
axes[1, 1].plot(time_ns, projected[:, 0], linewidth=0.5, label="PC1")
axes[1, 1].plot(time_ns, projected[:, 1], linewidth=0.5, label="PC2")
axes[1, 1].set_xlabel("Time (ns)")
axes[1, 1].set_ylabel("PC score")
axes[1, 1].set_title("PCA Projection Over Time")
axes[1, 1].legend()

# F: H-bond distribution
axes[1, 2].hist(n_hbonds, bins=30, color="steelblue", edgecolor="white", alpha=0.8)
axes[1, 2].axvline(x=n_hbonds.mean(), color="red", linestyle="--")
axes[1, 2].set_xlabel("H-bond Count")
axes[1, 2].set_title(f"H-bond Distribution (mean={n_hbonds.mean():.0f})")

plt.tight_layout()
plt.savefig("../results/md_analysis_summary.png", dpi=200, bbox_inches="tight")
plt.close()
print("   Saved: results/md_analysis_summary.png")

print(f"\n{'='*60}")
print("MD Analysis Complete!")
print(f"  RMSD stable? {'Yes' if rmsd.max() < 0.4 else 'Structural change detected'}")
print(f"  Most flexible region residues: {np.argsort(rmsf)[-5:]}")
print(f"  Top 2 PCs explain: {pca.explained_variance_ratio_.sum()*100:.0f}% of total motion")
print(f"  Average H-bonds: {n_hbonds.mean():.0f}")
print(f"{'='*60}")
