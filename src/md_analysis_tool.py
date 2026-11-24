#!/usr/bin/env python3
"""MD trajectory analysis toolkit.

Computes backbone RMSD, per-residue RMSF, radius of gyration, a PCA-based
free-energy landscape and hydrogen-bond statistics, then writes one six-panel
summary figure (plus an optional per-frame metrics table).

AMBER (.prmtop + .nc) and GROMACS (.tpr/.gro + .xtc) inputs are supported
through MDTraj.

Examples
--------
    python md_analysis_tool.py --top 1.prmtop --traj md.nc --out ../results
    python md_analysis_tool.py --top topol.tpr --traj traj.xtc \
        --ligand "resname LIG" --temperature 310 --label "GLB-CYP2C9"
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA

R_KJ = 0.0083144621  # gas constant, kJ mol-1 K-1


def load_trajectory(top: Path, traj: Path, stride: int) -> md.Trajectory:
    """Load a trajectory with an explicit topology and an optional stride."""
    if not top.exists():
        raise SystemExit(f"Topology not found: {top}")
    if not traj.exists():
        raise SystemExit(f"Trajectory not found: {traj}")
    t = md.load(str(traj), top=str(top), stride=stride)
    if t.n_frames < 2:
        raise SystemExit("Need at least two frames after striding.")
    return t


def backbone_slice(traj: md.Trajectory) -> md.Trajectory:
    """Keep protein backbone atoms (N, CA, C, O)."""
    indices = traj.topology.select("protein and backbone")
    if len(indices) == 0:
        raise SystemExit("No protein backbone atoms found in the topology.")
    return traj.atom_slice(indices)


def per_residue_rmsf(traj_bb: md.Trajectory, ref: md.Trajectory):
    """Per-atom RMSF averaged into one value per residue."""
    rmsf_atom = md.rmsf(traj_bb, ref)
    res_index = np.array([atom.residue.index for atom in traj_bb.topology.atoms])
    n_res = traj_bb.topology.n_residues
    counts = np.bincount(res_index, minlength=n_res)
    sums = np.bincount(res_index, weights=rmsf_atom, minlength=n_res)
    rmsf_res = sums / np.maximum(counts, 1)
    res_ids = np.array([res.resSeq for res in traj_bb.topology.residues])
    return res_ids, rmsf_res


def pca_free_energy(traj_bb: md.Trajectory, temperature: float):
    """Align every frame onto the first one, then build PC1/PC2 and the FEL."""
    aligned = traj_bb.superpose(traj_bb[0])
    xyz = aligned.xyz.reshape(aligned.n_frames, -1)
    pca = PCA(n_components=2, random_state=0)
    projection = pca.fit_transform(xyz)

    density = gaussian_kde(projection.T)(projection.T)
    rt = R_KJ * (273.15 + temperature)
    free_energy = -rt * np.log(density / density.max())
    free_energy = np.clip(free_energy, 0.0, 25.0)
    return projection, pca.explained_variance_ratio_, free_energy


def hbond_counts(traj: md.Trajectory, stride: int, ligand: str | None):
    """Hydrogen bonds per frame; with --ligand, only protein-ligand bridges.

    Wernet-Nilsson returns a list of (donor, hydrogen, acceptor) triplets per
    frame, so the whole trajectory is analysed in a single pass.
    """
    sub = traj[::stride]
    frames = md.wernet_nilsson(sub, periodic=False)
    if ligand:
        ligand_atoms = set(traj.topology.select(ligand).tolist())
        if not ligand_atoms:
            raise SystemExit(f"Selection matched no atoms: {ligand}")
        counts = np.array([
            sum(1 for donor, _, acceptor in frame
                if (donor in ligand_atoms) != (acceptor in ligand_atoms))
            for frame in frames
        ], dtype=float)
    else:
        counts = np.array([len(frame) for frame in frames], dtype=float)
    return counts, sub.time / 1000.0


def plot_summary(time_ns, rmsd, res_ids, rmsf_res, rg, projection, var_ratio,
                 free_energy, time_hb, hbonds, temperature, label, out_png: Path):
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    ax = axes[0, 0]
    ax.plot(time_ns, rmsd, lw=0.5, color="steelblue")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (nm)")
    ax.set_title(f"Backbone RMSD (max = {rmsd.max():.2f} nm)")

    ax = axes[0, 1]
    ax.plot(res_ids, rmsf_res, lw=0.6, color="steelblue")
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RMSF (nm)")
    ax.set_title("Per-residue RMSF")

    ax = axes[0, 2]
    ax.plot(time_ns, rg, lw=0.5, color="darkorange")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Rg (nm)")
    ax.set_title("Radius of gyration")

    ax = axes[1, 0]
    sc = ax.scatter(projection[:, 0], projection[:, 1], c=free_energy,
                    cmap="viridis", s=2, alpha=0.8)
    ax.set_xlabel(f"PC1 ({var_ratio[0] * 100:.0f}%)")
    ax.set_ylabel(f"PC2 ({var_ratio[1] * 100:.0f}%)")
    ax.set_title("Free-energy landscape")
    plt.colorbar(sc, ax=ax, label="Free energy (kJ/mol)")

    ax = axes[1, 1]
    ax.plot(time_ns, projection[:, 0], lw=0.5, label="PC1")
    ax.plot(time_ns, projection[:, 1], lw=0.5, label="PC2")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("PC score")
    ax.set_title("PCA projection")
    ax.legend()

    ax = axes[1, 2]
    ax.hist(hbonds, bins=30, color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(hbonds.mean(), color="red", ls="--")
    ax.set_xlabel("H-bond count")
    ax.set_ylabel("Frames")
    ax.set_title(f"H-bonds (mean = {hbonds.mean():.0f})")

    if label:
        fig.suptitle(label, fontsize=14)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--top", required=True, type=Path,
                        help="topology: .prmtop, .tpr, .gro or .pdb")
    parser.add_argument("--traj", required=True, type=Path,
                        help="trajectory: .nc, .xtc or .dcd")
    parser.add_argument("--out", default="../results", type=Path,
                        help="output directory")
    parser.add_argument("--stride", default=50, type=int,
                        help="keep every Nth frame when loading (default: 50)")
    parser.add_argument("--hb-stride", default=4, type=int,
                        help="extra stride for H-bond sampling (default: 4)")
    parser.add_argument("--temperature", default=300.0, type=float,
                        help="temperature for the free-energy scale (K)")
    parser.add_argument("--ligand", default=None,
                        help="MDTraj selection for the ligand, e.g. 'resname LIG'")
    parser.add_argument("--label", default="", help="figure title")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("MD trajectory analysis")
    print("=" * 60)

    print("\n[1] Loading trajectory ...")
    traj = load_trajectory(args.top, args.traj, args.stride)
    print(f"  Atoms: {traj.n_atoms}  Frames: {traj.n_frames}  "
          f"Length: {traj.time[-1] / 1000:.1f} ns")

    traj_bb = backbone_slice(traj)
    print(f"  Backbone atoms: {traj_bb.n_atoms}  Residues: {traj_bb.n_residues}")
    ref = traj_bb[0]
    time_ns = traj_bb.time / 1000.0

    print("\n[2] RMSD ...")
    rmsd = md.rmsd(traj_bb, ref)
    print(f"  {rmsd.min():.2f} - {rmsd.max():.2f} nm")

    print("[3] Per-residue RMSF ...")
    res_ids, rmsf_res = per_residue_rmsf(traj_bb, ref)
    top5 = np.argsort(rmsf_res)[-5:][::-1]
    print("  Most flexible: " + ", ".join(f"{res_ids[i]}({rmsf_res[i]:.2f} nm)"
                                          for i in top5))

    print("[4] Radius of gyration ...")
    rg = md.compute_rg(traj_bb)

    print("[5] PCA + free-energy landscape ...")
    projection, var_ratio, free_energy = pca_free_energy(traj_bb, args.temperature)
    print(f"  PC1 {var_ratio[0] * 100:.1f}%  PC2 {var_ratio[1] * 100:.1f}%")

    print("[6] Hydrogen bonds ...")
    hbonds, time_hb = hbond_counts(traj, args.hb_stride, args.ligand)
    print(f"  Mean {hbonds.mean():.1f} per frame")

    print("[7] Writing output ...")
    fig_path = args.out / "md_analysis_summary.png"
    plot_summary(time_ns, rmsd, res_ids, rmsf_res, rg, projection, var_ratio,
                 free_energy, time_hb, hbonds, args.temperature, args.label, fig_path)
    print(f"  Saved: {fig_path}")

    metrics_path = args.out / "md_analysis_metrics.csv"
    np.savetxt(
        metrics_path,
        np.column_stack([time_ns, rmsd, rg, projection[:, 0], projection[:, 1]]),
        delimiter=",", header="time_ns,rmsd_nm,rg_nm,pc1,pc2", comments="",
    )
    print(f"  Saved: {metrics_path}")

    print("\n" + "=" * 60)
    stability = "stable" if rmsd.max() < 0.4 else "structural drift detected"
    print(f"Backbone: {stability}")
    print(f"Top 2 PCs explain {var_ratio.sum() * 100:.0f}% of the motion")
    print("=" * 60)


if __name__ == "__main__":
    main()
