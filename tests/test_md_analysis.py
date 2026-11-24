"""Unit tests for the MD trajectory analysis toolkit (unittest only)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import mdtraj as md  # noqa: E402

import md_analysis_tool as mat  # noqa: E402


def make_trajectory(n_residues: int = 10, n_frames: int = 15, seed: int = 0) -> md.Trajectory:
    """Small synthetic protein: 4 backbone atoms per residue, random walk frames."""
    topology = md.Topology()
    chain = topology.add_chain()
    for index in range(n_residues):
        residue = topology.add_residue("ALA", chain, resSeq=index + 1)
        for name, element in (
            ("N", md.element.nitrogen),
            ("CA", md.element.carbon),
            ("C", md.element.carbon),
            ("O", md.element.oxygen),
        ):
            topology.add_atom(name, element, residue)

    rng = np.random.default_rng(seed)
    n_atoms = topology.n_atoms
    topology.create_standard_bonds()
    base = np.zeros((n_atoms, 3))
    base[:, 0] = np.arange(n_atoms) * 0.38
    xyz = (
        base[None]
        + np.cumsum(rng.normal(scale=0.01, size=(n_frames, 1, 3)), axis=0)
        + rng.normal(scale=0.03, size=(n_frames, n_atoms, 3))
    )
    return md.Trajectory(xyz, topology)


class RmsfTests(unittest.TestCase):
    def test_one_value_per_residue(self):
        traj = make_trajectory()
        backbone = mat.backbone_slice(traj)
        res_ids, rmsf = mat.per_residue_rmsf(backbone, backbone[0])
        self.assertEqual(len(res_ids), 10)
        self.assertEqual(len(rmsf), 10)
        self.assertTrue(np.all(rmsf >= 0.0))

    def test_residue_numbering_is_preserved(self):
        traj = make_trajectory()
        backbone = mat.backbone_slice(traj)
        res_ids, _ = mat.per_residue_rmsf(backbone, backbone[0])
        self.assertEqual(list(res_ids), list(range(1, 11)))


class PcaTests(unittest.TestCase):
    def test_projection_and_free_energy_ranges(self):
        traj = make_trajectory()
        backbone = mat.backbone_slice(traj)
        projection, variance, free_energy = mat.pca_free_energy(backbone, temperature=300.0)
        self.assertEqual(projection.shape, (traj.n_frames, 2))
        self.assertLessEqual(float(variance.sum()), 1.0 + 1e-9)
        self.assertGreaterEqual(float(free_energy.min()), 0.0)
        self.assertLessEqual(float(free_energy.max()), 25.0 + 1e-9)

    def test_alignment_removes_global_translation_and_rotation(self):
        traj = make_trajectory(n_frames=6)
        backbone = mat.backbone_slice(traj)
        aligned = backbone.superpose(backbone[0])
        # After superposition the first frame must match the reference exactly.
        self.assertLess(float(np.abs(aligned.xyz[0] - backbone.xyz[0]).max()), 1e-6)


class HbondTests(unittest.TestCase):
    def test_counts_are_zero_without_hydrogens(self):
        traj = make_trajectory()
        counts, times = mat.hbond_counts(traj, stride=1, ligand=None)
        self.assertEqual(len(counts), traj.n_frames)
        self.assertTrue(np.all(counts == 0.0))
        self.assertEqual(len(times), traj.n_frames)

    def test_ligand_selection_is_accepted(self):
        traj = make_trajectory()
        counts, _ = mat.hbond_counts(traj, stride=1, ligand="index 0")
        self.assertEqual(len(counts), traj.n_frames)

    def test_ligand_selection_that_matches_nothing_exits(self):
        traj = make_trajectory()
        with self.assertRaises(SystemExit):
            mat.hbond_counts(traj, stride=1, ligand="resname NOPE")


class LoadTests(unittest.TestCase):
    def test_missing_files_exit_with_message(self):
        with self.assertRaises(SystemExit):
            mat.load_trajectory(Path("missing.pdb"), Path("missing.nc"), stride=1)

    def test_backbone_selection_covers_all_atoms(self):
        traj = make_trajectory()
        self.assertEqual(mat.backbone_slice(traj).n_atoms, traj.n_atoms)


if __name__ == "__main__":
    unittest.main()
