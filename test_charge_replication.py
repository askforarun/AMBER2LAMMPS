"""Regression tests for the documented residual-charge benchmark."""

import json
import os
import unittest

from validation.charge_replication.benchmark import generate_results


HERE = os.path.dirname(os.path.abspath(__file__))


class ChargeReplicationBenchmarkTests(unittest.TestCase):

    def test_committed_results_are_reproducible(self):
        actual = generate_results()
        path = os.path.join(HERE, "validation", "charge_replication", "results.json")
        with open(path, encoding="utf-8") as stream:
            expected = json.load(stream)
        for key in (
            "fixture", "atoms_per_template", "target_charge_e", "interpretation"
        ):
            self.assertEqual(actual[key], expected[key])
        for key in (
            "source_topology_charge_e", "intermol_emitted_charge_e",
            "amber2lammps_raw_target_charge_e",
            "amber2lammps_neutral_target_charge_e",
            "amber2lammps_uniform_shift_per_atom_e",
        ):
            self.assertAlmostEqual(actual[key], expected[key], places=12)
        self.assertEqual(len(actual["replication"]), len(expected["replication"]))
        for actual_row, expected_row in zip(
            actual["replication"], expected["replication"]
        ):
            self.assertEqual(actual_row["copies"], expected_row["copies"])
            for key in actual_row.keys() - {"copies"}:
                self.assertAlmostEqual(actual_row[key], expected_row[key], places=9)

    def test_residual_scales_and_neutral_target_does_not(self):
        result = generate_results()
        residual = result["source_topology_charge_e"]
        for row in result["replication"]:
            self.assertAlmostEqual(
                row["raw_topology_total_e"], residual * row["copies"], places=12
            )
            self.assertLessEqual(
                abs(row["amber2lammps_neutral_target_total_e"]), 1.1e-6
            )
