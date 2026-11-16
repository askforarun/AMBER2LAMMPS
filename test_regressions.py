"""End-to-end regression tests for file safety, HMR masses, atom typing, and charges."""

import contextlib
import copy
import io
import json
import os
import shutil
import tempfile
import unittest

import parmed as pmd

import amber_to_lammps as a2l
import amber_to_lammps_mcp as a2l_mcp
from amber_to_lammps import amber2lammps
from validation.amber_protein_water_ions.prepare_shake_data import prepare as prepare_shake_data

HERE = os.path.dirname(os.path.abspath(__file__))
ETHANOL_PRMTOP = os.path.join(HERE, "ethanol.prmtop")
ETHANOL_PDB = os.path.join(HERE, "ethanol.pdb")
GLYCAM_PRMTOP = os.path.join(
    HERE, "validation", "amber_family_mix", "source", "glucose.prmtop"
)
MIXTURE_PDB = os.path.join(
    HERE, "validation", "amber_family_mix", "source", "ethanol_glucose.pdb"
)
DEBUG_FILES = ("pairs.txt", "bonds.txt", "angles.txt", "dihedrals.txt")


def read_sections(data_file):
    """Return {section_name: [non-empty lines]} for a LAMMPS data file."""
    sections = {}
    current = None
    with open(data_file) as f:
        for raw in f:
            line = raw.strip()
            if line in ("Masses", "Atoms", "Bonds", "Angles", "Dihedrals"):
                current = line
                sections[current] = []
            elif current and line:
                sections[current].append(line)
    return sections


class ConversionTestCase(unittest.TestCase):
    """Runs each test inside a fresh temporary working directory."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="a2l_test_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.old_cwd = os.getcwd()
        os.chdir(self.tmp)
        self.addCleanup(os.chdir, self.old_cwd)

    def convert(self, topology=ETHANOL_PRMTOP, charges=(0,), data_file="out.data",
                param_file="out.parm", **kwargs):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            amber2lammps(
                data_file, param_file, [topology], [1], ETHANOL_PDB, list(charges), **kwargs
            )
        return out.getvalue()


class McpConversionTests(ConversionTestCase):
    """The MCP-facing helper validates inputs and uses the conversion engine."""

    def test_mcp_helper_converts_and_returns_absolute_outputs(self):
        result = a2l_mcp.convert_amber_system(
            pdb_file=ETHANOL_PDB,
            topologies=[ETHANOL_PRMTOP],
            counts=[1],
            charges=[0],
            data_file="mcp.data",
            parameter_file="mcp.parm",
            type_labels=True,
        )
        self.assertEqual(result["status"], "success")
        self.assertEqual(
            os.path.realpath(result["data_file"]),
            os.path.realpath(os.path.join(self.tmp, "mcp.data")),
        )
        self.assertEqual(
            os.path.realpath(result["parameter_file"]),
            os.path.realpath(os.path.join(self.tmp, "mcp.parm")),
        )
        self.assertTrue(os.path.isfile(result["data_file"]))
        self.assertTrue(os.path.isfile(result["parameter_file"]))
        self.assertEqual(result["summary"]["atoms"], 9)
        self.assertEqual(result["summary"]["validation_status"], "not_run")
        self.assertTrue(os.path.isfile(result["summary"]["provenance_file"]))
        self.assertEqual(len(result["summary"]["provenance_sha256"]), 64)

    def test_provenance_records_hashes_without_claiming_validation(self):
        amber2lammps(
            "out.data", "out.parm", [ETHANOL_PRMTOP], [1], ETHANOL_PDB, [0],
            provenance_file="out.provenance.json",
        )
        with open("out.provenance.json", encoding="utf-8") as stream:
            record = json.load(stream)
        self.assertEqual(record["schema_version"], "1.0")
        self.assertEqual(record["conversion"]["status"], "success")
        self.assertEqual(record["conversion"]["counts"]["atoms"], 9)
        self.assertEqual(record["validation"]["status"], "not_run")
        self.assertEqual(len(record["inputs"]["pdb"]["sha256"]), 64)
        self.assertEqual(len(record["outputs"]["lammps_data"]["sha256"]), 64)

    def test_provenance_cannot_overwrite_inputs_or_outputs(self):
        with self.assertRaisesRegex(ValueError, "must not overwrite"):
            amber2lammps(
                "out.data", "out.parm", [ETHANOL_PRMTOP], [1], ETHANOL_PDB, [0],
                provenance_file="out.data",
            )
        self.assertFalse(os.path.exists("out.data"))

    def test_mcp_helper_rejects_mismatched_lists_before_writing(self):
        with self.assertRaisesRegex(ValueError, "must have the same length"):
            a2l_mcp.convert_amber_system(
                pdb_file=ETHANOL_PDB,
                topologies=[ETHANOL_PRMTOP],
                counts=[1, 2],
                charges=[0],
                data_file="mcp.data",
                parameter_file="mcp.parm",
            )
        self.assertFalse(os.path.exists("mcp.data"))


class Tip3pShakePreparationTests(ConversionTestCase):
    """AMBER water triangles become the topology required by LAMMPS SHAKE."""

    def test_hh_bonds_become_hoh_angles(self):
        fixture = os.path.join(HERE, "validation", "amber_protein_water_ions", "amber2lammps")
        count = prepare_shake_data(
            os.path.join(fixture, "system.data"),
            os.path.join(fixture, "system.parm"),
            "shake.data", "shake.parm",
        )
        self.assertEqual(count, 380)
        sections = read_sections("shake.data")
        self.assertFalse(any(line.split()[1] == "9" for line in sections["Bonds"]))
        self.assertEqual(
            sum(line.split()[1] == "12" for line in sections["Angles"]), 380
        )
        with open("shake.parm") as stream:
            self.assertIn("angle_coeff 12 0.0000 104.490", stream.read())


class AmberFamilyScalingTests(ConversionTestCase):
    """AMBER SCEE/SCNB values are detected and translated fail-closed."""

    def test_standard_amber_scaling_is_read_from_active_dihedrals(self):
        parm = pmd.load_file(ETHANOL_PRMTOP)
        scaling = a2l.topology_14_scaling(parm, "ethanol")
        self.assertEqual(scaling, a2l.SpecialBondScaling(scee=1.2, scnb=2.0))
        self.assertAlmostEqual(scaling.coulomb_weight, 1.0 / 1.2)
        self.assertAlmostEqual(scaling.lj_weight, 0.5)

    def test_improper_zero_placeholders_do_not_define_scaling(self):
        parm = pmd.load_file(os.path.join(HERE, "epon.prmtop"))
        scaling = a2l.topology_14_scaling(parm, "epon")
        self.assertEqual(scaling, a2l.SpecialBondScaling(scee=1.2, scnb=2.0))

    def test_mixed_gaff_and_glycam_style_scaling_has_no_global_profile(self):
        profiles = [
            a2l.SpecialBondScaling(scee=1.2, scnb=2.0),
            a2l.SpecialBondScaling(scee=1.0, scnb=1.0),
        ]
        self.assertIsNone(
            a2l.resolve_system_14_scaling(profiles, ["ligand.prmtop", "sugar.prmtop"])
        )
        self.assertEqual(a2l.distinct_14_scaling_profiles(profiles), profiles)

    def test_parameter_file_emits_topology_derived_special_bonds(self):
        self.convert()
        with open("out.parm") as stream:
            parameters = stream.read()
        self.assertIn(
            "special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.8333333333",
            parameters,
        )

    def test_gaff2_glycam_mixture_uses_hybrid_profiles_and_complete_matrix(self):
        amber2lammps(
            "mixed.data",
            "mixed.parm",
            [ETHANOL_PRMTOP, GLYCAM_PRMTOP],
            [1, 1],
            MIXTURE_PDB,
            [0.0000010021786493, 0.0],
        )
        with open("mixed.parm") as stream:
            parameters = stream.read()

        self.assertIn(
            "pair_style hybrid lj/cut/coul/long 9 9 lj/cut/coul/long 9 9",
            parameters,
        )
        self.assertIn(
            "pair_modify pair lj/cut/coul/long 1 special lj 0.0 0.0 0.5 "
            "special coul 0.0 0.0 0.8333333333",
            parameters,
        )
        self.assertIn(
            "pair_modify pair lj/cut/coul/long 2 special lj 0.0 0.0 1 "
            "special coul 0.0 0.0 1",
            parameters,
        )
        # Eleven atom types require 11 * 12 / 2 explicit I-J assignments.
        self.assertEqual(parameters.count("\npair_coeff "), 66)

    def test_exact_combined_pdb_count_does_not_emit_false_warning(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            amber2lammps(
                "mixed.data",
                "mixed.parm",
                [ETHANOL_PRMTOP, GLYCAM_PRMTOP],
                [1, 1],
                MIXTURE_PDB,
                [0.0000010021786493, 0.0],
            )
        self.assertNotIn("but multi-mode expects", output.getvalue())


class PdbBoxTests(ConversionTestCase):
    """CRYST1 bounds are used safely, with coordinate bounds as the fallback."""

    def pdb_with_prefix(self, prefix):
        path = os.path.join(self.tmp, "box.pdb")
        with open(ETHANOL_PDB) as source, open(path, "w") as target:
            target.write(prefix)
            target.write(source.read())
        return path

    def test_orthogonal_cryst1_sets_box_lengths(self):
        pdb = self.pdb_with_prefix(
            "CRYST1   30.000   31.000   32.000  90.00  90.00  90.00 P 1           1\n"
        )
        _atoms, _x, _y, _z, box = a2l.parse_pdb_coordinates(pdb)
        self.assertEqual(box, a2l.BoxBounds(0, 30.0, 0, 31.0, 0, 32.0))

    def test_missing_cryst1_uses_coordinate_fallback(self):
        _atoms, x, y, z, box = a2l.parse_pdb_coordinates(ETHANOL_PDB)
        self.assertIsNone(box)
        fallback = a2l.compute_box_bounds(x, y, z, 3.8)
        self.assertAlmostEqual(fallback.xlo, min(x) - 3.8)
        self.assertAlmostEqual(fallback.xhi, max(x) + 3.8)
        self.assertAlmostEqual(fallback.ylo, min(y) - 3.8)
        self.assertAlmostEqual(fallback.zhi, max(z) + 3.8)

    def test_nonorthogonal_cryst1_is_rejected(self):
        pdb = self.pdb_with_prefix(
            "CRYST1   30.000   31.000   32.000  90.00  90.00 120.00 P 1           1\n"
        )
        with self.assertRaisesRegex(ValueError, "Non-orthogonal CRYST1 cells are not supported"):
            a2l.parse_pdb_coordinates(pdb)

    def test_placeholder_cryst1_uses_coordinate_fallback(self):
        pdb = self.pdb_with_prefix(
            "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
        )
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            _atoms, _x, _y, _z, box = a2l.parse_pdb_coordinates(pdb)
        self.assertIsNone(box)
        self.assertIn("ignoring placeholder", output.getvalue())

    def test_cryst1_smaller_than_coordinate_span_is_rejected(self):
        pdb = self.pdb_with_prefix(
            "CRYST1    2.000    2.000    2.000  90.00  90.00  90.00 P 1           1\n"
        )
        with self.assertRaisesRegex(ValueError, "smaller than the coordinate span"):
            a2l.parse_pdb_coordinates(pdb)

    def test_malformed_cryst1_warns_without_verbose(self):
        pdb = self.pdb_with_prefix("CRYST1  not-a-cell\n")
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            _atoms, _x, _y, _z, box = a2l.parse_pdb_coordinates(pdb)
        self.assertIsNone(box)
        self.assertIn("Warning: Could not parse PDB CRYST1 record", output.getvalue())


class TypeLabelTests(ConversionTestCase):
    """Type labels are opt-in and preserve chemically distinct type names."""

    def test_default_output_has_no_type_label_sections(self):
        self.convert()
        with open("out.data") as f:
            data = f.read()
        self.assertNotIn("Type Labels", data)

    def test_typelabel_output_has_all_label_sections(self):
        self.convert(use_typelabels=True)
        with open("out.data") as f:
            data = f.read()
        for section in (
            "Atom Type Labels", "Bond Type Labels", "Angle Type Labels",
            "Dihedral Type Labels",
        ):
            self.assertIn(section, data)

    def test_invalid_and_colliding_atom_names_are_sanitized_uniquely(self):
        labels = a2l.build_atom_type_labels({"2oh": 1, "type_2oh": 2, "bad type": 3})
        self.assertEqual(labels["2oh"], "type_2oh")
        self.assertEqual(labels["type_2oh"], "type_2oh_2")
        self.assertEqual(labels["bad type"], "bad_type")
        for label in labels.values():
            self.assertFalse(label[0].isdigit() or label[0] in "*#")
            self.assertFalse(any(char.isspace() for char in label))

    def test_digit_leading_atom_type_produces_valid_atom_and_bonded_labels(self):
        parm = pmd.amber.AmberParm(ETHANOL_PRMTOP)
        for atom in parm.atoms:
            if atom.type == "oh":
                atom.type = "2oh"
        topology = os.path.join(self.tmp, "digit_type.prmtop")
        parm.write_parm(topology)
        self.convert(topology=topology, use_typelabels=True)

        labels = []
        in_labels = False
        with open("out.data") as f:
            for raw in f:
                line = raw.strip()
                if line == "Atom Type Labels":
                    in_labels = True
                elif line == "Masses":
                    break
                elif in_labels:
                    fields = line.split()
                    if len(fields) == 2 and fields[0].isdigit():
                        labels.append(fields[1])

        self.assertTrue(labels)
        for label in labels:
            self.assertFalse(label[0].isdigit() or label[0] in "*#")
            self.assertFalse(any(char.isspace() for char in label))

    def test_equal_parameters_keep_distinct_bond_labels(self):
        parm = pmd.amber.AmberParm(ETHANOL_PRMTOP)
        first = parm.bonds[0]
        second = next(
            bond for bond in parm.bonds
            if {bond.atom1.type, bond.atom2.type} != {first.atom1.type, first.atom2.type}
        )
        second.type = first.type

        numeric = a2l.process_connectivity([parm], [1], [ETHANOL_PRMTOP], False)
        labelled = a2l.process_connectivity([parm], [1], [ETHANOL_PRMTOP], True)
        self.assertGreater(labelled.bond_type_count, numeric.bond_type_count)

        expected = {
            "-".join(sorted((first.atom1.type, first.atom2.type))),
            "-".join(sorted((second.atom1.type, second.atom2.type))),
        }
        actual = {line.split(maxsplit=1)[1] for line in labelled.bond_type_label_lines}
        self.assertTrue(expected.issubset(actual))

    def test_same_bond_names_with_different_parameters_get_unique_ids_and_labels(self):
        first = pmd.amber.AmberParm(ETHANOL_PRMTOP)
        second = pmd.amber.AmberParm(ETHANOL_PRMTOP)
        second.bonds[0].type = copy.copy(second.bonds[0].type)
        second.bonds[0].type.k += 1.0

        labelled = a2l.process_connectivity(
            [first, second], [1, 1], [ETHANOL_PRMTOP, ETHANOL_PRMTOP], True
        )
        labels = [line.split(maxsplit=1)[1] for line in labelled.bond_type_label_lines]
        self.assertEqual(len(labels), len(set(labels)))

        first_id = labelled.bond_lines[0].split()[1]
        second_topology_start = len(first.bonds)
        second_id = labelled.bond_lines[second_topology_start].split()[1]
        self.assertNotEqual(first_id, second_id)


class DebugFileSafetyTests(ConversionTestCase):
    """Issue #1: conversions must never delete or overwrite unrelated files."""

    def write_user_files(self):
        for name in DEBUG_FILES:
            with open(name, "w") as f:
                f.write("user data\n")

    def assert_user_files_intact(self):
        for name in DEBUG_FILES:
            self.assertTrue(os.path.exists(name), f"{name} was deleted")
            with open(name) as f:
                self.assertEqual(f.read(), "user data\n", f"{name} was modified")

    def test_existing_debug_named_files_are_not_deleted(self):
        self.write_user_files()
        self.convert()
        self.assert_user_files_intact()

    def test_keep_temp_does_not_overwrite_existing_debug_named_files(self):
        self.write_user_files()
        self.convert(keep_temp=True)
        self.assert_user_files_intact()

    def test_default_run_creates_no_debug_files(self):
        self.convert()
        for name in DEBUG_FILES:
            self.assertFalse(os.path.exists(name), f"{name} written without --keep-temp")
            self.assertFalse(os.path.exists(f"out.{name}"), f"out.{name} written without --keep-temp")

    def test_keep_temp_writes_stem_prefixed_files_beside_data_file(self):
        os.mkdir("results")
        data_file = os.path.join("results", "system.data")
        output = self.convert(keep_temp=True, data_file=data_file,
                              param_file=os.path.join("results", "system.parm"))
        for name in DEBUG_FILES:
            path = os.path.join("results", f"system.{name}")
            self.assertTrue(os.path.exists(path), f"{path} missing with --keep-temp")
            self.assertIn(path, output)
            self.assertFalse(os.path.exists(name), f"bare {name} written to working directory")
        with open(os.path.join("results", "system.bonds.txt")) as f:
            self.assertTrue(f.readline().startswith("# system.bonds.txt generated"))

    def test_keep_temp_refuses_to_overwrite_parameter_file(self):
        with self.assertRaisesRegex(ValueError, "would overwrite an output file"):
            self.convert(keep_temp=True, data_file="out.data", param_file="out.bonds.txt")
        self.assertFalse(os.path.exists("out.data"), "data file written before the collision check")


class HydrogenMassRepartitioningTests(ConversionTestCase):
    """Issue #3: atoms sharing an AMBER type name but differing in mass get their own types."""

    def setUp(self):
        super().setUp()
        parm = pmd.amber.AmberParm(ETHANOL_PRMTOP)
        pmd.tools.actions.HMassRepartition(parm).execute()
        self.hmr_prmtop = os.path.join(self.tmp, "ethanol_hmr.prmtop")
        parm.write_parm(self.hmr_prmtop)
        self.hmr_parm = pmd.amber.AmberParm(self.hmr_prmtop)

    def test_hmr_converts_with_per_atom_masses(self):
        self.convert(topology=self.hmr_prmtop)
        sections = read_sections("out.data")

        type_mass = {}
        for line in sections["Masses"]:
            type_id, mass = line.split()[:2]
            type_mass[int(type_id)] = float(mass)

        atom_rows = sections["Atoms"]
        self.assertEqual(len(atom_rows), len(self.hmr_parm.atoms))
        for row, atom in zip(atom_rows, self.hmr_parm.atoms):
            type_id = int(row.split()[2])
            self.assertAlmostEqual(type_mass[type_id], atom.mass, places=3,
                                   msg=f"atom {atom.idx + 1} ({atom.type}) has wrong mass")

        total = sum(type_mass[int(row.split()[2])] for row in atom_rows)
        self.assertAlmostEqual(total, sum(a.mass for a in self.hmr_parm.atoms), places=2)

    def test_hmr_mass_variants_share_lj_parameters(self):
        self.convert(topology=self.hmr_prmtop)
        pair = {}
        with open("out.parm") as f:
            for line in f:
                if line.startswith("pair_coeff"):
                    fields = line.split()
                    name = line.split("#", 1)[1].strip()
                    pair[name] = (float(fields[3]), float(fields[4]))

        # Ethanol's two c3 carbons differ in mass under HMR: c3 and c3_top1.
        self.assertIn("c3", pair)
        self.assertIn("c3_top1", pair)
        self.assertEqual(pair["c3"], pair["c3_top1"])

    def test_hmr_connectivity_labels_use_renamed_atom_types(self):
        self.convert(topology=self.hmr_prmtop, use_typelabels=True)
        with open("out.data") as f:
            data = f.read()
        self.assertIn("c3-c3_top1", data)

    def test_standard_topology_has_one_type_per_amber_type(self):
        self.convert()
        parm = pmd.amber.AmberParm(ETHANOL_PRMTOP)
        sections = read_sections("out.data")
        self.assertEqual(len(sections["Masses"]), len({a.type for a in parm.atoms}))


class UnknownAtomTypeTests(unittest.TestCase):
    """Issue #5: unmapped atoms must raise instead of silently becoming type 1."""

    def test_unparsed_atoms_raise_during_topology_loading(self):
        def empty_details(parm, mask):
            return ""

        with self.assertRaisesRegex(ValueError, "Could not read atom type"):
            a2l.load_topology_context(
                [ETHANOL_PRMTOP], False, pmd.amber.AmberParm, empty_details
            )

    def test_unmapped_atom_raises_when_writing_atoms(self):
        parm = pmd.amber.AmberParm(ETHANOL_PRMTOP)
        n_atoms = len(parm.atoms)
        pdb_atoms = [{"x": 0.0, "y": 0.0, "z": 0.0} for _ in range(n_atoms)]
        spans = a2l.build_molecule_spans([parm], [1])
        type_remap = {idx: atom.type for idx, atom in enumerate(parm.atoms)}
        del type_remap[2]
        mapping = {name: i + 1 for i, name in enumerate(sorted(set(type_remap.values())))}

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "No LAMMPS atom type for atom 3"):
                a2l.write_atoms_section(
                    os.path.join(tmp, "out.data"), pdb_atoms, spans, [parm],
                    [type_remap], mapping, [0.0] * n_atoms, False,
                )


class ChargeWarningTests(ConversionTestCase):
    """Issue #2: keep the requested target, but always warn about large corrections."""

    def atom_charges(self):
        return [float(row.split()[3]) for row in read_sections("out.data")["Atoms"]]

    def test_large_correction_warns_without_verbose_and_honours_target(self):
        output = self.convert(charges=(1,))
        self.assertIn("Warning: topology 1", output)
        self.assertIn("--charges requests 1.000000 e", output)
        self.assertAlmostEqual(sum(self.atom_charges()), 1.0, places=6)

    def test_routine_neutralization_does_not_warn(self):
        parm = pmd.amber.AmberParm(ETHANOL_PRMTOP)
        residual = abs(sum(a.charge for a in parm.atoms))
        self.assertLess(residual, 0.01)

        output = self.convert(charges=(0,))
        self.assertNotIn("Warning", output)
        self.assertAlmostEqual(sum(self.atom_charges()), 0.0, places=6)


if __name__ == "__main__":
    unittest.main()
