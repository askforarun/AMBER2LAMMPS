"""End-to-end regression tests for file safety, HMR masses, atom typing, and charges."""

import contextlib
import io
import os
import shutil
import tempfile
import unittest

import parmed as pmd

import amber_to_lammps as a2l
from amber_to_lammps import amber2lammps

HERE = os.path.dirname(os.path.abspath(__file__))
ETHANOL_PRMTOP = os.path.join(HERE, "ethanol.prmtop")
ETHANOL_PDB = os.path.join(HERE, "ethanol.pdb")
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
