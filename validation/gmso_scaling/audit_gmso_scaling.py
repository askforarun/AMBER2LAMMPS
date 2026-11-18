#!/usr/bin/env python3
"""Audit whether ParmEd-to-GMSO import preserves AMBER 1-4 scaling.

This is a characterization script, not an AMBER2LAMMPS regression test.  It
reports the behavior of the installed GMSO version without assuming that the
same behavior applies to earlier or later releases.
"""

from __future__ import annotations

import argparse
import collections
import collections.abc
import json
import platform
from pathlib import Path

import numpy as np
import parmed

# GMSO 0.18.0 imports mBuild 0.10.9 from ``gmso.external``. That mBuild
# release still refers to NumPy's removed ``np.float`` alias. Keep this
# compatibility shim local to the audit so the published Conda packages can
# be exercised without modifying either installation.
if not hasattr(np, "float"):
    np.float = float
for _name in ("Iterable", "Mapping", "MutableMapping", "MutableSet", "Sequence"):
    if not hasattr(collections, _name):
        setattr(collections, _name, getattr(collections.abc, _name))

import gmso
from gmso.external.convert_parmed import from_parmed


ROOT = Path(__file__).resolve().parents[2]
FIXTURES = (
    ("GLYCAM06j-1 glucose", "validation/glycam_glucose/work/reference.prmtop"),
    ("ff14SB peptide + GAFF2 ethanol", "validation/protein_ligand/work/reference.prmtop"),
    ("GAFF2 ethanol + GLYCAM06j-1 glucose", "validation/amber_family_mix/amber_reference.prmtop"),
)


def amber_profiles(structure: parmed.Structure) -> list[dict[str, float]]:
    """Return active proper-dihedral SCEE/SCNB divisors and their weights."""
    divisors: set[tuple[float, float]] = set()
    for dihedral in structure.dihedrals:
        if dihedral.improper or dihedral.type is None:
            continue
        terms = dihedral.type if isinstance(dihedral.type, list) else [dihedral.type]
        for term in terms:
            divisors.add((float(term.scee), float(term.scnb)))
    return [
        {
            "scee_divisor": scee,
            "scnb_divisor": scnb,
            "coulomb_14_weight": 1.0 / scee,
            "lj_14_weight": 1.0 / scnb,
        }
        for scee, scnb in sorted(divisors)
    ]


def as_list(values):
    return np.asarray(values, dtype=float).tolist()


def audit_fixture(label: str, relative_path: str) -> dict:
    path = ROOT / relative_path
    structure = parmed.load_file(str(path))
    # GMSO 0.14.0 requires coordinates during ParmEd conversion. Coordinates
    # do not affect this topology-scaling audit.
    if structure.coordinates is None:
        structure.coordinates = np.zeros((len(structure.atoms), 3))
    topology = from_parmed(structure)
    profiles = amber_profiles(structure)
    imported_lj = as_list(topology.get_lj_scale())
    imported_coulomb = as_list(topology.get_electrostatics_scale())
    molecule_scales = {
        str(key): as_list(value)
        for key, value in topology.molecule_scaling_factors.items()
    }

    preserved = False
    if len(profiles) == 1:
        expected = profiles[0]
        preserved = bool(
            np.isclose(imported_lj[2], expected["lj_14_weight"])
            and np.isclose(imported_coulomb[2], expected["coulomb_14_weight"])
        )
    elif molecule_scales:
        # Multiple profiles require molecule-specific values. This deliberately
        # does not guess a molecule-to-profile mapping.
        preserved = False

    return {
        "label": label,
        "topology": relative_path,
        "atoms": len(structure.atoms),
        "amber_active_profiles": profiles,
        "gmso_global_lj_scaling": imported_lj,
        "gmso_global_electrostatics_scaling": imported_coulomb,
        "gmso_molecule_scaling_factors": molecule_scales,
        "preserved_exactly": preserved,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the JSON report here")
    args = parser.parse_args()

    report = {
        "scope": "ParmEd-to-GMSO AMBER 1-4 scaling import audit",
        "python": platform.python_version(),
        "gmso": getattr(gmso, "__version__", "unknown"),
        "parmed": getattr(parmed, "__version__", "unknown"),
        "numpy": np.__version__,
        "fixtures": [audit_fixture(*fixture) for fixture in FIXTURES],
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
