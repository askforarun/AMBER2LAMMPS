"""Machine-readable provenance records for AMBER2LAMMPS conversions."""

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
from datetime import datetime, timezone


SCHEMA_VERSION = "1.0"


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_version(distribution):
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def source_control():
    """Identify this source checkout without borrowing an enclosing repo."""
    try:
        module_root = os.path.dirname(os.path.realpath(__file__))
        top_level = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"], cwd=module_root, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        ).stdout.strip()
        if os.path.realpath(top_level) != module_root:
            return {"commit": None, "dirty": None, "repository_root": None}
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=module_root, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        ).stdout.strip() or None
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=module_root, check=True, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, text=True,
        ).stdout
        return {
            "commit": commit,
            "dirty": bool(status.strip()),
            "repository_root": module_root,
        }
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "dirty": None, "repository_root": None}


def file_record(path):
    absolute = os.path.abspath(path)
    return {
        "path": absolute,
        "sha256": sha256_file(absolute),
        "bytes": os.path.getsize(absolute),
    }


def _scaling_record(profile):
    if profile is None:
        return None
    return {
        "scee": profile.scee,
        "scnb": profile.scnb,
        "lammps_coulomb_14_weight": profile.coulomb_weight,
        "lammps_lj_14_weight": profile.lj_weight,
    }


def _lammps_contract(param_file, top_ctx, conn_ctx, c4_cutoff, c4_coulomb):
    """Describe the minimal LAMMPS setup needed to consume generated files."""
    profiles = [profile for profile in top_ctx.scaling_profiles if profile is not None]
    distinct = list(dict.fromkeys(profiles))
    mixed_scaling = len(distinct) > 1
    sets_pair_style = mixed_scaling or top_ctx.has_c4
    special_bonds = None
    if mixed_scaling:
        baseline = next(
            profile for profile in distinct
            if profile.lj_weight not in (0.0, 1.0)
            and profile.coulomb_weight not in (0.0, 1.0)
        )
        special_bonds = (
            f"special_bonds lj 0.0 0.0 {baseline.lj_weight:.10g} "
            f"coul 0.0 0.0 {baseline.coulomb_weight:.10g}"
        )
    elif top_ctx.system_14_scaling is not None:
        profile = top_ctx.system_14_scaling
        special_bonds = (
            f"special_bonds lj 0.0 0.0 {profile.lj_weight:.10g} "
            f"coul 0.0 0.0 {profile.coulomb_weight:.10g}"
        )

    packages = ["MOLECULE"]
    if conn_ctx.dihedral_count:
        packages.append("EXTRA-MOLECULE")
    if c4_coulomb == "long":
        packages.append("KSPACE")
    if top_ctx.has_c4:
        packages.append("LEPTON")
    return {
        "atom_style": "full",
        "include": os.path.basename(param_file),
        "sets_pair_style": sets_pair_style,
        "special_bonds": special_bonds,
        "required_packages": packages,
        "coulomb": c4_coulomb,
        "cutoff_angstrom": c4_cutoff,
        "requires_kspace": c4_coulomb == "long",
    }


def build_provenance(*, data_file, param_file, topologies, molecule_counts,
                     pdb_file, charges_target, buffer, keep_temp, c4_cutoff,
                     c4_coulomb, use_typelabels, top_ctx, conn_ctx,
                     expected_total_atoms, box_bounds, interface):
    """Build a provenance record from facts observed during conversion."""
    profiles = [_scaling_record(profile) for profile in top_ctx.scaling_profiles]
    distinct_profiles = {
        (item["scee"], item["scnb"]) for item in profiles if item is not None
    }
    if top_ctx.has_c4:
        nonbonded_strategy = "amber-12-6-4"
    elif len(distinct_profiles) > 1:
        nonbonded_strategy = "hybrid-per-scaling-profile"
    else:
        nonbonded_strategy = "single-scaling-profile"

    topology_records = []
    charge_warnings = []
    for index, (path, count, target, parm, profile) in enumerate(zip(
        topologies, molecule_counts, charges_target, top_ctx.all_parms, profiles
    )):
        raw_charge = float(sum(atom.charge for atom in parm.atoms))
        adjustment = float(target - raw_charge)
        if abs(adjustment) > 0.01:
            charge_warnings.append({
                "code": "large_charge_adjustment",
                "topology_index": index,
                "adjustment_per_molecule_e": adjustment,
            })
        topology_records.append({
            **file_record(path),
            "molecule_count": count,
            "atoms_per_molecule": len(parm.atoms),
            "raw_charge_per_molecule_e": raw_charge,
            "target_charge_per_molecule_e": float(target),
            "charge_adjustment_per_molecule_e": adjustment,
            "scaling_14": profile,
        })

    debug_files = []
    if keep_temp:
        stem = os.path.splitext(data_file)[0]
        debug_files = [file_record(f"{stem}.{name}") for name in (
            "pairs.txt", "bonds.txt", "angles.txt", "dihedrals.txt"
        )]

    source = source_control()
    return {
        "$schema": "https://raw.githubusercontent.com/askforarun/AMBER2LAMMPS/main/docs/provenance.schema.json",
        "schema_version": SCHEMA_VERSION,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "software": {
            "amber_to_lammps_version": package_version("amber-to-lammps") or "source-checkout",
            "amber_to_lammps_git_commit": source["commit"],
            "amber_to_lammps_git_dirty": source["dirty"],
            "amber_to_lammps_repository_root": source["repository_root"],
            "python": platform.python_version(),
            "parmed": package_version("ParmEd"),
            "numpy": package_version("numpy"),
            "platform": platform.platform(),
        },
        "invocation": {
            "interface": interface,
            "arguments": {
                "data_file": os.path.abspath(data_file),
                "parameter_file": os.path.abspath(param_file),
                "pdb_file": os.path.abspath(pdb_file),
                "topologies": [os.path.abspath(path) for path in topologies],
                "molecule_counts": list(molecule_counts),
                "target_charges_e": [float(value) for value in charges_target],
                "buffer_angstrom": buffer,
                "keep_debug_files": keep_temp,
                "c4_cutoff_angstrom": c4_cutoff,
                "c4_coulomb": c4_coulomb,
                "type_labels": use_typelabels,
            },
        },
        "inputs": {"pdb": file_record(pdb_file), "topologies": topology_records},
        "conversion": {
            "status": "success",
            "nonbonded_strategy": nonbonded_strategy,
            "combination_rule": "lorentz-berthelot-with-explicit-amber-overrides",
            "has_pair_specific_terms": top_ctx.has_pair_specific_terms,
            "has_c4_terms": top_ctx.has_c4,
            "warnings": charge_warnings,
            "counts": {
                "atoms": expected_total_atoms,
                "atom_types": len(top_ctx.atom_type_mapping),
                "bonds": conn_ctx.bond_count,
                "bond_types": conn_ctx.bond_type_count,
                "angles": conn_ctx.angle_count,
                "angle_types": conn_ctx.angle_type_count,
                "dihedrals": conn_ctx.dihedral_count,
                "dihedral_types": conn_ctx.dihedral_type_count,
            },
            "box_angstrom": {
                "xlo": box_bounds.xlo, "xhi": box_bounds.xhi,
                "ylo": box_bounds.ylo, "yhi": box_bounds.yhi,
                "zlo": box_bounds.zlo, "zhi": box_bounds.zhi,
            },
        },
        "lammps": _lammps_contract(
            param_file, top_ctx, conn_ctx, c4_cutoff, c4_coulomb
        ),
        "outputs": {
            "lammps_data": file_record(data_file),
            "lammps_parameters": file_record(param_file),
            "debug_files": debug_files,
        },
        "validation": {
            "status": "not_run",
            "note": "Conversion provenance is not evidence of sander or LAMMPS validation.",
        },
    }


def write_provenance(path, record):
    """Write canonical, stable JSON and return its absolute path and digest."""
    absolute = os.path.abspath(path)
    with open(absolute, "w", encoding="utf-8") as stream:
        json.dump(record, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return absolute, sha256_file(absolute)


def compact_summary(record, provenance_path=None, provenance_sha256=None):
    """Return the small summary intended for default MCP responses."""
    conversion = record["conversion"]
    return {
        "status": conversion["status"],
        "atoms": conversion["counts"]["atoms"],
        "topologies": len(record["inputs"]["topologies"]),
        "nonbonded_strategy": conversion["nonbonded_strategy"],
        "has_c4_terms": conversion["has_c4_terms"],
        "warnings": [
            {
                "code": item["code"],
                "topology": item["topology_index"],
                "delta_e": item["adjustment_per_molecule_e"],
            }
            for item in conversion["warnings"]
        ],
        "lammps": record["lammps"],
        "validation_status": record["validation"]["status"],
        "provenance_file": provenance_path,
        "provenance_sha256": provenance_sha256,
    }
