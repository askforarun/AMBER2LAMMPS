"""Local MCP server for AMBER2LAMMPS.

The MCP dependency is imported only when the server starts, keeping the core
converter and its Python 3.9+ installation independent of MCP.
"""

import contextlib
import io
import json
import os
import shutil
import subprocess
import threading
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from amber_to_lammps import amber2lammps
from amber_to_lammps_provenance import sha256_file, write_provenance


_CONVERSION_LOCK = threading.Lock()
_VALIDATION_LOCK = threading.Lock()


def _resolve_output(path: str, root: str) -> str:
    """Resolve an MCP output beneath its configured workspace."""
    root = os.path.realpath(root)
    candidate = path if os.path.isabs(path) else os.path.join(root, path)
    candidate = os.path.realpath(candidate)
    try:
        inside = os.path.commonpath((root, candidate)) == root
    except ValueError:
        inside = False
    if not inside:
        raise ValueError(
            f"Output path '{path}' is outside the MCP workspace '{root}'. "
            "Set AMBER2LAMMPS_MCP_ROOT when starting the server to choose another workspace."
        )
    return candidate


def _relative(path: str, root: str) -> str:
    return os.path.relpath(os.path.realpath(path), os.path.realpath(root))


def _compact_conversion(conversion, root, data_file, parameter_file,
                        provenance_file, debug_files):
    summary = conversion["summary"]
    return {
        "status": "success",
        "atoms": summary["atoms"],
        "topologies": summary["topologies"],
        "nonbonded_strategy": summary["nonbonded_strategy"],
        "warnings": summary["warnings"],
        "lammps": summary["lammps"],
        "output_directory": ".",
        "files": {
            "data": _relative(data_file, root),
            "parameters": _relative(parameter_file, root),
            "provenance": _relative(provenance_file, root),
            "debug": [_relative(path, root) for path in debug_files],
        },
        "validation_status": summary["validation_status"],
    }


def convert_amber_system(
    pdb_file: str,
    topologies: List[str],
    counts: List[int],
    charges: List[float],
    data_file: str,
    parameter_file: str,
    buffer: float = 3.8,
    type_labels: bool = False,
    keep_debug_files: bool = False,
    c4_cutoff: float = 9.0,
    c4_coulomb: str = "long",
    provenance_file: Optional[str] = None,
    overwrite: bool = False,
    allow_pdb_name_mismatch: bool = False,
    output_root: Optional[str] = None,
    detail: str = "summary",
) -> Dict[str, Any]:
    """Convert AMBER topology/PDB inputs and return a structured result.

    Outputs are confined to ``AMBER2LAMMPS_MCP_ROOT`` (the server's current
    directory by default). Existing files require explicit overwrite consent.
    """
    if not topologies:
        raise ValueError("At least one AMBER topology is required.")
    if len(topologies) != len(counts) or len(topologies) != len(charges):
        raise ValueError("topologies, counts, and charges must have the same length.")
    if any(count <= 0 for count in counts):
        raise ValueError("Every molecule count must be a positive integer.")
    if c4_coulomb not in ("long", "cut"):
        raise ValueError("c4_coulomb must be 'long' or 'cut'.")
    if detail not in ("summary", "full"):
        raise ValueError("detail must be 'summary' or 'full'.")

    inputs = [pdb_file, *topologies]
    missing = [path for path in inputs if not os.path.isfile(path)]
    if missing:
        raise ValueError("Input file(s) not found: " + ", ".join(missing))

    root = output_root or os.environ.get("AMBER2LAMMPS_MCP_ROOT", os.getcwd())
    if not os.path.isdir(root):
        raise ValueError(f"MCP workspace does not exist or is not a directory: {root}")
    data_file = _resolve_output(data_file, root)
    parameter_file = _resolve_output(parameter_file, root)
    if provenance_file is None:
        provenance_file = os.path.splitext(data_file)[0] + ".provenance.json"
    provenance_file = _resolve_output(provenance_file, root)

    planned_outputs = [data_file, parameter_file, provenance_file]
    if keep_debug_files:
        stem = os.path.splitext(data_file)[0]
        planned_outputs.extend(f"{stem}.{name}" for name in (
            "pairs.txt", "bonds.txt", "angles.txt", "dihedrals.txt"
        ))
    input_paths = {os.path.realpath(path) for path in inputs}
    collisions = [path for path in planned_outputs if os.path.realpath(path) in input_paths]
    if collisions:
        raise ValueError("Output path would overwrite an input file: " + ", ".join(collisions))
    output = io.StringIO()
    # redirect_stdout changes process-global state. MCP executes synchronous
    # tools in worker threads, so serialize until converter messages are
    # returned structurally instead of printed.
    with _CONVERSION_LOCK:
        # Check under the same lock as writing so simultaneous calls targeting
        # the same new path cannot both pass the overwrite guard.
        existing = [path for path in planned_outputs if os.path.exists(path)]
        if existing and not overwrite:
            raise FileExistsError(
                "Output file(s) already exist; set overwrite=true to replace them: "
                + ", ".join(existing)
            )
        with contextlib.redirect_stdout(output):
            conversion = amber2lammps(
                data_file=data_file,
                param_file=parameter_file,
                topologies=topologies,
                molecule_counts=counts,
                pdb_file=pdb_file,
                charges_target=charges,
                buffer=buffer,
                verbose=False,
                keep_temp=keep_debug_files,
                c4_cutoff=c4_cutoff,
                c4_coulomb=c4_coulomb,
                use_typelabels=type_labels,
                provenance_file=provenance_file,
                provenance_interface="mcp",
                allow_pdb_name_mismatch=allow_pdb_name_mismatch,
            )

    messages = [line for line in output.getvalue().splitlines() if line.strip()]
    debug_files = []
    if keep_debug_files:
        stem = os.path.splitext(data_file)[0]
        debug_files = [
            os.path.abspath(f"{stem}.{name}")
            for name in ("pairs.txt", "bonds.txt", "angles.txt", "dihedrals.txt")
        ]

    response = _compact_conversion(
        conversion, root, data_file, parameter_file, provenance_file, debug_files
    )
    if detail == "full":
        response["provenance"] = conversion["record"]
        response["messages"] = messages
    return response


def _parse_lammps_energies(output: str) -> Dict[str, float]:
    wanted = {
        "E_bond": "bond", "E_angle": "angle", "E_dihed": "dihedral",
        "E_vdwl": "vdw", "E_coul": "coulomb_short",
        "E_long": "coulomb_long", "PotEng": "potential",
    }
    lines = output.splitlines()
    for index, line in enumerate(lines[:-1]):
        headers = line.split()
        if not all(field in headers for field in wanted):
            continue
        for candidate in lines[index + 1:index + 5]:
            values = candidate.split()
            if len(values) != len(headers):
                continue
            try:
                parsed = dict(zip(headers, map(float, values)))
            except ValueError:
                continue
            return {name: parsed[field] for field, name in wanted.items()}
    return {}


def _check_lammps_run0(provenance_file: str, timeout: float = 60.0,
                       overwrite: bool = False, output_root: Optional[str] = None,
                       lammps_executable: Optional[str] = None) -> Dict[str, Any]:
    """Run a bounded LAMMPS smoke test described by conversion provenance."""
    if timeout <= 0 or timeout > 600:
        raise ValueError("timeout must be greater than 0 and no more than 600 seconds.")
    root = output_root or os.environ.get("AMBER2LAMMPS_MCP_ROOT", os.getcwd())
    provenance_path = _resolve_output(provenance_file, root)
    if not os.path.isfile(provenance_path):
        raise ValueError(f"Provenance file not found: {provenance_file}")
    with open(provenance_path, encoding="utf-8") as stream:
        record = json.load(stream)

    data_record = record["outputs"]["lammps_data"]
    parameter_record = record["outputs"]["lammps_parameters"]
    data_file = _resolve_output(data_record["path"], root)
    parameter_file = _resolve_output(parameter_record["path"], root)
    for path, expected in (
        (data_file, data_record["sha256"]),
        (parameter_file, parameter_record["sha256"]),
    ):
        if not os.path.isfile(path) or sha256_file(path) != expected:
            raise ValueError(f"Conversion output is missing or differs from provenance: {path}")

    configured = lammps_executable or os.environ.get("AMBER2LAMMPS_LAMMPS")
    if not configured:
        raise ValueError(
            "No LAMMPS executable configured; set AMBER2LAMMPS_LAMMPS when starting the server."
        )
    executable = configured if os.path.isabs(configured) else shutil.which(configured)
    if not executable or not os.path.isfile(executable) or not os.access(executable, os.X_OK):
        raise ValueError(f"Configured LAMMPS executable is not executable: {configured}")

    base = os.path.basename(provenance_path)
    if base.endswith(".provenance.json"):
        base = base[:-len(".provenance.json")]
    else:
        base = os.path.splitext(base)[0]
    directory = os.path.dirname(provenance_path)
    input_file = _resolve_output(os.path.join(directory, base + ".run0.input"), root)
    log_file = _resolve_output(os.path.join(directory, base + ".run0.log"), root)
    with _VALIDATION_LOCK:
        existing = [path for path in (input_file, log_file) if os.path.exists(path)]
        if existing and not overwrite:
            raise FileExistsError(
                "Validation file(s) already exist; set overwrite=true to replace them: "
                + ", ".join(existing)
            )

        contract = record["lammps"]
        counts = record["conversion"]["counts"]
        commands = ["units real", f"atom_style {contract['atom_style']}", "boundary p p p"]
        if counts["bonds"]:
            commands.append("bond_style harmonic")
        if counts["angles"]:
            commands.append("angle_style harmonic")
        if counts["dihedrals"]:
            commands.append("dihedral_style fourier")
        if not contract["sets_pair_style"]:
            cutoff = contract["cutoff_angstrom"]
            commands.append(
                f"pair_style lj/cut/coul/{contract['coulomb']} {cutoff:g} {cutoff:g}"
            )
        commands.extend([
            f"read_data {os.path.relpath(data_file, directory)}",
            f"include {os.path.relpath(parameter_file, directory)}",
        ])
        if contract["requires_kspace"]:
            commands.append("kspace_style pppm 1.0e-6")
        commands.extend([
            "thermo_style custom step ebond eangle edihed eimp evdwl ecoul elong etail pe",
            "thermo 1",
            "run 0",
        ])
        with open(input_file, "w", encoding="utf-8") as stream:
            stream.write("\n".join(commands) + "\n")

        status = "failed"
        error = None
        try:
            result = subprocess.run(
                [executable, "-in", input_file, "-log", "none"], cwd=directory,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                timeout=timeout, check=False,
            )
            output = result.stdout
            status = "passed" if result.returncode == 0 else "failed"
            if status == "failed":
                error = next(
                    (line.strip() for line in output.splitlines() if "ERROR" in line),
                    f"LAMMPS exited with status {result.returncode}",
                )
        except subprocess.TimeoutExpired as exc:
            output = exc.stdout or ""
            if isinstance(output, bytes):
                output = output.decode(errors="replace")
            error = f"LAMMPS run exceeded {timeout:g} seconds"
        with open(log_file, "w", encoding="utf-8") as stream:
            stream.write(output)

        energies = _parse_lammps_energies(output)
        record["validation"] = {
            "status": status,
            "kind": "lammps_run0",
            "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
            "energies_kcal_mol": energies,
            "error": error,
            "log": _relative(log_file, root),
        }
        write_provenance(provenance_path, record)
    response = {
        "status": status,
        "energies_kcal_mol": energies,
        "error": error,
        "files": {
            "input": _relative(input_file, root),
            "log": _relative(log_file, root),
            "provenance": _relative(provenance_path, root),
        },
    }
    return response


def create_server():
    """Create the MCP server, or explain how to install its optional extra."""
    try:
        from mcp.server import MCPServer
        from mcp.server.mcpserver.exceptions import ToolError
    except ImportError as exc:
        raise RuntimeError(
            "MCP support requires Python 3.11+ and the optional dependency. "
            "Install it with: python -m pip install 'amber-to-lammps[mcp]'"
        ) from exc

    server = MCPServer(
        "AMBER2LAMMPS",
        instructions=(
            "Convert local AMBER topology and PDB files into LAMMPS data and "
            "parameter files. Confirm output paths with the user before replacing files."
        ),
    )

    @server.tool()
    def convert(
        pdb_file: str,
        topologies: List[str],
        counts: List[int],
        charges: List[float],
        data_file: str,
        parameter_file: str,
        buffer: float = 3.8,
        type_labels: bool = False,
        keep_debug_files: bool = False,
        c4_cutoff: float = 9.0,
        c4_coulomb: str = "long",
        provenance_file: Optional[str] = None,
        overwrite: bool = False,
        allow_pdb_name_mismatch: bool = False,
        detail: str = "summary",
    ) -> Dict[str, Any]:
        """Convert AMBER files to LAMMPS files on the local machine.

        Provide one count and target molecular charge for each topology.
        A valid orthogonal PDB CRYST1 cell overrides buffer. Outputs must be
        inside AMBER2LAMMPS_MCP_ROOT. Existing files require overwrite=true.
        """
        try:
            return convert_amber_system(
                pdb_file=pdb_file,
                topologies=topologies,
                counts=counts,
                charges=charges,
                data_file=data_file,
                parameter_file=parameter_file,
                buffer=buffer,
                type_labels=type_labels,
                keep_debug_files=keep_debug_files,
                c4_cutoff=c4_cutoff,
                c4_coulomb=c4_coulomb,
                provenance_file=provenance_file,
                overwrite=overwrite,
                allow_pdb_name_mismatch=allow_pdb_name_mismatch,
                detail=detail,
            )
        except Exception as exc:
            raise ToolError(str(exc)) from exc

    @server.tool()
    def check_lammps_run0(
        provenance_file: str,
        timeout: float = 60.0,
        overwrite: bool = False,
    ) -> Dict[str, Any]:
        """Run a bounded LAMMPS run-0 check and return compact energies/errors."""
        try:
            return _check_lammps_run0(
                provenance_file=provenance_file, timeout=timeout, overwrite=overwrite
            )
        except Exception as exc:
            raise ToolError(str(exc)) from exc

    return server


def main() -> None:
    """Run the local MCP server over stdio."""
    create_server().run()


if __name__ == "__main__":
    main()
