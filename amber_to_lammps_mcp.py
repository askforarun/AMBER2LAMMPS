"""Local MCP server for AMBER2LAMMPS.

The MCP dependency is imported only when the server starts, keeping the core
converter and its Python 3.9+ installation independent of MCP.
"""

import contextlib
import io
import os
import threading
from typing import Any, Dict, List, Optional

from amber_to_lammps import amber2lammps


_CONVERSION_LOCK = threading.Lock()


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

    return {
        "status": "success",
        "data_file": os.path.abspath(data_file),
        "parameter_file": os.path.abspath(parameter_file),
        "debug_files": debug_files,
        "messages": messages,
        "summary": conversion["summary"],
    }


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
            )
        except Exception as exc:
            raise ToolError(str(exc)) from exc

    return server


def main() -> None:
    """Run the local MCP server over stdio."""
    create_server().run()


if __name__ == "__main__":
    main()
