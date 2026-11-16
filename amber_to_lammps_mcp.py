"""Local MCP server for AMBER2LAMMPS.

The MCP dependency is imported only when the server starts, keeping the core
converter and its Python 3.9+ installation independent of MCP.
"""

import contextlib
import io
import os
from typing import Any, Dict, List, Optional

from amber_to_lammps import amber2lammps


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
) -> Dict[str, Any]:
    """Convert AMBER topology/PDB inputs and return a structured result.

    Paths may be absolute or relative to the MCP server process. Existing
    output files with the requested names are replaced by the converter.
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

    if provenance_file is None:
        provenance_file = os.path.splitext(data_file)[0] + ".provenance.json"

    output = io.StringIO()
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
    ) -> Dict[str, Any]:
        """Convert AMBER files to LAMMPS files on the local machine.

        Provide one count and target molecular charge for each topology.
        A valid orthogonal PDB CRYST1 cell overrides buffer. Existing output
        files at data_file and parameter_file are replaced.
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
            )
        except Exception as exc:
            raise ToolError(str(exc)) from exc

    return server


def main() -> None:
    """Run the local MCP server over stdio."""
    create_server().run()


if __name__ == "__main__":
    main()
