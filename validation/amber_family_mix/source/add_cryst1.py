#!/usr/bin/env python3
"""Add the explicit periodic cell used by both validation engines."""

from pathlib import Path


path = Path(__file__).resolve().parent / "ethanol_glucose.pdb"
lines = path.read_text().splitlines()
lines = [line for line in lines if not line.startswith("CRYST1")]
insert_at = next(i for i, line in enumerate(lines) if line.startswith(("ATOM", "HETATM")))
lines.insert(
    insert_at,
    "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1",
)
path.write_text("\n".join(lines) + "\n")
print(f"Added a 30 Angstrom orthorhombic CRYST1 cell to {path.name}")
