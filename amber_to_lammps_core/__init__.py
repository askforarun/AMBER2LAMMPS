"""Internal, testable scientific core for AMBER2LAMMPS."""

from .model import SpecialBondScaling
from .topology import (
    UnsupportedScalingError,
    dihedral_phase_degrees,
    distinct_14_scaling_profiles,
    resolve_system_14_scaling,
    topology_14_scaling,
)

__all__ = [
    "SpecialBondScaling",
    "UnsupportedScalingError",
    "dihedral_phase_degrees",
    "distinct_14_scaling_profiles",
    "resolve_system_14_scaling",
    "topology_14_scaling",
]
