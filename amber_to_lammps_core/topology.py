"""Analysis of topology-level AMBER scaling semantics."""

from .model import SpecialBondScaling


class UnsupportedScalingError(ValueError):
    """Raised when AMBER 1-4 scaling cannot be represented exactly."""


def topology_14_scaling(parm, topology_name="topology"):
    """Return the unique active proper-dihedral SCEE/SCNB profile."""
    profiles = set()
    for dihedral in parm.dihedrals:
        if dihedral.improper or dihedral.ignore_end:
            continue
        terms = dihedral.type if isinstance(dihedral.type, (list, tuple)) else [dihedral.type]
        for term in terms:
            if term is None:
                continue
            scee = float(term.scee)
            scnb = float(term.scnb)
            if scee <= 0.0 or scnb <= 0.0:
                raise UnsupportedScalingError(
                    f"{topology_name} has non-positive SCEE/SCNB values on an active "
                    f"proper dihedral: SCEE={scee:g}, SCNB={scnb:g}."
                )
            profiles.add((round(scee, 10), round(scnb, 10)))
    if not profiles:
        return None
    if len(profiles) != 1:
        details = ", ".join(
            f"SCEE={scee:g}/SCNB={scnb:g}" for scee, scnb in sorted(profiles)
        )
        raise UnsupportedScalingError(
            f"{topology_name} uses multiple active 1-4 scaling profiles ({details}); "
            "one global LAMMPS special_bonds setting cannot represent them exactly."
        )
    scee, scnb = next(iter(profiles))
    return SpecialBondScaling(scee=scee, scnb=scnb)


def resolve_system_14_scaling(profiles, topologies):
    """Return the common active profile, or ``None`` for mixed/no profiles."""
    used = {}
    for topology, profile in zip(topologies, profiles):
        if profile is not None:
            used.setdefault(profile, []).append(topology)
    return next(iter(used), None) if len(used) == 1 else None


def distinct_14_scaling_profiles(profiles):
    """Return active scaling profiles in stable first-seen order."""
    return list(dict.fromkeys(profile for profile in profiles if profile is not None))


def dihedral_phase_degrees(term):
    """Return a ParmEd dihedral phase, whose public unit is degrees."""
    return float(term.phase)
