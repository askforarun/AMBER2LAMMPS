"""Force-field-neutral intermediate model objects."""

from dataclasses import dataclass


@dataclass(frozen=True)
class SpecialBondScaling:
    """AMBER divisors and their equivalent LAMMPS 1-4 weights."""

    scee: float
    scnb: float

    @property
    def coulomb_weight(self):
        return 1.0 / self.scee

    @property
    def lj_weight(self):
        return 1.0 / self.scnb
