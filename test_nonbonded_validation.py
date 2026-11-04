import unittest

from amber_to_lammps import (
    _ab_to_epsilon_sigma,
    find_pair_specific_nonbonded_terms,
)


class FakeParm:
    def __init__(self, acoef, bcoef, nb_index, ccoef=None, ntypes=2):
        self._ntypes = ntypes
        self.parm_data = {
            "LENNARD_JONES_ACOEF": acoef,
            "LENNARD_JONES_BCOEF": bcoef,
            "NONBONDED_PARM_INDEX": nb_index,
        }
        if ccoef is not None:
            self.parm_data["LENNARD_JONES_CCOEF"] = ccoef

    def ptr(self, name):
        if name.lower() != "ntypes":
            raise KeyError(name)
        return self._ntypes


def mixed_two_type_parm(ccoef=None):
    epsilon_1, epsilon_2 = 0.2, 0.1
    rmin_1, rmin_2 = 3.0, 4.0
    epsilon_12 = (epsilon_1 * epsilon_2) ** 0.5
    rmin_12 = 0.5 * (rmin_1 + rmin_2)
    acoef = [
        epsilon_1 * rmin_1**12,
        epsilon_12 * rmin_12**12,
        epsilon_2 * rmin_2**12,
    ]
    bcoef = [
        2.0 * epsilon_1 * rmin_1**6,
        2.0 * epsilon_12 * rmin_12**6,
        2.0 * epsilon_2 * rmin_2**6,
    ]
    return FakeParm(acoef, bcoef, [1, 2, 2, 3], ccoef)


class NonbondedValidationTests(unittest.TestCase):
    def test_standard_mixed_lj_is_supported(self):
        self.assertEqual(find_pair_specific_nonbonded_terms(mixed_two_type_parm()), [])

    def test_zero_c4_array_is_supported(self):
        self.assertEqual(find_pair_specific_nonbonded_terms(mixed_two_type_parm([0.0, 0.0, 0.0])), [])

    def test_nonzero_c4_is_detected(self):
        parm = mixed_two_type_parm([0.0, 12.5, 0.0])
        self.assertIn("12-6-4 C4", find_pair_specific_nonbonded_terms(parm)[0])

    def test_off_diagonal_override_is_detected(self):
        parm = mixed_two_type_parm()
        parm.parm_data["LENNARD_JONES_ACOEF"][1] *= 1.1
        problems = find_pair_specific_nonbonded_terms(parm)
        self.assertEqual(len(problems), 1)
        self.assertIn("off-diagonal LJ pair override", problems[0])

    def test_raw_ab_conversion(self):
        epsilon, sigma = _ab_to_epsilon_sigma(0.2 * 3.0**12, 2.0 * 0.2 * 3.0**6)
        self.assertAlmostEqual(epsilon, 0.2)
        self.assertAlmostEqual(sigma, 3.0 / 2.0 ** (1.0 / 6.0))


if __name__ == "__main__":
    unittest.main()
