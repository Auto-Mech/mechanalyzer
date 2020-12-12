""" test chemkin_io.calculator.thermo
"""

from __future__ import unicode_literals
from builtins import open
import os
import numpy as np
import chemkin_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


# Set paths
PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
SYNGAS_MECH_NAME = 'syngas_mechanism.txt'
SYNGAS_CSV_NAME = 'syngas_species.csv'

# Read mechanism and csvfiles
SYNGAS_MECH_STR = _read_file(
    os.path.join(DATA_PATH, SYNGAS_MECH_NAME))
SYNGAS_CSV_STR = _read_file(
    os.path.join(DATA_PATH, SYNGAS_CSV_NAME))

# Read species blocks
SYNGAS_THERMO_BLOCK = chemkin_io.parser.mechanism.thermo_block(
    SYNGAS_MECH_STR)
SYNGAS_BLOCK_STRS = chemkin_io.parser.thermo.data_strings(
    SYNGAS_THERMO_BLOCK)

# Set polynomial and temps for comparison tests
SPECIES_IDX = 10
SPECIES_POLYNOMIAL = SYNGAS_BLOCK_STRS[SPECIES_IDX]
TEMP1, TEMP2 = 500.0, 1000.0


def test__mechanism():
    """ test chemkin_io.calculator.thermo.mechanism
    """
    therm_dct = chemkin_io.calculator.thermo.mechanism(
        SYNGAS_THERMO_BLOCK, [TEMP1, TEMP2])
    print(therm_dct)


def test__enthalpy():
    """ test chemkin_io.calculator.thermo.enthalpy
    """
    ht1 = chemkin_io.calculator.thermo.enthalpy(
        SPECIES_POLYNOMIAL, TEMP1)
    ht2 = chemkin_io.calculator.thermo.enthalpy(
        SPECIES_POLYNOMIAL, TEMP2)

    ref_ht1 = -30.240360335437266
    ref_ht2 = -23.421760051163762

    assert np.isclose(ht1, ref_ht1)
    assert np.isclose(ht2, ref_ht2)


def test__entropy():
    """ test chemkin_io.calculator.thermo.entropy
    """
    st1 = chemkin_io.calculator.thermo.entropy(
        SPECIES_POLYNOMIAL, TEMP1)
    st2 = chemkin_io.calculator.thermo.entropy(
        SPECIES_POLYNOMIAL, TEMP2)

    ref_st1 = 0.06174090267135737
    ref_st2 = 0.07107861508623411

    assert np.isclose(st1, ref_st1)
    assert np.isclose(st2, ref_st2)


def test__gibbs():
    """ test chemkin_io.calculator.thermo.gibbs
    """
    gt1 = chemkin_io.calculator.thermo.gibbs(
        SPECIES_POLYNOMIAL, TEMP1)
    gt2 = chemkin_io.calculator.thermo.gibbs(
        SPECIES_POLYNOMIAL, TEMP2)

    ref_gt1 = -61.110811671115954
    ref_gt2 = -94.50037513739787

    assert np.isclose(gt1, ref_gt1)
    assert np.isclose(gt2, ref_gt2)


def test__heat_capacity():
    """ test chemkin_io.calculator.thermo.heat_capacity
    """
    cp1 = chemkin_io.calculator.thermo.heat_capacity(
        SPECIES_POLYNOMIAL, TEMP1)
    cp2 = chemkin_io.calculator.thermo.heat_capacity(
        SPECIES_POLYNOMIAL, TEMP2)

    ref_cp1 = 0.011992997783769053
    ref_cp2 = 0.01494647663556653

    assert np.isclose(cp1, ref_cp1)
    assert np.isclose(cp2, ref_cp2)


if __name__ == '__main__':
    test__mechanism()
    test__enthalpy()
    test__entropy()
    test__gibbs()
    test__heat_capacity()
