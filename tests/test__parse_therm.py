""" test chemkin_io.parser.thermo
"""

from __future__ import unicode_literals
from builtins import open
import os
import numpy
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

# Set polynomial for comparison tests
SPECIES_IDX = 10
SPECIES_POLYNOMIAL = SYNGAS_BLOCK_STRS[SPECIES_IDX]


def test__data_strings():
    """ test chemkin_io.parser.thermo.data_strings
    """
    thm_strs = chemkin_io.parser.thermo.data_strings(
        SYNGAS_THERMO_BLOCK)
    assert len(thm_strs) == 19


def test__data_block():
    """ test chemkin_io.parser.thermo.data_block
    """
    blocks = chemkin_io.parser.thermo.data_block(
        SYNGAS_THERMO_BLOCK)
    spc_block = blocks[SPECIES_IDX]

    ref_block = (
        'H2O2(11)',
        (200.0, 6000.0, 1000.0),
        (4.31515, -0.000847391, 1.76404e-05, -2.26763e-08,
         9.0895e-12, -17706.7, 3.27373),
        (4.57977, 0.00405326, -1.29845e-06, 1.98211e-10,
         -1.13969e-14, -18007.2, 0.664971)
    )

    assert len(spc_block) == len(ref_block)
    assert spc_block[0] == ref_block[0]
    assert numpy.allclose(spc_block[1], ref_block[1])
    assert numpy.allclose(spc_block[2], ref_block[2])
    assert numpy.allclose(spc_block[3], ref_block[3])


def test__species_name():
    """ test chemkin_io.parser.thermo.species_name
    """
    name = chemkin_io.parser.thermo.species_name(
        SPECIES_POLYNOMIAL)
    ref_name = 'H2O2(11)'
    assert name == ref_name


def test__temp_common_default():
    """ test chemkin_io.parser.thermo.temp_common_default
    """
    temp_common = chemkin_io.parser.thermo.temp_common_default(
        SYNGAS_THERMO_BLOCK)
    ref_temp_common = 1000.0
    assert numpy.isclose(temp_common, ref_temp_common)


def test__temperatures():
    """ test chemkin_io.parser.thermo.temperatures
    """
    temps = chemkin_io.parser.thermo.temperatures(
        SPECIES_POLYNOMIAL)
    ref_temps = (200.0, 6000.0, 1000.0)
    assert numpy.allclose(temps, ref_temps)


def test__low_coefficients():
    """ test chemkin_io.parser.thermo.low_coefficients
    """
    low_cfts = chemkin_io.parser.thermo.low_coefficients(
        SPECIES_POLYNOMIAL)
    ref_low_cfts = (4.31515, -0.000847391, 1.76404e-05, -2.26763e-08,
                    9.0895e-12, -17706.7, 3.27373)
    assert numpy.allclose(low_cfts, ref_low_cfts)


def test__high_coefficients():
    """ test chemkin_io.parser.thermo.high_coefficients
    """
    high_cfts = chemkin_io.parser.thermo.high_coefficients(
        SPECIES_POLYNOMIAL)
    ref_high_cfts = (4.57977, 0.00405326, -1.29845e-06, 1.98211e-10,
                     -1.13969e-14, -18007.2, 0.664971)
    assert numpy.allclose(high_cfts, ref_high_cfts)


if __name__ == '__main__':
    test__data_strings()
    test__data_block()
    test__temp_common_default()
    test__species_name()
    test__temperatures()
    test__low_coefficients()
    test__high_coefficients()
