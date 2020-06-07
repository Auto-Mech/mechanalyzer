""" test chemkin_io.parser.mechanism
"""

from __future__ import unicode_literals
from builtins import open
import os
import chemkin_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


# Set paths
PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
HEPTANE_MECH_NAME = 'heptane_mechanism.txt'
SYNGAS_MECH_NAME = 'syngas_mechanism.txt'

# Read mechanism files
HEPTANE_MECH_STR = _read_file(
    os.path.join(DATA_PATH, HEPTANE_MECH_NAME))
SYNGAS_MECH_STR = _read_file(
    os.path.join(DATA_PATH, SYNGAS_MECH_NAME))

# Read species blocks
HEPTANE_SPECIES_BLOCK = chemkin_io.parser.mechanism.species_block(
    HEPTANE_MECH_STR)
SYNGAS_SPECIES_BLOCK = chemkin_io.parser.mechanism.species_block(
    SYNGAS_MECH_STR)


def test__names():
    """ test chemkin_io.parser.species.names
    """
    spc_names1 = chemkin_io.parser.species.names(
        HEPTANE_SPECIES_BLOCK)
    spc_names2 = chemkin_io.parser.species.names(
        SYNGAS_SPECIES_BLOCK)
    assert len(spc_names1) == 1268
    assert len(spc_names2) == 19


if __name__ == '__main__':
    test__names()
