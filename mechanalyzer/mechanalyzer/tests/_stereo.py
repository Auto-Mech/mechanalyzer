""" test the stereochem
"""

import os
import chemkin_io
import mechanalyzer


# Get mechanism information
def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
MECH_NAME = 'ste.txt'
CSV_NAME = 'ste.csv'

MECH_STR = _read_file(
    os.path.join(DATA_PATH, MECH_NAME))
CSV_STR = _read_file(
    os.path.join(DATA_PATH, CSV_NAME))

# Get the species
SPC_DCT = mechanalyzer.parser.spc.build_spc_dct(CSV_STR, 'csv')

# Get the reactions
REACTION_BLOCK = chemkin_io.parser.mechanism.reaction_block(
    MECH_STR)
REACTION_DCT = chemkin_io.parser.reaction.data_dct(
    REACTION_BLOCK)


def test__():
    """ test
    """

    STE = mechanalyzer.parser.expand_mech_stereo(
        REACTION_DCT, SPC_DCT)
    print('STE', STE)
    # STE_MECH_DCT, STE_SPC_DCT = mechanalyzer.parser.expand_mech_stereo(
    #     REACTION_DCT, SPC_DCT)


if __name__ == '__main__':
    test__()
