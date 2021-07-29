""" species dict function tests
"""

import mechanalyzer
from mechanalyzer import par


# Species info
INCHI1 = 'InChI=1S/C2H5/c1-2/h1H2,2H3'
CHARGE1 = 0
MULT1 = 2

INCHI2 = 'InChI=1S/H'
CHARGE2 = 0
MULT2 = 2

INCHI3 = 'InChI=1S/C2H6/c1-2/h1-2H3'
CHARGE3 = 0
MULT3 = 1

INCHI4 = 'InChI=1S/C2H4/c1-2/h1-2H2'
CHARGE4 = 1  # Fake charge
MULT4 = 3  # Fake triplet for testing

INCHI5 = 'InChI=1S/H2/h1H'
CHARGE5 = 0
MULT5 = 1

INCHI6 = 'FakeIch6'
CHARGE6 = 0
MULT6 = 2

INCHI7 = 'FakeIch7'
CHARGE7 = 1
MULT7 = 2


SPC_DCT = {
    'C2H5': {
        par.SPC.INCHI: INCHI1,
        par.SPC.CHARGE: CHARGE1,
        par.SPC.MULT: MULT1
    },
    'H': {
        par.SPC.INCHI: INCHI2,
        par.SPC.CHARGE: CHARGE2,
        par.SPC.MULT: MULT2
    },
    'C2H6': {
        par.SPC.INCHI: INCHI3,
        par.SPC.CHARGE: CHARGE3,
        par.SPC.MULT: MULT3
    },
    'C2H4': {
        par.SPC.INCHI: INCHI4,
        par.SPC.CHARGE: CHARGE4,
        par.SPC.MULT: MULT4
    },
    'H2': {
        par.SPC.INCHI: INCHI5,
        par.SPC.CHARGE: CHARGE5,
        par.SPC.MULT: MULT5
    },
}

HEADERS = ('smiles', 'inchi', 'mult', 'charge')


def test__csv_str():
    """ test mechanalyzer.parser.spc
    """

    csv_str = mechanalyzer.parser.spc.dct_to_csv_str(SPC_DCT, HEADERS)
    print(csv_str)


def test__fml_dct():
    """ test mechanalyzer.parser.spc
    """

    fml_count_dct = mechanalyzer.parser.spc.build_fml_count_dct(
        SPC_DCT)
    print(fml_count_dct)


if __name__ == '__main__':
    test__csv_str()
    test__fml_dct()
