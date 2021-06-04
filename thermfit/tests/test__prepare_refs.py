""" Test prepare refs
"""

import thermfit


# Overall species dict
SPC_DCT = {
    'C2H6': {
        'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3',
        'mult': 1,
        'charge': 0},
    'C2H5': {
        'inchi': 'InChI=1S/C2H5/c1-2/h1H2,2H3',
        'mult': 2,
        'charge': 0},
    'H': {
        'inchi': 'InChI=1S/H',
        'mult': 2,
        'charge': 0},
    'H2': {
        'inchi': 'InChI=1S/H2/h1H',
        'mult': 1,
        'charge': 0},
    'CH4': {
        'inchi': 'InChI=1S/CH4/h1H4',
        'mult': 1,
        'charge': 0},

}

# Info for species refs
SPC_NAMES = ('C2H6',)

# Scheme information
SCHEME = 'basic'


def test__():
    """ test thermfit.._basis.prepare_refs
    """

    dct1, dct2 = thermfit.prepare_refs(
        SCHEME, SPC_DCT, SPC_NAMES,
        repeats=False, parallel=False, zrxn=None)

    print(dct1)
    print(dct2)


if __name__ == '__main__':
    test__()
