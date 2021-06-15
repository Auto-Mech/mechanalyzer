""" Test prepare refs
"""

import automol.reac
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
    'ts_1_1_1': {
        'inchi': '',
        'mult': 2,
        'charge': 0},

}

# Info for species
SPC_NAMES = ('C2H6',)
TS_NAMES = ('ts_1_1_1',)

# Info for TS
HABS_ZRXN_STR = """
reaction class: hydrogen abstraction
forward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  15: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  16: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  17: {symbol: X, implicit_hydrogen_valence: 0, stereo_parity: null}
  18: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  19: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  20: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  21: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  22: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  23: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  24: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  25: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  26: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  27: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  6-9: {order: 1, stereo_parity: null}
  6-10: {order: 1, stereo_parity: null}
  6-11: {order: 1, stereo_parity: null}
  7-12: {order: 1, stereo_parity: null}
  7-13: {order: 1, stereo_parity: null}
  7-14: {order: 1, stereo_parity: null}
  12-15: {order: 1, stereo_parity: null}
  15-16: {order: 0.9, stereo_parity: null}
  16-17: {order: 0, stereo_parity: null}
  16-18: {order: 0.1, stereo_parity: null}
  18-19: {order: 1, stereo_parity: null}
  18-20: {order: 1, stereo_parity: null}
  18-21: {order: 1, stereo_parity: null}
  19-22: {order: 1, stereo_parity: null}
  19-23: {order: 1, stereo_parity: null}
  19-24: {order: 1, stereo_parity: null}
  20-25: {order: 1, stereo_parity: null}
  20-26: {order: 1, stereo_parity: null}
  20-27: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
- [18, 19, 20, 21, 22, 23, 24, 25, 26, 27]
backward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  15: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  16: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  17: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  18: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  19: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  20: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  21: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  22: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  23: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  24: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  25: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  26: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  1-6: {order: 1, stereo_parity: null}
  2-3: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  2-9: {order: 1, stereo_parity: null}
  3-10: {order: 1, stereo_parity: null}
  3-11: {order: 0.9, stereo_parity: null}
  11-16: {order: 0.1, stereo_parity: null}
  12-15: {order: 1, stereo_parity: null}
  12-18: {order: 1, stereo_parity: null}
  12-19: {order: 1, stereo_parity: null}
  12-20: {order: 1, stereo_parity: null}
  13-15: {order: 1, stereo_parity: null}
  13-21: {order: 1, stereo_parity: null}
  13-22: {order: 1, stereo_parity: null}
  13-23: {order: 1, stereo_parity: null}
  14-15: {order: 1, stereo_parity: null}
  14-17: {order: 1, stereo_parity: null}
  14-24: {order: 1, stereo_parity: null}
  14-25: {order: 1, stereo_parity: null}
  15-26: {order: 1, stereo_parity: null}
  16-17: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
- [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]
"""
HABS_ZRXN = automol.reac.from_string(HABS_ZRXN_STR)


# Scheme information
SCHEME = 'basic'


def test__species():
    """ test thermfit.._basis.prepare_refs
    """

    dct1, dct2 = thermfit.prepare_refs(
        SCHEME, SPC_DCT, SPC_NAMES,
        repeats=False, parallel=False, zrxn=None)

    print(dct1)
    print(dct2)


def test__transition_state():
    """ test thermfit.._basis.prepare_refs
    """

    dct1, dct2 = thermfit.prepare_refs(
        SCHEME, SPC_DCT, TS_NAMES,
        repeats=False, parallel=False, zrxn=HABS_ZRXN)

    print(dct1)
    print(dct2)

    dct1, dct2 = thermfit.prepare_refs(
        'cbh0', SPC_DCT, TS_NAMES,
        repeats=False, parallel=False, zrxn=HABS_ZRXN)

    print(dct1)
    print(dct2)

    dct1, dct2 = thermfit.prepare_refs(
        'cbh1', SPC_DCT, TS_NAMES,
        repeats=False, parallel=False, zrxn=HABS_ZRXN)

    print(dct1)
    print(dct2)

    dct1, dct2 = thermfit.prepare_refs(
        'cbh1_0', SPC_DCT, TS_NAMES,
        repeats=False, parallel=False, zrxn=HABS_ZRXN)

    print(dct1)
    print(dct2)


if __name__ == '__main__':
    test__species()
    test__transition_state()
