""" test thermfit.cbh._basic
"""

import numpy
import automol.reac
import thermfit.cbh


# Species
C3H8_ICH = 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'
C3H7OH_ICH = 'InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3'
C3H7NH2_ICH = 'InChI=1S/C3H9N/c1-2-3-4/h2-4H2,1H3'
C3H7CL_ICH = 'InChI=1S/C3H7Cl/c1-2-3-4/h2-3H2,1H3'
CH3CH2SCH3_ICH = 'InChI=1S/C3H8S/c1-3-4-2/h3H2,1-2H3'
HOCH2CH2SCH3_ICH = 'InChI=1S/C3H8OS/c1-5-3-2-4/h4H,2-3H2,1H3'

# TS ZRXN OBJECTS
HYD_ABS_ZRXN_STR = """
reaction class: hydrogen abstraction
forward TS atoms:
  1: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: X, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  2-4: {order: 1, stereo_parity: null}
  2-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  4-7: {order: 0.9, stereo_parity: null}
  4-8: {order: 1, stereo_parity: null}
  4-9: {order: 1, stereo_parity: null}
  7-10: {order: 0, stereo_parity: null}
  7-11: {order: 0.1, stereo_parity: null}
  11-12: {order: 1, stereo_parity: null}
  11-13: {order: 1, stereo_parity: null}
  11-14: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
- [11, 12, 13, 14]
backward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 0.9, stereo_parity: null}
  5-6: {order: 0.1, stereo_parity: null}
  6-7: {order: 1, stereo_parity: null}
  6-9: {order: 1, stereo_parity: null}
  6-10: {order: 1, stereo_parity: null}
  7-8: {order: 1, stereo_parity: null}
  7-11: {order: 1, stereo_parity: null}
  7-12: {order: 1, stereo_parity: null}
  8-13: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5]
- [6, 7, 8, 9, 10, 11, 12, 13]
"""

# Set scheme variable
SCHEME = 'basic'


def test__species():
    """ test thermfit.cbh._spc
    """

    c3h8_basis = thermfit.cbh.species_basis(C3H8_ICH, SCHEME)
    c3h7oh_basis = thermfit.cbh.species_basis(C3H7OH_ICH, SCHEME)
    c3h7nh2_basis = thermfit.cbh.species_basis(C3H7NH2_ICH, SCHEME)
    c3h7cl_basis = thermfit.cbh.species_basis(C3H7CL_ICH, SCHEME)
    ch3ch2sch3_basis = thermfit.cbh.species_basis(CH3CH2SCH3_ICH, SCHEME)
    hoch2ch2sch3_basis = thermfit.cbh.species_basis(HOCH2CH2SCH3_ICH, SCHEME)

    assert c3h8_basis[0] == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4')
    assert numpy.allclose(c3h8_basis[1],
                          numpy.array([-2.,  3.]))

    assert c3h7oh_basis[0] == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/H2O/h1H2')
    assert numpy.allclose(c3h7oh_basis[1],
                          numpy.array([-3.,  3.,  1.]))

    assert c3h7nh2_basis[0] == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/H3N/h1H3')
    assert numpy.allclose(c3h7nh2_basis[1],
                          numpy.array([-3.,  3.,  1.]))

    assert c3h7cl_basis[0] == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/ClH/h1H')
    assert numpy.allclose(c3h7cl_basis[1],
                          numpy.array([-3.,  3.,  1.]))

    assert ch3ch2sch3_basis[0] == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/O2S/c1-3-2',
        'InChI=1S/H2O/h1H2')
    assert numpy.allclose(ch3ch2sch3_basis[1],
                          numpy.array([0.,  3.,  1., -2.]))

    assert hoch2ch2sch3_basis[0] == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/H2O/h1H2',
        'InChI=1S/O2S/c1-3-2')
    assert numpy.allclose(hoch2ch2sch3_basis[1],
                          numpy.array([-1.,  3., -1.,  1.]))


def test__transition_states():
    """ test thermfit.cbh._ts
    """

    habs_zrxn = automol.reac.from_string(HYD_ABS_ZRXN_STR)
    habs_basis = thermfit.cbh.ts_basis(habs_zrxn, SCHEME)

    assert habs_basis[0] == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/H2O/h1H2')
    assert numpy.allclose(habs_basis[1],
                          numpy.array([-2.5, 3.0, 1.0]))


if __name__ == '__main__':
    test__species()
    test__transition_states()
