""" Script for calculating CBH values
"""

import automol.reac
import thermfit.cbh


# Species ZMA
CH3CH2OH_ICH = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
CH3CH2OH_ZMA = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('C', (0, None, None), ('R1', None, None),
     (2.8621866421132123, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.0691950317120837, 1.9320905931404335, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.0665277363750003, 1.9353637583977303, 2.108686322109069)),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4'),
     (2.069152234029259, 1.930351996269131, 4.218481783495319)),
    ('O', (1, 0, 2), ('R5', 'A5', 'D5'),
     (2.683747067528887, 1.9194991329920528, 1.0492575691901376)),
    ('H', (1, 0, 5), ('R6', 'A6', 'D6'),
     (2.067751813667005, 1.9377163905267163, 4.182292202036946)),
    ('H', (1, 0, 5), ('R7', 'A7', 'D7'),
     (2.0671092193398275, 1.928341045304867, 2.082936263972801)),
    ('H', (5, 1, 0), ('R8', 'A8', 'D8'),
     (1.8376095698733521, 1.8770195234414855, 5.238498010242486)))
# Example for high-level CBH schemes
# Radical Example

# TS ZRXN OBJECTS
ZRXN_STR = """
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


def test__species():
    """ test
    """
    basis = thermfit.cbh.species_cbh_basis(CH3CH2OH_ICH, 'cbh0')
    print(basis)

    assert basis == (
        ('InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2', 'InChI=1S/H2/h1H'),
        (2, 1, -2.0)
    )


def test__transition_state():
    """ test
    """
    zrxn = automol.reac.from_string(ZRXN_STR)
    basis = thermfit.cbh.ts_cbh_basis(zrxn, 'cbh0')
    print(basis)

    assert basis == (
        ['InChI=1S/H2O/h1H2', 'InChI=1S/CH4/h1H4',
         (('InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3'),
          ('InChI=1S/CH4/h1H4', 'InChI=1S/CH3/h1H3')),
         'InChI=1S/H2/h1H'],
        [1.0, 1.0, 1.0, -2.0]
    )


if __name__ == '__main__':
    test__species()
    test__transition_state()
