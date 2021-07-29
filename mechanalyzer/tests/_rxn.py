""" test automol.reac
"""

import automol
import mechanalyzer


# Define mechanism start and series for reaction generation
MECH_SPC_DCT = {
    'C4H7O2(1)': {
        'smiles': 'O=CCCC[O]',
        'inchi': 'InChI=1S/C4H7O2/c5-3-1-2-4-6/h3H,1-2,4H2',
        'charge': 0,
        'mult': 2
    }
}
MECH_RXN_DCT = {}

RSERIES = (
    ('all', (automol.par.ReactionClass.Typ.HYDROGEN_MIGRATION,)),
)

RADICALS = (
    ('H', 'InChI=1S/H'),
    ('OH', 'InChI=1S/HO/h1H'),
    ('CH3', 'InChI=1S/CH3/h1H3')
)
RAD_NAMES = tuple(rad[0] for rad in RADICALS)
RAD_ICHS = tuple(rad[1] for rad in RADICALS)


def test__mech_build():
    """ test building the mechanism
    """

    spc_dct, rxn_dct = mechanalyzer.builder.rxn.build_mechanism(
        MECH_SPC_DCT, MECH_RXN_DCT, rxn_series=RSERIES)
    print('species')
    for name, dct in spc_dct.items():
        print(name)
        print(dct)
        print('\n')
    print('\nreaction')
    for name, dct in rxn_dct.items():
        print(name)
        print(dct)
        print('\n')


if __name__ == '__main__':
    test__mech_build()
