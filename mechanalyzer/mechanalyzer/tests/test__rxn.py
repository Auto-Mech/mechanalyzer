"""
test the code
"""

import mechanalyzer


INI_NAME = 'C4H10'
INI_ICH = 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3'

RADICALS = (
    ('H', 'InChI=1S/H'),
    ('OH', 'InChI=1S/HO/h1H'),
    ('CH3', 'InChI=1S/CH3/h1H3')
)
RAD_NAMES = tuple(rad[0] for rad in RADICALS)
RAD_ICHS = tuple(rad[1] for rad in RADICALS)

INI_ICHS = (INI_ICH,) + RAD_ICHS


def test__mech_build():
    """ test building the mechanism
    """

    # print('\nmechanism')
    mechanism, spc_dct = mechanalyzer.builder.rxn.build_mechanism(
        INI_ICHS)

    mech_str, spc_str = mechanalyzer.builder.rxn.mechanism_strs(
        mechanism, spc_dct)

    print('\n(mech_str)')
    print(mech_str)
    print('\n\n(spc_str)')
    print(spc_str)


if __name__ == '__main__':
    test__mech_build()
