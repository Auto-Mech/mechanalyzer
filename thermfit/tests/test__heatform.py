""" test thermfit.heatform
"""

import numpy
import thermfit.heatform


# Formula dicts for various species
C3H7OH_ICH = 'InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3'
C3H7OH_H0 = -200.0
C3H7OH_BASIS = ('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2')
C3H7OH_BASIS_H0 = (-0.48, -100.0, -230.0)
C3H7OH_CFFS = numpy.array([-3., 3., 1.])

# Transition state names
TSNAME = 'CH4+O=CH3+OH'
RXN_ICHS = (('InChI=1S/CH4/h1H4', 'InChI=1S/HO/h1H'),
            ('InChI=1S/CH3/h1H3', 'InChI=1S/H2O/h1H2'))

# Temps for database
TEMP1 = 0
TEMP2 = 298

# Database sets to read
REF_SET1 = 'ANL0'
REF_SET2 = 'ATcT'


def test__hform():
    """ test thermfit.heatform.calc_hform_0k
    """

    # test all ref sets
    hf0k = thermfit.heatform.calc_hform_0k(
        C3H7OH_H0, C3H7OH_BASIS_H0, C3H7OH_BASIS, C3H7OH_CFFS, REF_SET1)
    print(hf0k)

    # test for reaction


def test__ref_enthalpy():
    """ test thermfit.heatform.reference_enthalpy
    """

    # Species
    hf0k = thermfit.heatform.reference_enthalpy(
        C3H7OH_H0, REF_SET1, TEMP1, rxn=False)
    print(hf0k)

    # Transition State
    db_rxn_ich = thermfit.heatform.format_reaction_inchi(RXN_ICHS)
    hf0k = thermfit.heatform.reference_enthalpy(
        db_rxn_ich, REF_SET1, TEMP1, rxn=True)
    print(hf0k)


if __name__ == '__main__':
    test__hform()
    test__ref_enthalpy()
