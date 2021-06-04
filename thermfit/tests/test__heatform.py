""" test thermfit.heatform
"""

import numpy
import thermfit.heatform


# Formula dicts for various species
C3H7OH_H0 = -200.0
C3H7OH_BASIS = ('InChI=1S/H2/h1H', 'InChI=1S/CH4/h1H4', 'InChI=1S/H2O/h1H2')
C3H7OH_BASIS_H0 = (-0.48, -100.0, -230.0)
C3H7OH_CFFS = numpy.array([-3.,  3.,  1.])
REF_SET1 = 'ANL0'


def test__hform():
    """ test thermfit.heatform.calc_hform_0k
    """

    # test all ref sets
    hf0k = thermfit.heatform.calc_hform_0k(
       C3H7OH_H0, C3H7OH_BASIS_H0, C3H7OH_BASIS, C3H7OH_CFFS, REF_SET1)
    print(hf0k)


if __name__ == '__main__':
    test__hform()
