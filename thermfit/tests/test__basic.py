""" test thermfit.cbh._basic
"""

import numpy
import thermfit.cbh


# Formula dicts for various species
C3H8_ICH = 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'
C3H7OH_ICH = 'InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3'
C3H7NH2_ICH = 'InChI=1S/C3H9N/c1-2-3-4/h2-4H2,1H3'
C3H7CL_ICH = 'InChI=1S/C3H7Cl/c1-2-3-4/h2-3H2,1H3'
CH3CH2SCH3_ICH = 'InChI=1S/C3H8S/c1-3-4-2/h3H2,1-2H3'
HOCH2CH2SCH3_ICH = 'InChI=1S/C3H8OS/c1-5-3-2-4/h4H,2-3H2,1H3'


def test__basis():
    """ test thermfit.cbh._basic.basic_basis
    """

    c3h8_basis = thermfit.cbh.basic_basis(C3H8_ICH)
    c3h7oh_basis = thermfit.cbh.basic_basis(C3H7OH_ICH)
    c3h7nh2_basis = thermfit.cbh.basic_basis(C3H7NH2_ICH)
    c3h7cl_basis = thermfit.cbh.basic_basis(C3H7CL_ICH)
    ch3ch2sch3_basis = thermfit.cbh.basic_basis(CH3CH2SCH3_ICH)
    hoch2ch2sch3_basis = thermfit.cbh.basic_basis(HOCH2CH2SCH3_ICH)

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


if __name__ == '__main__':
    test__basis()
