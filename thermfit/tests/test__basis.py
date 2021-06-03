""" test thermfit.basis
"""

import numpy
import thermfit.basis


# Formula dicts for various species
C3H8_FML = {'H': 8, 'C': 3}
C3H7OH_FML = {'H': 8, 'C': 3, 'O': 1}
C3H7NH2_FML = {'H': 9, 'C': 3, 'N': 1}
C3H7CL_FML = {'H': 7, 'C': 3, 'Cl': 1}
CH3CH2SCH3_FML = {'H': 8, 'C': 3, 'S': 1}
HOCH2CH2SCH3_FML = {'H': 8, 'C': 3, 'O': 1, 'S': 1}


def test__():
    """ test basis
    """

    c3h8_basis = thermfit.basis.basis_species(C3H8_FML)
    c3h7oh_basis = thermfit.basis.basis_species(C3H7OH_FML)
    c3h7nh2_basis = thermfit.basis.basis_species(C3H7NH2_FML)
    c3h7cl_basis = thermfit.basis.basis_species(C3H7CL_FML)
    ch3ch2sch3_basis = thermfit.basis.basis_species(CH3CH2SCH3_FML)
    hoch2ch2sch3_basis = thermfit.basis.basis_species(HOCH2CH2SCH3_FML)

    assert c3h8_basis == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4')
    assert c3h7oh_basis == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/H2O/h1H2')
    assert c3h7nh2_basis == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/H3N/h1H3')
    assert c3h7cl_basis == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/ClH/h1H')
    assert ch3ch2sch3_basis == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/O2S/c1-3-2',
        'InChI=1S/H2O/h1H2')
    assert hoch2ch2sch3_basis == (
        'InChI=1S/H2/h1H',
        'InChI=1S/CH4/h1H4',
        'InChI=1S/H2O/h1H2',
        'InChI=1S/O2S/c1-3-2')

    c3h8_cffs = thermfit.basis.basis_coefficients(
        c3h8_basis, C3H8_FML)
    c3h7oh_cffs = thermfit.basis.basis_coefficients(
        c3h7oh_basis, C3H7OH_FML)
    c3h7nh2_cffs = thermfit.basis.basis_coefficients(
        c3h7nh2_basis, C3H7NH2_FML)
    c3h7cl_cffs = thermfit.basis.basis_coefficients(
        c3h7cl_basis, C3H7CL_FML)
    ch3ch2sch3_cffs = thermfit.basis.basis_coefficients(
        ch3ch2sch3_basis, CH3CH2SCH3_FML)
    hoch2ch2sch3_cffs = thermfit.basis.basis_coefficients(
        hoch2ch2sch3_basis, HOCH2CH2SCH3_FML)

    assert numpy.allclose(c3h8_cffs, numpy.array([-2.,  3.]))
    assert numpy.allclose(c3h7oh_cffs, numpy.array([-3.,  3.,  1.]))
    assert numpy.allclose(c3h7nh2_cffs, numpy.array([-3.,  3.,  1.]))
    assert numpy.allclose(c3h7cl_cffs, numpy.array([-3.,  3.,  1.]))
    assert numpy.allclose(ch3ch2sch3_cffs, numpy.array([0.,  3.,  1., -2.]))
    assert numpy.allclose(hoch2ch2sch3_cffs, numpy.array([-1.,  3., -1.,  1.]))


if __name__ == '__main__':
    test__()
