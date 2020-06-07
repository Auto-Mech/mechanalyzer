""" test chemkin_io.parser.mechanism
"""

from __future__ import unicode_literals
from builtins import open
import os
import numpy
import chemkin_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


# Set paths
PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
FAKE1_MECH_NAME = 'fake1_mech.txt'

# Read mechanism files
FAKE1_MECH_STR = _read_file(
    os.path.join(DATA_PATH, FAKE1_MECH_NAME))

# Build the reactions blocks and data strings
FAKE1_REACTION_BLOCK = chemkin_io.parser.mechanism.reaction_block(
    FAKE1_MECH_STR)
FAKE1_REACTION_STRS = chemkin_io.parser.reaction.data_strings(
    FAKE1_REACTION_BLOCK)
FAKE1_REACTION_DCT = chemkin_io.parser.reaction.data_dct(
    FAKE1_REACTION_BLOCK)

# Set the reactions to elements of the string, use dct to have DUPs together
GROUPED_REACTION_STRS = list(FAKE1_REACTION_DCT.values())
HIGHP_REACTION = GROUPED_REACTION_STRS[0]
DUP_HIGHP_REACTION = GROUPED_REACTION_STRS[1]
LINDEMANN_REACTION = GROUPED_REACTION_STRS[2]
TROE1_REACTION = GROUPED_REACTION_STRS[3]
TROE2_REACTION = GROUPED_REACTION_STRS[4]
PLOG_REACTION = GROUPED_REACTION_STRS[5]
DUP_PLOG_REACTION = GROUPED_REACTION_STRS[6]
CHEBYSHEV_REACTION = GROUPED_REACTION_STRS[7]
FAKE_HIGHP_REACTION2 = GROUPED_REACTION_STRS[9]


def test__data_objs():
    """ test chemkin_io.parser.reaction.data_strings
        test chemkin_io.parser.reaction.data_dct
    """

    rxn_strs = chemkin_io.parser.reaction.data_strings(
        FAKE1_REACTION_BLOCK)
    rxn_dct = chemkin_io.parser.reaction.data_dct(
        FAKE1_REACTION_BLOCK)

    assert len(rxn_strs) == 13
    assert len(rxn_dct) == 10


def test__reactant_names():
    """ test chemkin_io.parser.reaction.reactant_names
    """
    names = chemkin_io.parser.reaction.reactant_names(PLOG_REACTION)
    test_names = ('HOCO',)
    assert names == test_names


def test__product_names():
    """ test chemkin_io.parser.reaction.product_names
    """
    names = chemkin_io.parser.reaction.product_names(PLOG_REACTION)
    test_names = ('CO', 'OH')
    assert names == test_names


def test__high_p_parameters():
    """ test chemkin_io.parser.reaction.high_p_parameters
        test chemkin_io.parser.reaction.are_highp_fake
    """
    params1 = chemkin_io.parser.reaction.high_p_parameters(
        HIGHP_REACTION)
    params2 = chemkin_io.parser.reaction.high_p_parameters(
        DUP_HIGHP_REACTION)
    params3 = chemkin_io.parser.reaction.high_p_parameters(
        LINDEMANN_REACTION)
    params4 = chemkin_io.parser.reaction.high_p_parameters(
        TROE1_REACTION)
    params5 = chemkin_io.parser.reaction.high_p_parameters(
        PLOG_REACTION)
    params6 = chemkin_io.parser.reaction.high_p_parameters(
        DUP_PLOG_REACTION)
    params7 = chemkin_io.parser.reaction.high_p_parameters(
        CHEBYSHEV_REACTION)
    params8 = chemkin_io.parser.reaction.high_p_parameters(
        FAKE_HIGHP_REACTION2)

    ref_params1 = [[24100000000000.0, 0.0, 3970.0]]
    ref_params2 = [[1740000000000.0, 0.0, 318.0],
                   [75900000000000.0, 0.0, 7269.0]]
    ref_params3 = [[4650000000000.0, 0.44, 0.0]]
    ref_params4 = [[2000000000000.0, 0.9, 48749.0]]
    ref_params5 = [[6.3e+32, -5.96, 32470.0]]
    ref_params6 = [[1770000000000.0, 0.16, 4206],
                   [1760000000000.0, 0.25, 4305]]
    ref_params7 = [[1.0, 0.0, 0.0]]
    ref_params8 = [[1.0, 0.0, 0.0],
                   [1.0, 0.0, 0.0]]

    assert numpy.allclose(params1, ref_params1)
    assert numpy.allclose(params2, ref_params2)
    assert numpy.allclose(params3, ref_params3)
    assert numpy.allclose(params4, ref_params4)
    assert numpy.allclose(params5, ref_params5)
    assert numpy.allclose(params6, ref_params6)
    assert numpy.allclose(params7, ref_params7)
    assert numpy.allclose(params8, ref_params8)

    arefake1 = chemkin_io.parser.reaction.are_highp_fake(
        params1)
    arefake2 = chemkin_io.parser.reaction.are_highp_fake(
        params2)
    arefake3 = chemkin_io.parser.reaction.are_highp_fake(
        params3)
    arefake4 = chemkin_io.parser.reaction.are_highp_fake(
        params4)
    arefake5 = chemkin_io.parser.reaction.are_highp_fake(
        params5)
    arefake6 = chemkin_io.parser.reaction.are_highp_fake(
        params6)
    arefake7 = chemkin_io.parser.reaction.are_highp_fake(
        params7)
    arefake8 = chemkin_io.parser.reaction.are_highp_fake(
        params8)

    assert not arefake1
    assert not arefake2
    assert not arefake3
    assert not arefake4
    assert not arefake5
    assert not arefake6
    assert arefake7
    assert arefake8


def test__low_p_parameters():
    """ test chemkin_io.parser.reaction.low_p_parameters
    """
    params1 = chemkin_io.parser.reaction.low_p_parameters(
        LINDEMANN_REACTION)
    params2 = chemkin_io.parser.reaction.low_p_parameters(
        TROE1_REACTION)
    params3 = chemkin_io.parser.reaction.low_p_parameters(
        TROE2_REACTION)

    ref_params1 = [[1.737e+19, -1.23, 0.0]]
    ref_params2 = [[2.49e+24, -2.3, 48749.0]]
    ref_params3 = [[2.49e+24, -2.4, 58.75]]

    assert numpy.allclose(params1, ref_params1)
    assert numpy.allclose(params2, ref_params2)
    assert numpy.allclose(params3, ref_params3)


def test__troe_parameters():
    """ test chemkin_io.parser.reaction.troe_parameters
    """
    params1 = chemkin_io.parser.reaction.troe_parameters(
        TROE1_REACTION)
    params2 = chemkin_io.parser.reaction.troe_parameters(
        TROE2_REACTION)

    ref_params1 = [0.43, 1e-30, 1e+30, None]
    ref_params2 = [0.58, 30.0, 90000.0, 90000.0]

    assert numpy.allclose(params1[0:3], ref_params1[0:3])
    assert params1[3] is None and ref_params1[3] is None
    assert numpy.allclose(params2, ref_params2)


def test__buffer_enhance_factors():
    """ test chemkin_io.parser.reaction.buffer_enhance_factors
    """
    fct_dct1 = chemkin_io.parser.reaction.buffer_enhance_factors(
        TROE1_REACTION)
    fct_dct2 = chemkin_io.parser.reaction.buffer_enhance_factors(
        TROE2_REACTION)

    ref_fct_dct1 = {
        'H2O': 7.65,
        'CO2': 1.6,
        'N2': 1.5,
        'O2': 1.2,
        'HE': 0.65,
        'H2O2': 7.7,
        'H2': 3.7,
        'CO': 2.8
    }

    ref_fct_dct2 = {
        'H2O': 6.63,
        'H2': 3.27,
        'N2': 1.33,
        'CO': 2.8,
        'H2O2': 6.61,
        'CO2': 1.6,
        'O2': 1.2
    }

    fct_inf1 = zip(fct_dct1.items(), ref_fct_dct1.items())
    for (key, vals), (ref_key, ref_vals) in fct_inf1:
        assert key == ref_key
        assert numpy.isclose(vals, ref_vals)
    fct_inf2 = zip(fct_dct2.items(), ref_fct_dct2.items())
    for (key, vals), (ref_key, ref_vals) in fct_inf2:
        assert key == ref_key
        assert numpy.isclose(vals, ref_vals)


def test__plog_parameters():
    """ test chemkin_io.parser.reaction.plog_parameters
    """
    params1 = chemkin_io.parser.reaction.plog_parameters(
        PLOG_REACTION)
    params2 = chemkin_io.parser.reaction.plog_parameters(
        DUP_PLOG_REACTION)

    ref_params1 = {
        0.001: [[1.55e-08, 2.93, 8768.0]],
        0.003: [[1770.0, 0.34, 18076.0]],
        0.0296: [[20200000000000.0, -1.87, 22755.0]],
        0.0987: [[1.68e+18, -3.05, 24323.0]],
        0.2961: [[2.5e+24, -4.63, 27067.0]],
        0.9869: [[4.54e+26, -5.12, 27572.0]],
        2.9607: [[7.12e+28, -5.6, 28535.0]],
        9.869: [[5.48e+29, -5.7, 28899.0]],
        29.607: [[9.89e+31, -6.19, 30518.0]],
        98.69: [[5.74e+33, -6.53, 32068.0]],
        296.07: [[2.61e+33, -6.29, 32231.0]],
        986.9: [[6.3e+32, -5.96, 32470.0]]
    }

    ref_params2 = {
        0.01: [[7.88E+20, -2.67, 6742.0],
               [13600000000.0, 0.62, -277.6]],
        0.1: [[7.72E+20, -2.67, 6713.0],
              [14200000000.0, 0.62, -247.7]],
        0.316: [[9.87E+20, -2.7, 6724.0],
                [16600000000.0, 0.6, -162.5]],
        1.0: [[7.10E+20, -2.65, 6489.0],
              [20200000000.0, 0.58, 38.4]],
        3.16: [[4.50E+20, -2.53, 6406.0],
               [9750000000.0, 0.67, 248.0]],
        10.0: [[1.76E+23, -3.22, 8697.0],
               [7340000000.0, 0.72, 778.1]],
        31.6: [[3.14E+25, -3.77, 11530.0],
               [1570000000.0, 0.92, 1219.0]],
        100.0: [[1.02E+26, -3.8, 13910.0],
                [78500000.0, 1.28, 1401.0]]
    }

    par_inf1 = zip(params1.items(), ref_params1.items())
    for (key, vals), (ref_key, ref_vals) in par_inf1:
        assert numpy.isclose(key, ref_key)
        assert numpy.allclose(vals, ref_vals)
    par_inf2 = zip(params2.items(), ref_params2.items())
    for (key, vals), (ref_key, ref_vals) in par_inf2:
        assert numpy.isclose(key, ref_key)
        assert numpy.allclose(vals, ref_vals)


def test__chebyshev_parameters():
    """ test chemkin_io.parser.reaction.chebyshev_parameters
    """
    params = chemkin_io.parser.reaction.chebyshev_parameters(
        CHEBYSHEV_REACTION)

    ref_params = {
        't_limits': [300.0, 2200.0],
        'p_limits': [0.01, 98.702],
        'alpha_dim': [6, 4],
        'alpha_elm': [[8.684, 8.684, 1.879e-15, 1.879e-15],
                      [-0.2159, -0.2159, 2.929e-17, 2.929e-17],
                      [-1.557e-15, -1.557e-15, -8.346e-31, -8.346e-31],
                      [0.2159, 0.2159, -2.929e-17, -2.929e-17],
                      [-2.684, -2.684, -1.879e-15, -1.879e-15],
                      [0.2159, 0.2159, -2.929e-17, -2.929e-17]]
    }

    assert params['t_limits'] == ref_params['t_limits']
    assert params['p_limits'] == ref_params['p_limits']
    assert params['alpha_dim'] == ref_params['alpha_dim']
    for (row, ref_row) in zip(params['alpha_elm'], ref_params['alpha_elm']):
        assert numpy.allclose(row, ref_row)


def test__ratek_fit_info():
    """ test chemkin_io.parser.reaction.ratek_fit_info
    """

    params = chemkin_io.parser.reaction.ratek_fit_info(DUP_PLOG_REACTION)

    ref_params = {
        1.0: {'temps': [100, 500], 'mean_err': 1.1, 'max_err': 11.1},
        10.0: {'temps': [100, 700], 'mean_err': 2.2, 'max_err': 22.2},
        100.0: {'temps': [100, 1000], 'mean_err': 3.3, 'max_err': 33.3},
        'High': {'temps': [100, 1500], 'mean_err': 4.4, 'max_err': 44.4}
    }

    assert set(params.keys()) == set(ref_params.keys())
    for key in params:
        assert numpy.allclose(params[key]['temps'], params[key]['temps'])
        assert numpy.isclose(params[key]['mean_err'], params[key]['mean_err'])
        assert numpy.isclose(params[key]['max_err'], params[key]['max_err'])


if __name__ == '__main__':
    test__data_objs()
    test__reactant_names()
    test__product_names()
    test__high_p_parameters()
    test__low_p_parameters()
    test__troe_parameters()
    test__buffer_enhance_factors()
    test__plog_parameters()
    test__chebyshev_parameters()
    test__ratek_fit_info()
