""" test chemkin_io.calculator.mechanism
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
TROE_REACTION = GROUPED_REACTION_STRS[3]
PLOG_REACTION = GROUPED_REACTION_STRS[5]
DUP_PLOG_REACTION = GROUPED_REACTION_STRS[6]
CHEBYSHEV_REACTION = GROUPED_REACTION_STRS[7]

# Set temperatures and pressures
T_REF = 1.0
UNITS = ('cal/mole', 'moles')
TEMPS = numpy.array([500.0, 1000.0, 1500.0, 2000.0])
PRESSURES = [1.0, 5.0, 10.0]
PRESSURES2 = [0.0100, 0.0700, 0.987]
PRESSURES3 = [0.1, 0.5, 2]
PRESSURES4 = ['high']


def test__high_p_rate_constants():
    """ test chemkin_io.calculator.rates.reaction
        for a reaction with only high-pressure params
    """
    ktp_dct1 = chemkin_io.calculator.rates.reaction(
        HIGHP_REACTION, UNITS, T_REF, TEMPS, pressures=PRESSURES4)
    ktp_dct2 = chemkin_io.calculator.rates.reaction(
        DUP_HIGHP_REACTION, UNITS, T_REF, TEMPS, pressures=PRESSURES4)

    assert ktp_dct1.keys() == ktp_dct2.keys() == ['high']
    rates1 = ktp_dct1['high']
    rates2 = ktp_dct2['high']

    ref_rates1 = numpy.array(
        [4.43369728e+11, 3.26882402e+12, 6.36209341e+12, 8.87573427e+12])
    ref_rates2 = numpy.array(
        [1.31390851e+12, 3.43989296e+12, 8.18869558e+12, 1.37943670e+13])

    assert numpy.allclose(rates1, ref_rates1, atol=0.0001)
    assert numpy.allclose(rates2, ref_rates2, atol=0.0002)


def test__lindemann_rate_constants():
    """ test chemkin_io.calculator.rates.reaction
        for a reaction with high-pressure and low-pressure params
    """

    ktp_dct = chemkin_io.calculator.rates.reaction(
        LINDEMANN_REACTION, UNITS, T_REF, TEMPS, pressures=PRESSURES)

    assert ktp_dct.keys() == PRESSURES
    rates1 = ktp_dct[1.0]
    rates2 = ktp_dct[5.0]
    rates3 = ktp_dct[10.0]

    # [5.29221873e+13, 2.99126335e+13, 1.52069760e+13, 8.61077092e+12]
i   # [6.68892451e+13, 6.70212791e+13, 4.98980340e+13, 3.41335883e+13]
    # [6.91711751e+13, 7.93217841e+13, 6.98028298e+13, 5.42239489e+13]


def test__troe_rate_constants():
    """ test chemkin_io.calculator.rates.reaction
        for a reaction with only high-pressure, low-pressure, and Troe params
    """
    ktp_dct = chemkin_io.calculator.rates.reaction(
        TROE_REACTION, UNITS, T_REF, TEMPS, pressures=PRESSURES)
    # [2.10968415e-07 8.38272437e+03 2.10974013e+07 8.26369002e+08]
    # [2.38384566e-07 1.42242506e+04 4.13908769e+07 1.94896432e+09]
    # [2.45064585e-07 1.63156371e+04 5.39828882e+07 2.59555337e+09]


def test__plog_rate_constants():
    """ test chemkin_io.calculator.rates.reaction
        for a reaction with only high-pressure and PLog params
    """
    ktp_dct = chemkin_io.calculator.rates.reaction(
        PLOG_REACTION, UNITS, T_REF, TEMPS, pressures=PRESSURES2)
    ktp_dct = chemkin_io.calculator.rates.reaction(
        DUP_PLOG_REACTION, UNITS, T_REF, TEMPS, pressures=PRESSURES3)
    # [2.19605076e-03 3.82242886e+01 8.58251750e+02 3.79239652e+03]
    # [1.15498396e-01 2.90865625e+03 5.31130047e+04 1.80240066e+05]
    # [6.12150518e+00 1.86775646e+05 2.38952477e+06 5.53241147e+06]
    # [9.14719356e+11 1.42260637e+12 1.70609757e+12 1.90158713e+12]
    # [8.37486568e+11 1.39876603e+12 1.68929974e+12 1.88214742e+12]
    # [6.64709527e+11 1.36151762e+12 1.69077959e+12 1.89116549e+12]


def test__chebyshev_rate_constants():
    """ test chemkin_io.calculator.rates.reaction
        for a reaction with only high-pressure and Chebyshev params
    """
    ktp_dct = chemkin_io.calculator.rates.reaction(
        CHEBYSHEV_REACTION, UNITS, T_REF, TEMPS, pressures=PRESSURES)
    print('\nChebyshev rate_constants')
    for key, val in ktp_dct.items():
        print(key)
        print(val)
    # [1.28848565e+06 5.35177959e+10 7.50345079e+09 1.88855517e+07]
    # [1.75978967e+08 3.00571254e+14 2.12088277e+13 6.59216199e+09]
    # [1.46257915e+09 1.23807607e+16 6.49967375e+14 8.20719048e+10]


def test__mechanism():
    """ test chemkin_io.calculator.reaction.mechanism
    """
    ktp_dct = chemkin_io.calculator.rates.mechanism(
        FAKE1_REACTION_BLOCK, UNITS, T_REF, TEMPS, pressures=PRESSURES)
    for spc, ktp in ktp_dct.items():
        print(spc)
        print(ktp)


def test__branching_fractions():
    """ test chemkin_io.calculator.reaction.branching_fractions
    """
    print('\n\n')
    mech_dct = chemkin_io.calculator.rates.mechanism(
        FAKE1_REACTION_BLOCK, UNITS, T_REF, TEMPS, pressures=PRESSURES)
    tot_dct, branch_dct = chemkin_io.calculator.rates.branching_fractions(
        mech_dct, PRESSURES)
    for spc, ktp in tot_dct.items():
        print(spc)
        print(ktp)
    print('\n\n')
    for spc, ktp in branch_dct.items():
        print(spc)
        print(ktp)


if __name__ == '__main__':
    test__high_p_rate_constants()
    test__lindemann_rate_constants()
    test__troe_rate_constants()
    test__plog_rate_constants()
    test__chebyshev_rate_constants()
    test__mechanism()
    test__branching_fractions()
