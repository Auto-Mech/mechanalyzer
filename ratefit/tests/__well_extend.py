""" test
"""

import ratefit


with open('data/c3h3.inp', 'r') as fobj:
    INP_STR = fobj.read()
with open('data/c3h3.out', 'r') as fobj:
    OUT_STR = fobj.read()
with open('data/c3h3.aux', 'r') as fobj:
    AUX_STR = fobj.read()
with open('data/c3h3.logf', 'r') as fobj:
    LOG_STR = fobj.read()

LUMP_PRESSURE = 1.0  # fix
LUMP_TEMP = 1200


def test__():
    """ test ratefit.fit._wellextend.well_lumped_input_file
    """

    inp_str = ratefit.fit.well_lumped_input_file(
        INP_STR, OUT_STR, AUX_STR, LOG_STR,
        LUMP_PRESSURE, LUMP_TEMP)

    print('INIT\n')
    print('------------------------------------------')
    print(INP_STR)
    print('\n------------------------------------------')
    print('\n\nLUMP\n')
    print(inp_str)


test__()
