""" Functions for dealing with a ktp dictionary
"""

import numpy


def read(dct, rxn, pressure, val):
    """ Read the entries of a ktp dictionary. Right now assumes val exists
    """

    assert pressure == 'high' or isinstance(pressure, float)
    assert val in ('temps', 'rates')

    # Deal with pressure floats
    if pressure == 'high':
        pval = 'high'
    else:
        for _pressure in (p for p in dct[rxn].keys() if p != 'high'):
            if numpy.isclose(pressure, _pressure, atol=1.0e-4):
                pval = _pressure

    # Set final val index for reading rates versus temp array
    val_idx = 1 if val == 'rates' else 0

    return dct[rxn][pval][val_idx]
    # a lot of checks below, will be slow...
    # rxn_dct = dct.get(None)
    # if rxn_dct is not None:
    # p_dct = rxn_dct.get(pressure)
    # if p_dct is not None:
    # value = p_dct[validx]
    # return value


def check_p_t(pressures, temps):
    """ Enforces certain rules regarding the pressure and temperature arrays.

        :param pressures: array of pressures
        :type pressures: numpy.ndarray
        :param temps: array of temps
        :type temps: numpy.ndarray
    """

    # Check that the dimensionality of the temps array is either 1 or 2
    temp_dim = numpy.ndim(temps)
    assert temp_dim in (1, 2), (
        f'The dimensionality of temps is {temp_dim};',
        'it should be either 1 or 2'
    )

    # If temps is 2-D, enforce that the number of values in
    # each temp array matches the number of pressures
    if temp_dim == 2:
        len_temps = numpy.shape(temps)[1]
        len_pressures = len(pressures)
        assert len_pressures == len_temps, (
            f'# of pressures is {len_pressures},',
            f'while # of temps in each array is {len_temps}'
        )
