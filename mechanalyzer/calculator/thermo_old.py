""" Calculates thermodynamic quantities using thermo data strings
"""


import numpy as np
from ioformat import phycon
from chemkin_io.parser import thermo as thm_parser


# functions which calculate quantiies using data from the thermo section #
def mechanism(block_str, temps, rval=phycon.RC):
    """ Parses the all the reactions data string in the thermo block
        in a mechanism file for their NASA polynomials and
        uses them to calculate thermochemical values: H(T), Cp(T), S(T), G(T).

        :param block_str: string of Reaction block of ChemKin input
        :type block_str: str
        :param temps: temperatures to calculate Thermo quantities (K)
        :type temps: list(float)
        :return mech_thermo_dct: dct of thermo data [H(T), Cp(T), S(T), G(T)]
        :rtype: dict[spc: [[H(T)], [Cp(T)], [S(T)], [G(T)]]]
    """

    nasa_dct = thm_parser.data_dct(block_str)

    mech_thermo_dct = {}
    for name, thermo_dstr in nasa_dct.items():
        # print('dstr\n', thermo_dstr)
        h_t, cp_t, s_t, g_t, = [], [], [], []
        for temp in temps:
            h_t.append(enthalpy(thermo_dstr, temp, rval=rval))
            cp_t.append(heat_capacity(thermo_dstr, temp, rval=rval))
            s_t.append(entropy(thermo_dstr, temp, rval=rval))
            g_t.append(gibbs(thermo_dstr, temp, rval=rval))

        mech_thermo_dct[name] = [h_t, cp_t, s_t, g_t]

    return mech_thermo_dct


def enthalpy(thm_dstr, temp, rval=phycon.RC):
    """ Calculate the Enthalpy [H(T)] of a species using the
        coefficients of its NASA polynomial.

        :param thm_dstr: string containing NASA polynomial of species
        :type thm_dstr: str
        :param temp: temperature to calculate Enthalpy (K)
        :type temp: float
        :return h_t: value for the Enthalpy (???)
        :rtype: float
    """

    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        h_t = (
            cfts[0] +
            ((cfts[1] * temp) / 2.0) +
            ((cfts[2] * temp**2) / 3.0) +
            ((cfts[3] * temp**3) / 4.0) +
            ((cfts[4] * temp**4) / 5.0) +
            (cfts[5] / temp)
        )
        h_t *= (rval * temp)
    else:
        h_t = None

    return h_t


def heat_capacity(thm_dstr, temp, rval=phycon.RC):
    """ Calculate the Heat Capacity [Cp(T)] of a species using the
        coefficients of its NASA polynomial.

        :param thm_dstr: string containing NASA polynomial of species
        :type thm_dstr: str
        :param temp: temperature to calculate heat capacity (K)
        :type temp: float
        :return cp_t: value for the Heat Capacity
        :rtype: float
    """
    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        cp_t = (
            cfts[0] +
            (cfts[1] * temp) +
            (cfts[2] * temp**2) +
            (cfts[3] * temp**3) +
            (cfts[4] * temp**4)
        )
        cp_t *= rval
    else:
        cp_t = None

    return cp_t


def entropy(thm_dstr, temp, rval=phycon.RC):
    """ Calculate the Entropy [S(T)] of a species using the
        coefficients of its NASA polynomial.

        :param thm_dstr: string containing NASA polynomial of species
        :type thm_dstr: str
        :param float temp: temperature to calculate Entropy
        :return s_t: value for the Entropy
        :rtype: float
    """
    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        s_t = (
            (cfts[0] * np.log(temp)) +
            (cfts[1] * temp) +
            ((cfts[2] * temp**2) / 2.0) +
            ((cfts[3] * temp**3) / 3.0) +
            ((cfts[4] * temp**4) / 4.0) +
            (cfts[6])
        )
        s_t *= rval
    else:
        s_t = None

    return s_t


def gibbs(thm_dstr, temp, rval=phycon.RC):
    """ Calculate the Gibbs Free Energy [G(T)] of a species using the
        coefficients of its NASA polynomial.

        :param thm_dstr: string containing NASA polynomial of species
        :type thm_dstr: str
        :param temp: temperature to calculate Gibbs Free Energy
        :type temp: float
        :return g_t: value for the Gibbs Free Energy
        :rtype: float
    """

    h_t = enthalpy(thm_dstr, temp, rval=rval)
    s_t = entropy(thm_dstr, temp, rval=rval)
    if h_t is not None and s_t is not None:
        g_t = h_t - (s_t * temp)
    else:
        g_t = None

    return g_t


def _coefficients_for_specific_temperature(thm_dstr, temp):
    """ Parse out the coefficients of a NASA polynomial from
        a ChemKin-formatted string. The input temperature value
        determines whether the low- or high-temperature coefficients
        are read from the string.

        :param thm_dstr: String containing NASA polynomial of species
        :type thm_dstr: str
        :param temp: Temperature used to read the coefficients
        :type temp: float
        :return cfts: low- or high-temperature coefficients of NASA polynomial
        :rtype: list(float)
    """

    temps = thm_parser.temperatures(thm_dstr)
    if temps[0] <= temp <= temps[1]:
        cfts = thm_parser.low_coefficients(thm_dstr)
    elif temps[1] < temp < temps[2]:
        cfts = thm_parser.high_coefficients(thm_dstr)
    else:
        cfts = None

    return cfts
