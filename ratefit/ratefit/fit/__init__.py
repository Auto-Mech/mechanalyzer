"""
Functions to fit rate constants to various functional forms
"""

from ratefit.fit._pdep import pressure_dependent_ktp_dct
from ratefit.fit._pdep import assess_pressure_dependence
from ratefit.fit._err import fitting_error_dct
from ratefit.fit._err import fitting_errors
from ratefit.fit._util import get_valid_tk
from ratefit.fit._util import filter_ktp_dct
from ratefit.fit._util import invert_ktp_dct
from ratefit.fit._util import set_a_conversion_factor

from ratefit.fit._fit import fit_ktp_dct


from ratefit.fit import arrhenius
from ratefit.fit import troe
from ratefit.fit import chebyshev


__all__ = [
    'pressure_dependent_ktp_dct',
    'assess_pressure_dependence',
    'fitting_error_dct',
    'fitting_errors',
    'get_valid_tk',
    'filter_ktp_dct',
    'invert_ktp_dct',
    'set_a_conversion_factor',
    'fit_ktp_dct',
    'arrhenius',
    'troe',
    'chebyshev',
]
