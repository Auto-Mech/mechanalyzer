"""
Module deal with rate constant functional forms; Either
 (1) fits a set of rate constants [k(T, P)] to
     various functional forms
 (2) calculates k(T, P) using a functional form given the
     fitting parameters are known
"""

from ratefit.fxns import single_arrhenius
from ratefit.fxns import double_arrhenius
from ratefit.fxns import lindemann
from ratefit.fxns import troe
from ratefit.err import assess_pressure_dependence
from ratefit.err import calc_sse_and_mae
from ratefit import fit


__all__ = [
    'single_arrhenius',
    'double_arrhenius',
    'lindemann',
    'troe',
    'assess_pressure_dependence',
    'calc_sse_and_mae',
    'fit'
]
