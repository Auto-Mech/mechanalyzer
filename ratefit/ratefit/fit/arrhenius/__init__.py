"""
Functions to fit rate constants to Single or Double Arrhenius Functions
Performs fits either using SciPy or SJK's dsarrfit code
"""

from ratefit.fit.arrhenius._fit import single
from ratefit.fit.arrhenius._fit import double


__all__ = [
    'single',
    'double',
]
