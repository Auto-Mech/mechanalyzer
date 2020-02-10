"""
Functions to fit rate constants to Single or Double Arrhenius Functions
Performs fits either using SciPy or SJK's dsarrfit code
"""

from ratefit.fit.troe.fit import single


__all__ = [
    'single',
]
