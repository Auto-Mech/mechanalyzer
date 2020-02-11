"""
Module deal with rate constant functional forms; Either
 (1) fits a set of rate constants [k(T, P)] to
     various functional forms
 (2) calculates k(T, P) using a functional form given the
     fitting parameters are known
"""

from ratefit import calc
from ratefit import fit


__all__ = [
    'calc',
    'fit'
]
