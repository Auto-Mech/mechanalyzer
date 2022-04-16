"""
calculates derived quantities from the strings of the species
"""

from mechanalyzer.calculator import rates
from mechanalyzer.calculator import thermo
from mechanalyzer.calculator import combine
from mechanalyzer.calculator import compare
from mechanalyzer.calculator import statmodels
from mechanalyzer.calculator import bf
from mechanalyzer.calculator._prompt import prompt_dissociation_ktp_dct
from mechanalyzer.calculator._prompt import multipes_prompt_dissociation_ktp_dct

__all__ = [
    'rates',
    'thermo',
    'combine',
    'compare',
    'statmodels',
    'bf',
    'prompt_dissociation_ktp_dct',
    'multipes_prompt_dissociation_ktp_dct'
]
