"""
Test the rate plotting functionality for comparing two mechanisms
"""

import tempfile
import numpy as np
#import chemkin_io
import mechanalyzer
import sys


# Set up the thermo data
KTP_DCT = {
    (('H2O2', 'H'), ('H2O', 'HO')):
        {'mech1': {
            'high': np.array([2.6700349e+35, 1.9685318e+36, 3.8313422e+36])},
         'mech2': {
             'high': np.array([5.3289909e+35, 3.9288955e+36, 7.6467869e+36])}},
    (('H2O2', 'HO'), ('H2O', 'HO2')):
        {'mech1': {
            'high': np.array([7.9125421e+35, 2.0715520e+36, 4.9313478e+36])},
         'mech2': {
             'high': np.array([1.2165173e+36, 2.1188423e+36, 3.8957365e+36])}},
    (('H', 'O2'), ('HO2',)):
        {'mech1': {
            'high': np.array([4.3127140e+37, 5.8506465e+37, 6.9933299e+37]),
            1.0: np.array([3.18704869e+37, 1.8013809e+37, 9.1578552e+36]),
            4.0: np.array([3.96279900e+37, 3.7456918e+37, 2.6299543e+37]),
            5.0: np.array([4.02816459e+37, 4.03611586e+37, 3.0049299e+37])},
         'mech2': {
             'high': np.array([2.4577833e+37, 3.3342394e+37, 3.9854461e+37]),
             1.0: np.array([2.33652347e+37, 2.50639834e+37, 2.0178190e+37]),
             4.0: np.array([2.42630350e+37, 3.07992175e+37, 3.2042978e+37]),
             5.0: np.array([2.43253478e+37, 3.12763354e+37, 3.3350312e+37])}},
    (('H2O2',), ('H', 'HO2')):
        {'mech1': {
            'high': np.array([1.5926047e+17, 1.3393927e+28, 6.86564074e+31]),
            1.0: np.array([1.27048152e+17, 5.04819473e+27, 1.27051523e+31]),
            4.0: np.array([1.41955665e+17, 8.10184558e+27, 2.27765507e+31]),
            5.0: np.array([1.43558545e+17, 8.56604416e+27, 2.49261693e+31])},
         'mech2': {
             'high': np.array([5.5741166e+17, 4.6878747e+28, 2.40297426e+32]),
             1.0: np.array([3.15267828e+17, 8.42049048e+27, 1.51635839e+31]),
             4.0: np.array([4.23443851e+17, 1.51569420e+28, 3.68688300e+31]),
             5.0: np.array([4.36337785e+17, 1.65913353e+28, 4.12639542e+31])}},
}
TEMPS = np.array([500.0, 1000.0, 1500.0])

# Set paths to make the plots
#PLOT_PATH = tempfile.mkdtemp()
#print(PLOT_PATH)
PLOT_PATH = sys.argv[1]



def test__plot_rates():
    """ test chemkin_io.plotter.rates
    """
    mechanalyzer.plotter.rates.build(KTP_DCT, TEMPS, dir_prefix=PLOT_PATH)


if __name__ == '__main__':
    test__plot_rates()
