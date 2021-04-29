""" Test read MESS and fit rates
"""

import os
import tempfile
import shutil
import ratefit


# Set path to data files
PATH = os.path.dirname(os.path.realpath(__file__))
MESS_OUT_FILE = os.path.join(PATH, 'data', 'c4h9o3_rate.out')
TMP_DIR = tempfile.mkdtemp()
print('Temp Run Dir:', TMP_DIR)

shutil.copy(MESS_OUT_FILE, os.path.join(TMP_DIR, 'rate.out'))

# Data
PES_FORMULA = 'CH4'
FIT_METHOD1 = 'arrhenius'
FIT_METHOD2 = 'chebyshev'

# influences the direction of the reactions obtained from MESS if forw and rev not taken
# no explicit control inside of mechdriver
LABEL_DCT = {
    'WELL1': 'w1',
    'WELL2': 'w2',
    'WELL3': 'w3',
    'WELL4': 'w4',
    'BIMOL1': 'p1',
    'BIMOL2': 'p2',
    'BIMOL3': 'p3',
    'BIMOL4': 'p4',
    'BIMOL5': 'p5',
    'BIMOL6': 'p6',
    'BIMOL7': 'p7'
}
FIT_TEMPS = (300, 400, 500, 600, 700, 800, 900, 1000, 1100,
             1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000)


def test__full_arrfit():
    """ test ratefit.fit.ktp_dct
    """

    arrfit_dct = ratefit.fit.fit_ktp_dct(
        TMP_DIR, FIT_METHOD1, label_dct=LABEL_DCT)
    print('\n\n\n-------\n\n\n')
    print('CHEMKIN STR\n\n')
    arrfit_str = ''
    for key, val in arrfit_dct.items():
        arrfit_str += val
        arrfit_str += '\n\n'
    print(arrfit_str)


def test__full_chebfit():
    """ test ratefit.fit.ktp_dct
    """

    chebfit_dct = ratefit.fit.fit_ktp_dct(
        TMP_DIR, FIT_METHOD2, label_dct=LABEL_DCT)
    print('\n\n\n-------\n\n\n')
    print('CHEMKIN STR\n\n')
    chebfit_str = ''
    for key, val in chebfit_dct.items():
        chebfit_str += val
        chebfit_str += '\n\n'
    print(chebfit_str)


if __name__ == '__main__':
    test__full_arrfit()
    # test__full_chebfit()
