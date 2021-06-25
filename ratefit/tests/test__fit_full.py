""" Test read MESS and fit rates

    NEED TO MAKE A REAL TEST
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
FIT_METHOD1 = 'arrhenius'

# influences direction of reactions from MESS if forw and rev not taken
# no explicit control inside of mechdriver
LABEL_DCT = {
    'w1': 'WELL1',
    'w2': 'WELL2',
    'w3': 'WELL3',
    'w4': 'WELL4',
    'p1': 'BIMOL1',
    'p2': 'BIMOL2',
    'p3': 'BIMOL3',
    'p4': 'BIMOL4',
    'p5': 'BIMOL5',
    'p6': 'BIMOL6',
    'p7': 'BIMOL7'
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
    for _, val in arrfit_dct.items():
        arrfit_str += val
        arrfit_str += '\n\n'
    print(arrfit_str)


if __name__ == '__main__':
    test__full_arrfit()
