import time
import numpy
from ratefit.fit import plog
from ratefit.fit import err

TEMPS = numpy.linspace(400, 900, 11)
TEMPS_SHORT = numpy.linspace(450, 900, 10)  # used with KTS1 since no 400 K val

# Define rates at 0.01, 0.1, 1, 10, and 100 atm 
KTS1 = numpy.array([5.45E-05, 0.000126326, 0.000302054, 0.00107086, 
                    0.00553726, 0.0416265, 0.369654, 3.08305, 
                    20.4691, 130.603,])  # leaving off 400 K value
KTS2 = numpy.array([0.00373048, 0.0124119, 0.0307059, 0.066723, 
                    0.162508, 0.489855, 2.14679, 13.2429, 
                    87.8471, 462.577, 2226.93,])
KTS3 = numpy.array([0.508358, 2.4336, 7.23671, 16.2493, 
                    33.3004, 68.2562, 161.142, 479.797, 
                    1768.32, 6495.67, 23780,])
KTS4 = numpy.array([12.9012, 101.173, 430.225, 1199.49, 
                    2600.21, 4880.51, 8806.79, 16446.2, 
                    33401.4, 73633.2, 186461,])
KTS5 = numpy.array([51.9327, 655.99, 4402.47, 18140.5, 
                    51929.1, 114477, 213629, 362666, 
                    589511, 946634, 1.84E+06,])
KTP_DCT = {
    0.01:  (TEMPS_SHORT, KTS1),
    0.1:   (TEMPS, KTS2),
    1.0:   (TEMPS, KTS3),
    10.0:  (TEMPS, KTS4),
    100.0: (TEMPS, KTS5),}


def test_plog():

    params, err_dct = plog.get_params(KTP_DCT, dbltol=15, dbl_iter=1)
    max_err = err.get_max_err(err_dct)
    assert max_err < 25  # error in %


if __name__ == '__main__':
    test_plog()
