"""
Test the ratefit pressure dependence checker
"""

import numpy as np
import ratefit


TK_DCT = {
    0.1:    [[1200.0, 1500.0], [1.020, 1.050]],
    1.0:    [[1000.0, 1200.0, 1500.0], [1.0, 1.030, 1.060]],
    10.0:   [[1200.0, 1500.0], [1.040, 2.080]],
    100.0:  [[1200.0, 1300.0, 1500.0], [1.050, 1.50, 2.090]],
    'high': [[1200.0, 1500.0], [1.060, 3.040]]
}
TEMP_COMPARE1 = [1200.0]
TEMP_COMPARE2 = [1200.0, 1500.0]
TOL = 20.0
PLOW = None
PHIGH = None

np.set_printoptions(precision=15)


def test__assess_pressure_dependence():
    """ test ratefit.err.assess_pressure_dependence
    """

    is_pdependent1 = ratefit.calc.assess_pressure_dependence(
        TK_DCT, TEMP_COMPARE1, tolerance=TOL, plow=PLOW, phigh=PHIGH)
    is_pdependent2 = ratefit.calc.assess_pressure_dependence(
        TK_DCT, TEMP_COMPARE2, tolerance=TOL, plow=PLOW, phigh=PHIGH)

    assert not is_pdependent1
    assert is_pdependent2


if __name__ == '__main__':
    test__assess_pressure_dependence()
