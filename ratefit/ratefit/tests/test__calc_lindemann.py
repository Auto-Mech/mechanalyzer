"""
Test the ratefit rate constant calculators for Lindemann and Troe expressions
"""

import numpy
import ratefit


# K(T) values at high- and low-pressure limits (T = [300, 3000]K, 300 K Steps)
HIGHP_KTS = numpy.array(
    [3.364246691071372e+14, 6.303626542377624e+14, 9.092106004369605e+14,
     1.178705457518608e+15, 1.441457594471313e+15, 1.698960838677903e+15,
     1.952180784432845e+15, 2.201791555928543e+15, 2.448289643769715e+15,
     2.692055417308714e+15])
LOWP_KTS = numpy.array(
    [4.957573597808111e+18, 1.010824653785529e+18, 3.983435704402545e+17,
     2.056812288924121e+17, 1.231626977413509e+17, 8.099919515048618e+16,
     5.683130672639794e+16, 4.180898071689080e+16, 3.189099736360676e+16,
     2.503025834637141e+16])

# A_HIGH, N_HIGH, EA_HIGH = 2.000e+12, 0.900, 4.87490
# A_LOW, N_LOW, EA_LOW = 2.490e24, -2.300, 4.87490

# Pressure and Temperature Range and Parameter Setting
PRESSURES = numpy.array([0.1, 2.0, 5.0, 10.0])
TEMPS = numpy.arange(300.0, 3300.0, 300.0)
T_REF = 1.0

# Set Reference Data for Comparison
REF_LIND_KTPS = {
    0.1: numpy.array([1.90012212e+13, 2.04642017e+12, 5.39063794e+11, 2.08842501e+11,
                      1.00055327e+11, 5.48373687e+10, 3.29794522e+10, 2.12293442e+10,
                      1.43940899e+10, 1.01677501e+10]),
    2.0: numpy.array([1.83310605e+14, 3.85505323e+13, 1.06611782e+13, 4.16283619e+12,
                      1.99847087e+12, 1.09607519e+12, 6.59377398e+11, 4.24509116e+11,
                      2.87849644e+11, 2.03340410e+11]),
    5.0: numpy.array([2.52171901e+14, 8.82782026e+13, 2.61922597e+13, 1.03522489e+13,
                      4.98580852e+12, 2.73753882e+12, 1.64760874e+12, 1.06096595e+12,
                      7.19497221e+11, 5.08293436e+11]),
    10.0: numpy.array([2.88268239e+14, 1.54868128e+14, 5.09176979e+13, 2.05242390e+13,
                       9.93724544e+12, 5.46626983e+12, 3.29243871e+12, 2.12090992e+12,
                       1.43857168e+12, 1.01639496e+12])
}


def test__calc():
    """ test ratefit.calc.lindemann
    """

    lind_ktps = ratefit.calc.lindemann(HIGHP_KTS, LOWP_KTS, TEMPS, PRESSURES)

    assert numpy.allclose(tuple(lind_ktps.keys()), PRESSURES)

    assert numpy.allclose(lind_ktps[0.1], REF_LIND_KTPS[0.1], atol=0.01)
    assert numpy.allclose(lind_ktps[2.0], REF_LIND_KTPS[2.0], atol=0.01)
    assert numpy.allclose(lind_ktps[5.0], REF_LIND_KTPS[5.0], atol=0.01)
    assert numpy.allclose(lind_ktps[10.0], REF_LIND_KTPS[10.0], atol=0.01)


if __name__ == '__main__':
    test__calc()
