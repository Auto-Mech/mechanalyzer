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

# Set Troe Parameters used to Modify the Lindeman calcs
TROE_ALPHA, TROE_T3, TROE_T1, TROE_T2 = 6.0e-1, 1.0e3, 7.0, 1.7e3

# Pressure and Temperature Range and Parameter Setting
PRESSURES = numpy.array([0.1, 2.0, 5.0, 10.0])
TEMPS = numpy.arange(300.0, 3300.0, 300.0)
T_REF = 1.0

# Set Reference Data for Comparison
REF_TROE_3PAR_KTPS = {
    0.1: numpy.array(
        [9.01050498e+12, 1.14424645e+12, 2.94814049e+11, 1.06633311e+11,
         4.65907324e+10, 2.29034703e+10, 1.21838431e+10, 6.85018641e+09,
         4.00860012e+09, 2.41596987e+09]),
    2.0: numpy.array(
        [5.43538329e+13, 1.33486434e+13, 3.89709460e+12, 1.44805772e+12,
         6.30903451e+11, 3.05193075e+11, 1.58593498e+11, 8.67036892e+10,
         4.91838285e+10, 2.86739835e+10]),
    5.0: numpy.array(
        [8.34523647e+13, 2.48240318e+13, 7.90510043e+12, 3.01068770e+12,
         1.31964256e+12, 6.37327836e+11, 3.29377943e+11, 1.78696748e+11,
         1.00461305e+11, 5.79982923e+10]),
    10.0: numpy.array(
        [1.13703761e+14, 3.77000554e+13, 1.30243198e+13, 5.09444837e+12,
         2.25136399e+12, 1.08810396e+12, 5.60704044e+11, 3.02702223e+11,
         1.69142296e+11, 9.69910174e+10])
}
REF_TROE_4PAR_KTPS = {
    0.1: numpy.array(
        [9.10804237e+12, 1.32381727e+12, 4.08085808e+11, 1.72813491e+11,
         8.71750733e+10, 4.93305843e+10, 3.03035430e+10, 1.97990317e+10,
         1.35713589e+10, 9.66617670e+09]),
    2.0: numpy.array(
        [5.49814978e+13, 1.69633215e+13, 6.60737791e+12, 3.07653037e+12,
         1.62249315e+12, 9.39965627e+11, 5.85382158e+11, 3.85776465e+11,
         2.65960784e+11, 1.90195810e+11]),
    5.0: numpy.array(
        [8.42788325e+13, 3.22567471e+13, 1.45617858e+13, 7.20296797e+12,
         3.90278656e+12, 2.29203205e+12, 1.43833189e+12, 9.52289745e+11,
         6.58503695e+11, 4.71883493e+11]),
    10.0: numpy.array(
        [1.14711507e+14, 4.91180769e+13, 2.54925292e+13, 1.34565982e+13,
         7.50640433e+12, 4.47219704e+12, 2.82864928e+12, 1.88160867e+12,
         1.30503129e+12, 9.37076153e+11])
}


def test__calc_three_param():
    """ test ratefit.calc._rates.troe
    """

    troe_3par_ktps = ratefit.calc.troe(
        HIGHP_KTS, LOWP_KTS, TEMPS, PRESSURES,
        TROE_ALPHA, TROE_T3, TROE_T1, ts2=None)

    assert numpy.allclose(tuple(troe_3par_ktps.keys()), PRESSURES)

    assert numpy.allclose(troe_3par_ktps[0.1], REF_TROE_3PAR_KTPS[0.1], atol=0.01)
    assert numpy.allclose(troe_3par_ktps[2.0], REF_TROE_3PAR_KTPS[2.0], atol=0.01)
    assert numpy.allclose(troe_3par_ktps[5.0], REF_TROE_3PAR_KTPS[5.0], atol=0.01)
    assert numpy.allclose(troe_3par_ktps[10.0], REF_TROE_3PAR_KTPS[10.0], atol=0.01)


def test__calc_four_param():
    """ test ratefit.calc._rates.troe
    """

    troe_4par_ktps = ratefit.calc.troe(
        HIGHP_KTS, LOWP_KTS, TEMPS, PRESSURES,
        TROE_ALPHA, TROE_T3, TROE_T1, ts2=TROE_T2)

    assert numpy.allclose(tuple(troe_4par_ktps.keys()), PRESSURES)

    assert numpy.allclose(troe_4par_ktps[0.1], REF_TROE_4PAR_KTPS[0.1], atol=0.01)
    assert numpy.allclose(troe_4par_ktps[2.0], REF_TROE_4PAR_KTPS[2.0], atol=0.01)
    assert numpy.allclose(troe_4par_ktps[5.0], REF_TROE_4PAR_KTPS[5.0], atol=0.01)
    assert numpy.allclose(troe_4par_ktps[10.0], REF_TROE_4PAR_KTPS[10.0], atol=0.01)


if __name__ == '__main__':
    test__calc_three_param()
    test__calc_four_param()
