""" test ratefit.calc.lowp_limit
"""

import numpy
import ratefit


# K(T) at low-pressure limits (T = [300, 3000]K, 300 K Steps)
LOWP_KTS = numpy.array(
    [4.957573597808111e+18, 1.010824653785529e+18, 3.983435704402545e+17,
     2.056812288924121e+17, 1.231626977413509e+17, 8.099919515048618e+16,
     5.683130672639794e+16, 4.180898071689080e+16, 3.189099736360676e+16,
     2.503025834637141e+16])

# A_HIGH, N_HIGH, EA_HIGH = 2.000e+12, 0.900, 4.87490
# A_LOW, N_LOW, EA_LOW = 2.490e24, -2.300, 4.87490

# Temperature Range and Parameter Setting
PRESSURES = numpy.array([0.5, 1.0, 2.0, 5.0, 10.0])
TEMPS = numpy.arange(300.0, 3300.0, 300.0)

REF_LOWP_KTPS = {
    0.5: numpy.array(
        [1.00693247e+14, 1.02654267e+13, 2.69691795e+12, 1.04439755e+12,
         5.00311362e+11, 2.74195694e+11, 1.64900047e+11, 1.06147744e+11,
         7.19708728e+10, 5.08389426e+10]),
    1.0: numpy.array(
        [2.01386495e+14, 2.05308534e+13, 5.39383590e+12, 2.08879510e+12,
         1.00062272e+12, 5.48391388e+11, 3.29800094e+11, 2.12295489e+11,
         1.43941746e+11, 1.01677885e+11]),
    2.0: numpy.array(
        [4.02772990e+14, 4.10617069e+13, 1.07876718e+13, 4.17759020e+12,
         2.00124545e+12, 1.09678278e+12, 6.59600187e+11, 4.24590977e+11,
         2.87883491e+11, 2.03355770e+11]),
    5.0: numpy.array(
        [1.00693247e+15, 1.02654267e+14, 2.69691795e+13, 1.04439755e+13,
         5.00311362e+12, 2.74195694e+12, 1.64900047e+12, 1.06147744e+12,
         7.19708728e+11, 5.08389426e+11]),
    10.0: numpy.array(
        [2.01386495e+15, 2.05308534e+14, 5.39383590e+13, 2.08879510e+13,
         1.00062272e+13, 5.48391388e+12, 3.29800094e+12, 2.12295489e+12,
         1.43941746e+12, 1.01677885e+12])
}


def test__calc():
    """ test ratefit.calc.rate.lowp_limit
        test ratefit.calc.rate.lowp_limit_one_pressure
    """

    lowp_ktps = ratefit.calc.lowp_limit(LOWP_KTS, TEMPS, PRESSURES)

    assert numpy.allclose(tuple(lowp_ktps.keys()), PRESSURES)

    assert numpy.allclose(lowp_ktps[0.5], REF_LOWP_KTPS[0.5], atol=0.01)
    assert numpy.allclose(lowp_ktps[1.0], REF_LOWP_KTPS[1.0], atol=0.01)
    assert numpy.allclose(lowp_ktps[2.0], REF_LOWP_KTPS[2.0], atol=0.01)
    assert numpy.allclose(lowp_ktps[5.0], REF_LOWP_KTPS[5.0], atol=0.01)
    assert numpy.allclose(lowp_ktps[10.0], REF_LOWP_KTPS[10.0], atol=0.01)
