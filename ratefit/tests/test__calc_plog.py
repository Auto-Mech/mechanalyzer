"""
Test the ratefit rate constant calculators
"""

import numpy
import ratefit

# Arrhenius Parameters at each pressure for single and double fits
# (cm3/s, None, cal/mol)
PLOG_SGL_PARAM_DCT = {
    0.03: (2.88300e+15, -2.583, 1244.0),
    0.10: (1.14500e+16, -2.602, 1498.0),
    0.30: (3.31100e+16, -2.594, 1702.0),
    1.00: (1.30500e+17, -2.611, 1980.0),
    3.00: (3.87400e+17, -2.606, 2213.0),
    10.0: (1.52900e+18, -2.623, 2521.0),
    30.0: (4.66900e+18, -2.621, 2789.0),
    100.: (1.72500e+19, -2.630, 3128.0)
}
PLOG_DBL_PARAM_DCT = {
    0.03: (9.91035e+15, -2.79795, 1679.0, 7.89749e+13, -2.25289, 5770.0),
    0.10: (3.93595e+16, -2.81695, 1933.0, 3.13654e+14, -2.27190, 831.0),
    0.30: (1.13816e+17, -2.80895, 2137.0, 9.06993e+14, -2.26390, 1035.0),
    1.00: (4.48596e+17, -2.82595, 2415.0, 3.57483e+15, -2.28090, 1313.0),
    3.00: (1.33169e+18, -2.82095, 2648.0, 1.06122e+16, -2.27590, 1546.0),
    10.0: (5.25597e+18, -2.83795, 2956.0, 4.18845e+16, -2.29290, 1854.0),
    30.0: (1.60498e+19, -2.83595, 3224.0, 1.27900e+17, -2.29090, 2122.0),
    100.: (5.92972e+19, -2.84495, 3563.0, 4.72536e+17, -2.29990, 2461.0)
}

# Temperature Range and Parameter Setting
PRESSURES = numpy.array([0.1, 0.5, 1.0, 5.0, 10.0, 50.0])
TEMPS = numpy.arange(500.0, 2750.0, 250.0)
T_REF = 1.0

# Set Reference Data (K(T)s from 500-2500 K, 250 K Steps)
REF_SGL_PLOG_KTPS = {
    0.1: (TEMPS, numpy.array(
        [2.40615343e+08, 1.38480280e+08, 8.42217726e+07,
         5.47949766e+07, 3.77017890e+07, 2.71234897e+07,
         2.02225018e+07, 1.55211487e+07, 1.22014976e+07])),
    0.5: (TEMPS, numpy.array(
        [9.04908728e+08, 5.80383777e+08, 3.72652196e+08,
         2.50477382e+08, 1.76130070e+08, 1.28697719e+08,
         9.70804091e+07, 7.51916782e+07, 5.95418369e+07])),
    1.0: (TEMPS, numpy.array(
        [1.59646640e+09, 1.07613194e+09, 7.07767291e+08,
         4.82394612e+08, 3.42260057e+08, 2.51634151e+08,
         1.90660699e+08, 1.48163745e+08, 1.17625407e+08])),
    5.0: (TEMPS, numpy.array(
        [5.80435137e+09, 4.41628149e+09, 3.08530947e+09,
         2.18017995e+09, 1.58441261e+09, 1.18495834e+09,
         9.09376087e+08, 7.13718057e+08, 5.71106482e+08])),
    10.0: (TEMPS, numpy.array(
        [1.00717033e+10, 8.10060825e+09, 5.81372338e+09,
         4.17301575e+09, 3.06350104e+09, 2.30721639e+09,
         1.77962807e+09, 1.40205382e+09, 1.12519966e+09])),
    50.0: (TEMPS, numpy.array(
        [3.49793326e+10, 3.22781535e+10, 2.48095880e+10,
         1.85540953e+10, 1.39979466e+10, 1.07494117e+10,
         8.41296362e+09, 6.70334623e+09, 5.42842105e+09]))
}
REF_DBL_PLOG_KTPS = {
    0.1: (TEMPS, numpy.array(
        [2.40706563e+08, 1.38458087e+08, 8.42529745e+07,
         5.48066437e+07, 3.76971595e+07, 2.71142770e+07,
         2.02159862e+07, 1.55205017e+07, 1.22075183e+07])),
    0.5: (TEMPS, numpy.array(
        [9.05251773e+08, 5.80290797e+08, 3.72790280e+08,
         2.50530730e+08, 1.76108453e+08, 1.28654011e+08,
         9.70491321e+07, 7.51885441e+07, 5.95712159e+07])),
    1.0: (TEMPS, numpy.array(
        [1.59707264e+09, 1.07596027e+09, 7.08030032e+08,
         4.82497685e+08, 3.42218280e+08, 2.51548859e+08,
         1.90599399e+08, 1.48157667e+08, 1.17683522e+08])),
    5.0: (TEMPS, numpy.array(
        [5.80655635e+09, 4.41557705e+09, 3.08645480e+09,
         2.18064579e+09, 1.58421925e+09, 1.18455675e+09,
         9.09083781e+08, 7.13688853e+08, 5.71388725e+08])),
    10.0: (TEMPS, numpy.array(
        [1.00755441e+10, 8.09932871e+09, 5.81589069e+09,
         4.17391394e+09, 3.06313192e+09, 2.30643801e+09,
         1.77905872e+09, 1.40199855e+09, 1.12575740e+09])),
    50.0: (TEMPS, numpy.array(
        [3.49926908e+10, 3.22730713e+10, 2.48188492e+10,
         1.85580981e+10, 1.39962671e+10, 1.07457907e+10,
         8.41027650e+09, 6.70308549e+09, 5.43111470e+09]))
}


def test__plog_from_single_arrhenius():
    """ test ratefit.calc.plog
    """

    plog_ktps = ratefit.calc.plog(PLOG_SGL_PARAM_DCT, T_REF, TEMPS, PRESSURES)

    assert numpy.allclose(tuple(plog_ktps.keys()), PRESSURES)

    assert numpy.allclose(
        plog_ktps[0.1][1], REF_SGL_PLOG_KTPS[0.1][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[0.5][1], REF_SGL_PLOG_KTPS[0.5][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[1.0][1], REF_SGL_PLOG_KTPS[1.0][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[5.0][1], REF_SGL_PLOG_KTPS[5.0][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[10.0][1], REF_SGL_PLOG_KTPS[10.0][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[50.0][1], REF_SGL_PLOG_KTPS[50.0][1], atol=0.001)


def test__plog_from_double_arrhenius():
    """ test ratefit.calc.plog
    """

    plog_ktps = ratefit.calc.plog(PLOG_DBL_PARAM_DCT, T_REF, TEMPS, PRESSURES)

    assert set(PRESSURES) <= set(plog_ktps.keys())
    assert numpy.allclose(
        plog_ktps[0.1][1], REF_DBL_PLOG_KTPS[0.1][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[0.5][1], REF_DBL_PLOG_KTPS[0.5][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[1.0][1], REF_DBL_PLOG_KTPS[1.0][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[5.0][1], REF_DBL_PLOG_KTPS[5.0][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[10.0][1], REF_DBL_PLOG_KTPS[10.0][1], atol=0.001)
    assert numpy.allclose(
        plog_ktps[50.0][1], REF_DBL_PLOG_KTPS[50.0][1], atol=0.001)
