"""
Test the ratefit rate constant calculators
"""

import numpy
import pandas
import ratefit


# Arrhenius Parameters for single and double fits (cm3/s, None, cal/mol)
SGL_A1, SGL_N1, SGL_EA1 = 2.33400e+07, 2.084, 11104

DBL_A1, DBL_N1, DBL_EA1 = 8.03039e+07, 1.869, 11539
DBL_A2, DBL_N2, DBL_EA2 = 6.38505e+05, 2.414, 10434

# Temperature Range and Parameter Setting
TEMPS = numpy.arange(500.0, 2600.0, 100.0)
T_REF = 1.0

# Set Reference Data for Comparison
DATA = numpy.array([
    [1.378150410638410E+08, 1.379756660315570E+08],
    [1.297836416161970E+09, 1.297804739145440E+09],
    [6.769080072288560E+09, 6.769503762000790E+09],
    [2.425098083508950E+10, 2.425882930117080E+10],
    [6.735608561972120E+10, 6.738903275406440E+10],
    [1.560888894758010E+11, 1.561706159302750E+11],
    [3.164069841563090E+11, 3.165505453448880E+11],
    [5.792159567474590E+11, 5.793995414948790E+11],
    [9.791388652530330E+11, 9.792822687500520E+11],
    [1.553294914048580E+12, 1.553245242125170E+12],
    [2.340215195297280E+12, 2.339748206448680E+12],
    [3.378946826539680E+12, 3.377787484338700E+12],
    [4.708351078718470E+12, 4.706214302634070E+12],
    [6.366577818492750E+12, 6.363228486366410E+12],
    [8.390689299332350E+12, 8.386019320369610E+12],
    [1.081640479955560E+13, 1.081052443247610E+13],
    [1.367793986577610E+13, 1.367128024647480E+13],
    [1.700791768884750E+13, 1.700134160880170E+13],
    [2.083733419920020E+13, 2.083225242294330E+13],
    [2.519556223651740E+13, 2.519405297926480E+13],
    [3.011038339418150E+13, 3.011531294245150E+13]
])
REF_ARR_KTS = pandas.DataFrame(
    data=DATA, index=TEMPS, columns=('Single', 'Double'))


def test__single_arrhenius():
    """ test ratefit.calc._rates.arrhenius
        test ratefit.calc._rates.single_arrhenius
    """

    calc_ks1 = ratefit.calc.arrhenius(
        ((SGL_A1, SGL_N1, SGL_EA1),), T_REF, TEMPS)
    calc_ks2 = ratefit.calc.single_arrhenius(
        SGL_A1, SGL_N1, SGL_EA1, T_REF, TEMPS)

    assert numpy.allclose(calc_ks1, REF_ARR_KTS['Single'], atol=0.01)
    assert numpy.allclose(calc_ks2, REF_ARR_KTS['Single'], atol=0.01)


def test__double_arrhenius():
    """ test ratefit.calc._rates.arrhenius
        test ratefit.calc._rates.double_arrhenius
    """

    calc_ks1 = ratefit.calc.arrhenius(
        ((DBL_A1, DBL_N1, DBL_EA1), (DBL_A2, DBL_N2, DBL_EA2)), T_REF, TEMPS)
    calc_ks2 = ratefit.calc.double_arrhenius(
        DBL_A1, DBL_N1, DBL_EA1, DBL_A2, DBL_N2, DBL_EA2, T_REF, TEMPS)

    assert numpy.allclose(calc_ks1, REF_ARR_KTS['Double'], atol=0.01)
    assert numpy.allclose(calc_ks2, REF_ARR_KTS['Double'], atol=0.01)
