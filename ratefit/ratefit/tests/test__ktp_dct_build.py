""" ratefit.fit._util
"""

import numpy
import ratefit


# Dicts with fake rates containing several different situations for k(T)s
# that can be obtained from master equation output

# Handles negatives and undefined rates
KTP_DCT1 = {
    # Weird negative, non-physical k(T)s and undefined k(T)s
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (-1.0e3, 2.0e3, -3.0e3, None, None, None)),
    # Just negative k(T) at first position and undefined k(T)s
    0.1: ((200., 400., 600., 800., 1000., 1200.),
          (-1.0e4, 2.0e4, 3.0e4, None, None, None)),
    # Positive k(T) and undefined k(T)s
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (1.0e5, 2.0e5, 3.0e5, None, None, None)),
    # Positive k(T) and undefined k(T)s
    10.0: ((200., 400., 600., 800., 1000., 1200.),
           (1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6)),
    # Positive k(T) and undefined k(T)s
    100.0: ((200., 400., 600., 800., 1000., 1200.),
            (1.0e7, 2.0e7, 3.0e7, 4.0e7, 5.0e7, 6.0e7))
}

# Handles very small rate constants that may not wish to be modelled
KTP_DCT2 = {
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (1.0e-25, 5.0e-25, 9.0e-24, 1.5e-23, 6.5e-22, 1.0e-20)),
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (8.0e-25, 3.0e-24, 2.0e-22, 6.0e-20, 2.5e-18, 5.0e-16)),
}

# Dict that should come back empty after filter
KTP_DCT3 = {
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (-1.0e-25, -5.0e-25, -9.0e-24, -1.5e-23, -6.5e-22, -1.0e-20)),
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (-8.0e-25, -3.0e-24, -2.0e-22, -6.0e-20, -2.5e-18, -5.0e-16)),
}

# Normal k(T,P) dict used for other functions
KTP_DCT4 = {
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (1.0e+5, 2.0e+5, 3.0e+5, 4.0e+5, 5.0e+5, 6.0e+5)),
    10.0: ((200., 400., 600., 800., 1000., 1200.),
           (1.0e+6, 2.0e+6, 3.0e+6, 4.0e+6, 5.0e+6, 6.0e+6)),
    'high': ((200., 400., 600., 800., 1000., 1200.),
             (1.0e+10, 2.0e+10, 3.0e+10, 4.0e+10, 5.0e+10, 6.0e+10))
}


def test__filter_ktp():
    """ test ratefit.fit.util.get_valid_tk for no rates
    """

    # Check the temperatures and undefined (these are wrong)
    ref_dct1 = {
        1.0: (numpy.array([200., 400., 600.]),
              numpy.array([1.0e5, 2.0e5, 3.0e5])),
        10.0: (numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
               numpy.array([1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6])),
        100.0: (numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
                numpy.array([1.0e7, 2.0e7, 3.0e7, 4.0e7, 5.0e7, 6.0e7]))
    }

    ref_dct2 = {
        0.01: (numpy.array([400.]),
               numpy.array([2.0e3])),
        0.1: (numpy.array([400., 600.]),
              numpy.array([2.0e5, 3.0e5])),
        1.0: (numpy.array([400., 600.]),
              numpy.array([2.0e5, 3.0e5])),
        10.0: (numpy.array([400., 600., 800.]),
               numpy.array([2.0e5, 3.0e5, 4.0e5])),
        100.0: (numpy.array([400., 600., 800.]),
                numpy.array([2.0e5, 3.0e5, 4.0e5]))}

    filt_ktp_dct_temp1 = ratefit.fit.filter_ktp_dct(
        KTP_DCT1, bimol=False, tmin=None, tmax=None)

    filt_ktp_dct_temp2 = ratefit.fit.filter_ktp_dct(
        KTP_DCT1, bimol=False, tmin=400.0, tmax=800.0)
    print(filt_ktp_dct_temp1)
    print(filt_ktp_dct_temp2)
    import sys
    sys.exit()
    assert filt_ktp_dct_temp1 == ref_dct1
    assert filt_ktp_dct_temp2 == ref_dct2

    # Check for small rate constants (based on bimol flag)
    ref_dct3 = {
        0.01: (numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
               numpy.array([1.0e-25, 5.0e-25, 9.0e-24, 1.5e-23, 6.5e-22, 1.0e-20])),
        1.0: (numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
              numpy.array([8.0e-25, 3.0e-24, 2.0e-22, 6.0e-20, 2.5e-18, 5.0e-16]))}

    ref_dct4 = {
        0.01: (numpy.array([600.,  800., 1000., 1200.]),
               numpy.array([9.0e-24, 1.5e-23, 6.5e-22, 1.0e-20])),
        1.0: (numpy.array([400.,  600.,  800., 1000., 1200.]),
              numpy.array([3.0e-24, 2.0e-22, 6.0e-20, 2.5e-18, 5.0e-16]))}

    filt_ktp_dct_smallk1 = ratefit.fit.filter_ktp_dct(
        KTP_DCT2, bimol=False, tmin=None, tmax=None)
    filt_ktp_dct_smallk2 = ratefit.fit.filter_ktp_dct(
        KTP_DCT2, bimol=True, tmin=None, tmax=None)

    assert filt_ktp_dct_smallk1 == ref_dct3
    assert filt_ktp_dct_smallk2 == ref_dct4

    # Dict that should come back empty
    ref_dct5 = {}

    filt_ktp_dct_emptyret = ratefit.fit.filter_ktp_dct(
        KTP_DCT3, bimol=False, tmin=None, tmax=None)

    assert filt_ktp_dct_emptyret == ref_dct5


def test__invert():
    """ test ratefit.fit.invert_ktp_dct
    """

    ref_inv_ktp_dct = {
        200.: ((1.0, 10.0, 'high'), (1.0e+5, 1.0e+6, 1.0e+10)),
        400.: ((1.0, 10.0, 'high'), (2.0e+5, 2.0e+6, 2.0e+10)),
        600.: ((1.0, 10.0, 'high'), (3.0e+5, 3.0e+6, 3.0e+10)),
        800.: ((1.0, 10.0, 'high'), (4.0e+5, 4.0e+6, 4.0e+10)),
        1000.0: ((1.0, 10.0, 'high'), (5.0e+5, 5.0e+6, 5.0e+10)),
        1200.0: ((1.0, 10.0, 'high'), (6.0e+5, 6.0e+6, 6.0e+10))
    }

    inv_ktp_dct = ratefit.fit.invert_ktp_dct(KTP_DCT4)

    assert inv_ktp_dct == ref_inv_ktp_dct


def test__pull_highp():
    """
    """



if __name__ == '__main__':
    test__filter_ktp()
    test__invert()
