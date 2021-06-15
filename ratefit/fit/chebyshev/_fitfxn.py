""" fit rate constants to Arrhenius expressions
"""

import numpy as np
from scipy.special import eval_chebyt


RCOND = -1
# RCOND = -1 if int(np.__version__.split('.')[1]) < 14 else None, break docs


def reaction(ktp_dct, temps, tdeg=6, pdeg=4, a_conv_factor=1.0):
    """ Fits T,P-dependent rate constants [k(T,P)]s to a
        a Chebyshev functional expression.

        :param ktp_dct: k(T,P)s
        :type ktp_dct: k(T,P) dictionary
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: np.ndarray
        :param tdeg: degree of the temp component of Chebyshev polynomial
        :type tdeg: int
        :param pdeg: degree of the pressure component of Chebyshev polynomial
        :type pdeg: int
        :rtype: (np.ndarray, tuple(float), tuple(float))
    """

    # Get the pressures from the ktp_dct
    pressures = tuple(pressure for pressure in ktp_dct.keys()
                      if pressure != 'high')

    # Determine the number and range of temperatures and pressures
    tnum, pnum = len(temps), len(pressures)
    tmin, tmax = min(temps), max(temps)
    pmin, pmax = min(pressures), max(pressures)

    # Calculate the reduced temperatures and pressures
    tred = []
    for temp in temps:
        tred.append(
            (2.0 * temp**(-1) - tmin**(-1) - tmax**(-1)) /
            (tmax**(-1) - tmin**(-1))
        )

    pred = []
    for pressure in pressures:
        pred.append(
            (2.0 * np.log10(pressure) - np.log10(pmin) - np.log10(pmax)) /
            (np.log10(pmax) - np.log10(pmin))
        )

    # Build a numpy array for the fits
    ktps = conv_dct_to_array(ktp_dct, temps, a_conv_factor=a_conv_factor)

    # Create matrices for fits
    amat = np.zeros((tnum * pnum, tdeg * pdeg), np.float64)
    bvec = np.zeros((tnum * pnum), np.float64)
    nzero = 0
    for tidx1, temp in enumerate(tred):
        for pidx1, press in enumerate(pred):
            for tidx2 in range(tdeg):
                for pidx2 in range(pdeg):
                    idx1 = (pidx1 * tnum + tidx1)
                    idx2 = (pidx2 * tdeg + tidx2)
                    amat[idx1, idx2] = (
                        eval_chebyt(tidx2, temp) * eval_chebyt(pidx2, press)
                    )
            if ktps[tidx1, pidx1] is not None:
                bvec[idx1] = np.log10(ktps[tidx1, pidx1])
            else:
                bvec[idx1] = None
                nzero += 1

    nnonzero = tnum * pnum - nzero
    idxp = -1
    amatp = np.zeros((nnonzero, tdeg * pdeg), np.float64)
    bvecp = np.zeros((nnonzero), np.float64)
    for idx in range(tnum*pnum):
        if not np.isnan(bvec[idx]):
            idxp += 1
            bvecp[idxp] = bvec[idx]
            for idx2 in range(tdeg*pdeg):
                amatp[idxp, idx2] = amat[idx, idx2]

    # Perform least-squares fit to get alpha coefficients
    theta = np.linalg.lstsq(amatp, bvecp, rcond=RCOND)[0]

    alpha = np.zeros((tdeg, pdeg), np.float64)
    for tidx2 in range(tdeg):
        for pidx2 in range(pdeg):
            alpha[tidx2, pidx2] = theta[pidx2 * tdeg + tidx2]

    return alpha, (tmin, tmax), (pmin, pmax)


def conv_dct_to_array(ktp_dct, temps, a_conv_factor=1.0):
    """ Convert a numpy
    """

    # Calculate the fitting parameters from the filtered T,k lists
    mat_rows = []
    for tk_arr in ktp_dct.values():

        # Set the temperatures and rate constants
        k_temps = list(tk_arr[0])
        k_rate_constants = list(tk_arr[1])

        mat_row = []
        for temp in temps:
            if temp in k_temps:
                idx = k_temps.index(temp)
                val = k_rate_constants[idx]
            else:
                val = None
            mat_row.append(val)

        mat_rows.append(mat_row)

    # Build array and invert from P,T to T,P
    pt_array = np.array(mat_rows)
    tp_array = np.transpose(pt_array)

    # Convert units
    tp_array = tp_array * a_conv_factor

    return tp_array
