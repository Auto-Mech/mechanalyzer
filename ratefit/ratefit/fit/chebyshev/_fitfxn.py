""" fit rate constants to Arrhenius expressions
"""

import numpy as np
from scipy.special import eval_chebyt


# RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K
RCOND = -1 if int(np.__version__.split('.')[1]) < 14 else None


def reaction(temps, ktp_dct, tdeg=6, pdeg=4, a_conv_factor=1):
    """ Fits T,P-dependent rate constants [k(T,P)]s to a
        a Chebyshev functional expression.

        :param alpha: Chebyshev coefficient matrix
        :type alpha: np.ndarray
        :param tdeg: degree of the temp component of Chebyshev polynomial
        :type tdeg: int
        :param pdeg: degree of the pressure component of Chebyshev polynomial
        :type pdeg: int
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: np.ndarray
        :param pressure: Pressure used to calculate k(T,P)s
        :type pressure: float
        :return alpha: alpha fitting coefficients
        :rtype np.ndarray
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
    A = np.zeros((tnum * pnum, tdeg * pdeg), np.float64)
    b = np.zeros((tnum * pnum), np.float64)
    nzero = 0
    for t1, T in enumerate(tred):
        for p1, P in enumerate(pred):
            for t2 in range(tdeg):
                for p2 in range(pdeg):
                    idx1, idx2 = (p1 * tnum + t1), (p2 * tdeg + t2)
                    A[idx1, idx2] = eval_chebyt(t2, T) * eval_chebyt(p2, P)
            if ktps[t1, p1] is not None:
                b[idx1] = np.log10(ktps[t1, p1])
            else:
                b[idx1] = None
                nzero += 1

    nnonzero = tnum * pnum - nzero
    idxp = -1
    Ap = np.zeros((nnonzero, tdeg * pdeg), np.float64)
    bp = np.zeros((nnonzero), np.float64)
    for idx in range(tnum*pnum):
        if not np.isnan(b[idx]):
            idxp += 1
            bp[idxp] = b[idx]
            for idx2 in range(tdeg*pdeg):
                Ap[idxp,idx2] = A[idx,idx2]

        
    # Perform least-squares fit to get alpha coefficients
    theta = np.linalg.lstsq(Ap, bp, rcond=RCOND)[0]

    alpha = np.zeros((tdeg, pdeg), np.float64)
    for t2 in range(tdeg):
        for p2 in range(pdeg):
            alpha[t2, p2] = theta[p2 * tdeg + t2]

    return alpha, (tmin, tmax), (pmin, pmax)


def conv_dct_to_array(ktp_dct, temps, a_conv_factor=1):
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
