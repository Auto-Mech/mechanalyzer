""" Fits rate constants to a PLOG form
"""

from ratefit.fit import arr
from autoreact.params import RxnParams

def get_params(ktp_dct, dbltol=15, dbl_iter=1, tref=1.0):
    """ Gets the fitting parameters for a PLOG fit to rate constant data.
        Also gets the errors of that fit. Performs either a single or double
        Arrhenius fit at each pressure

        :param ktp_dct: rate constants to be fitted; should be P-dependent
        :type ktp_dct: dict {pressure: (temps, kts)}
        :param dbltol: the tolerance at which to reject a single fit
        :type dbltol: float
        :param dbl_iter: max number of iterations for double fitting
        :type dbl_iter: int
        :param tref: reference temp for the modified Arrhenius form
        :type tref: float
        :return params: fitted Arrhenius parameters
        :rtype: autoreact.RxnParams object
        :return err_dct: fitting errors
        :rtype: dict {pressure: (temps, errs)}
    """

    # Get the pressures
    pressures = [pressure for pressure in ktp_dct
                 if pressure != 'high']

    # Run the Arrhenius fitter for each pressure
    plog_dct = {}
    err_dct = {}
    for pressure in pressures:
        fake_ktp_dct = {}
        fake_ktp_dct['high'] = ktp_dct[pressure]  # store with 'high' for arr
        fake_params, fake_err_dct = arr.get_params(
            fake_ktp_dct, dbltol=dbltol, dbl_iter=dbl_iter, tref=tref)
        arr_params = fake_params.arr  # get the Arrhenius parameters
        plog_dct[pressure] = arr_params  # store in the plog_dct
        err_dct[pressure] = fake_err_dct['high']

    params = RxnParams(plog_dct=plog_dct)  # instantiate RxnParams

    return params, err_dct
