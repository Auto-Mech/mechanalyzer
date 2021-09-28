from ratefit.fit import arr
from autoreact.params import RxnParams

def get_params(ktp_dct, dbl_tol=15, dbl_iter=30):

    # Get the pressures
    pressures = [pressure for pressure in ktp_dct
                 if pressure != 'high']

    # Run the Arrhenius fitter for each pressure
    plog_dct = {}
    for pressure in pressures:
        fake_ktp_dct = {}
        fake_ktp_dct['high'] = ktp_dct[pressure]
        fake_params = arr.get_params(fake_ktp_dct, dbl_tol=dbl_tol,
                                     dbl_iter=dbl_iter)
        arr_params = fake_params.arr  # get the Arrhenius parameters
        plog_dct[pressure] = arr_params  # store 

    params = RxnParams(plog_dct=plog_dct)  # instantiate RxnParams

    return params
