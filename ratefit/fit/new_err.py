import random
from statistics import mean
from mechanalyzer.calculator import rates

def get_err_dct(ref_ktp_dct, params):

    # Get pressures and temps and calculate a ktp_dct using the fit params
    pressures = [pressure
                 for pressure in ref_ktp_dct.keys()
                 if pressure != 'high']  # don't pass 'high'
#    unique_temps = []
#    for temps, _ in ref_ktp_dct.values():
#        for temp in temps:
#            if temps not in unique_temps:
#                unique_temps.append(temp)
    temps = random.choice(list(ref_ktp_dct.values()))[0]  # get a random set
    fit_ktp_dct = rates.eval_params(params, pressures, temps)

    # Calculate the err_dct
    err_dct = {}
    for pressure, (temps, _) in fit_ktp_dct.items():
        fit_kts = fit_ktp_dct[pressure][1]
        ref_kts = ref_ktp_dct[pressure][1]
        errs = 100 * (fit_kts - ref_kts) / ref_kts
        err_dct[pressure] = (temps, errs)

    return err_dct


def get_max_err(err_dct):
    """ Get the singular max error from an err_dct

    """

    max_err = 0
    for pressure, (temps, errs) in err_dct.items():
        if max(abs(errs)) > max_err:
            max_err = max(abs(errs))

    return max_err


def get_max_errs(err_dct):
    """ Get the max error for each pressure in an err_dct

    """
    
    max_errs = []
    for pressure, (temps, errs) in err_dct.items():
        max_errs.append(max(abs(errs)))

    return max_errs


def get_mean_errs(err_dct):
    """ Get the mean error for each pressure in an err_dct

    """
    
    mean_errs = []
    for pressure, (temps, errs) in err_dct.items():
        mean_errs.append(mean(abs(errs)))

    return mean_errs

