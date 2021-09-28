""" This is Clayton's new version
"""

import copy
import numpy
from ratefit.fit import arr
from ratefit.fit import plog
from ratefit.fit import cheb

DEFAULT_PDEP = {
    'assess_temps': (500.0, 1000.0),
    'tol': 20.0,
    'plow': None,
    'phigh': None,
    'pval': 1.0,
}
DEFAULT_ARR = {  # also used for PLOG fitting
    'dbl_tol': 15.0,
    'dbl_iter': 30
}
DEFAULT_TROE = {
    'params': ('ts1', 'ts2', 'ts3', 'alpha'),
    'tol': 20.0
}
DEFAULT_CHEB = {
    'tdeg': 6,
    'pdeg': 4,
    'tol': 20.0
}
NICKNAMES = {  # used for printing messages
    'arr': 'Arrhenius',
    'plog': 'PLOG',
    'cheb': 'Chebyshev',
    'troe': 'Troe'
}


def fit_ktp_dct(ktp_dct, fit_method,
                pdep_dct=DEFAULT_PDEP,
                arrfit_dct=DEFAULT_ARR,
                chebfit_dct=DEFAULT_CHEB,
                troefit_dct=DEFAULT_TROE,
                fit_temps=None, fit_pressures=None):
    """ Fits a single ktp_dct to some desired form

        :param ktp_dct: dictionary of k(T,P) values 
        :type ktp_dct: dct {P1: (temps1, kts1), P2: ...} 
        :param fit_method: desired fit form; 'arr', 'plog', 'cheb', or 'troe'
        :type fit_method: str
        :param arrfit_dct: instructions for Arrhenius fitting (also for PLOG)
        :type arrfit_dct: dct
        :param chebfit_dct: instructions for Chebyshev fitting
        :type chebfit_dct: dct
        :param troefit_dct: instructions for Troe fitting
        :type troefit_dct: dct
        :param fit_temps: limits of allowed temperatures; k(T)s outside these 
            limits will be thrown out 
        :type fit_temps: tuple (tmin, tmax)
        :param fit_pressures: limits of allowed pressures; pressures outside
            these limits will be thrown out
        :type fit_pressures: tuple (pmin, pmax)
    """

    # Filter out undesired temperatures and pressures, if given
    ktp_dct = filter_ktp_dct(ktp_dct, fit_temps=fit_temps, 
                             fit_pressures=fit_pressures)

    # Get the pdep_version of the ktp_dct
    assess_temps = pdep_dct['assess_temps']
    tol = pdep_dct['tol']
    plow = pdep_dct['plow']
    phigh = pdep_dct['phigh']
    pval = pdep_dct['pval']
    pdep_ktp_dct = get_pdep_ktp_dct(ktp_dct, assess_temps=assess_temps,
                                    tol=tol, plow=plow, phigh=phigh, pval=pval)

    # Check the ktp_dct and the fit_method to see how to fit rates
    actual_fit_method = assess_fit_method(pdep_ktp_dct, fit_method) 

    # Get the desired fits as instances of the RxnParams class
    if actual_fit_method is None:  # occurs when the ktp_dct is empty
        params = None
    elif actual_fit_method == 'arr':
        dbl_tol = arrfit_dct['dbl_tol']
        dbl_iter = arrfit_dct['dbl_iter']
        params = arr.get_params(pdep_ktp_dct, dbl_tol=dbl_tol, 
                                dbl_iter=dbl_iter)
    elif actual_fit_method == 'plog':
        dbl_tol = arrfit_dct['dbl_tol']
        dbl_iter = arrfit_dct['dbl_iter']
        params = plog.get_params(pdep_ktp_dct, dbl_tol=dbl_tol, 
                                 dbl_iter=dbl_iter)
    elif actual_fit_method == 'cheb':
        tdeg = chebfit_dct['tdeg']
        pdeg = chebfit_dct['pdeg']
        tol = chebfit_dct['tol']
        params = cheb.get_params(pdep_ktp_dct, tdeg=tdeg, pdeg=pdeg, 
                                 tol=tol)
    elif actual_fit_method == 'troe':
        pass

    return params


def assess_fit_method(pdep_ktp_dct, fit_method):
    """ Checks a pdep_ktp_dct (i.e., already filtered for pressure dependence)
        to see if the selected fit method is acceptable. If something is amiss,
        corrects the fit method to the proper form.
    """

    if pdep_ktp_dct:
        pressures = list(pdep_ktp_dct.keys())
        npressures = len(pressures)
        if npressures == 1 or (npressures == 2 and 'high' in pressures):
            actual_fit_method = 'arr'
        else:
            # If the fit_method is Arrhenius but there are multiple pressures
            if fit_method == 'arr':
                actual_fit_method = 'plog'
            else:
                actual_fit_method = fit_method
    else:
        actual_fit_method = None

    # If Chebyshev was requested and granted, check that it is viable to do
    # Chebyshev based on the number of k(T) values at each pressure
    if actual_fit_method == 'cheb':
        # Note: this check ignores the 'high' entry
        cheb_viable = cheb.check_viability(pdep_ktp_dct)
        if not cheb_viable:
            actual_fit_method = 'plog'

    # Print message to say what fitting will be done
    if actual_fit_method != fit_method:  # if the fit method was changed
        if actual_fit_method == 'plog':  
            # If changed to PLOG from Arrhenius
            if fit_method == 'arr':
                print('\nPressure dependence detected. Fitting to PLOG form...')
            elif fit_method == 'cheb':
                print('\nRates invalid for Chebyshev. Fitting to PLOG form...')
        else:  # if changed to Arrhenius from anything else
            print(f'\nNot enough pressures for {NICKNAMES[fit_method]} form.')
            print('\nFitting k(T,P)s to Arrhenius form...')
    elif fit_method is None:  # if the ktp_dct was empty
        print('\nNo valid k(T,Ps)s to fit. Skipping to next reaction...')
    else:  # if the fit method was unchanged  
        print(f'\nFitting k(T,P)s to {NICKNAMES[fit_method]} form...')

    return actual_fit_method


def get_pdep_ktp_dct(ktp_dct, 
                     assess_temps=DEFAULT_PDEP['assess_temps'], 
                     tol=DEFAULT_PDEP['tol'],
                     plow=DEFAULT_PDEP['plow'],
                     phigh=DEFAULT_PDEP['phigh'],
                     pval=DEFAULT_PDEP['pval']):
    """ Returns a ktp_dct with either (i) pressure-dependent rates or (ii)
        pressure-independent rates, depending on what's in the ktp_dct
    """

    # Assess the pressure dependence of the rate constants
    rxn_is_pdependent = assess_pdep(ktp_dct, assess_temps, tol=tol, plow=plow,
                                    phigh=phigh)

    if rxn_is_pdependent:
        print('Reaction found to be pressure dependent.',
              'Fitting all k(T)s from all pressures.')
        pdep_ktp_dct = copy.deepcopy(ktp_dct)
    else:
        # If pval is in the dct, return kts at that pressure
        if pval in ktp_dct:
            pdep_ktp_dct = {'high': ktp_dct[pval]}
            print(f'No pressure dependence. Grabbing k(T)s at {pval} atm.')
        # Otherwise, get kts at some other pressure
        else:
            print(f'No pressure dependence, but no k(T)s at {pval} atm.')
            pressures = [pressure for pressure in ktp_dct
                         if pressure != 'high']
            if pressures == []:
                pdep_ktp_dct = {'high': ktp_dct['high']}
                print('No numerical pressures. Grabbing "high" kts.')
            else:
                pressure = pressures[-1]
                pdep_ktp_dct = {pressure: ktp_dct[pressure]}
                print(f'Grabbing kts at {pressure} atm.')

    return pdep_ktp_dct


def assess_pdep(ktp_dct, 
                assess_temps=DEFAULT_PDEP['assess_temps'], 
                tol=DEFAULT_PDEP['tol'],
                plow=DEFAULT_PDEP['plow'],
                phigh=DEFAULT_PDEP['phigh']):
    """ Assesses if there are significant differences between k(T,P) values
        at low and high pressure, signaling pressure dependence
    """

    # Get list of pressures, ignoring the high-pressure limit rates
    pressures = [pressure for pressure in ktp_dct
                 if pressure != 'high']

    # Set the low- and high-pressure if not specified by user
    if plow is None and pressures != []:
        plow = min(pressures)
    if phigh is None and pressures != []:
        phigh = max(pressures)

    # Check % difference for k(T, P) vals
    is_pdep = False
    if plow in ktp_dct and phigh in ktp_dct:  # won't get run if only 'high'

        # Loop over temps to examine for large % dif in k(T) at low- and high-P
        for assess_temp in assess_temps:
            # For low and high P, find the idx for assess_temp
            temps_low = ktp_dct[plow][0]
            temps_high = ktp_dct[phigh][0]
            temp_low_match = numpy.where(
                numpy.isclose(temps_low, assess_temp))[0]
            temp_high_match = numpy.where(
                numpy.isclose(temps_high, assess_temp))[0]
            if temp_low_match.size > 0 and temp_high_match.size > 0:
                temp_low_idx = temp_low_match[0]
                temp_high_idx = temp_high_match[0]
                # Grab k value for the appropriate temp and pressure
                kt_low = ktp_dct[plow][1][temp_low_idx]
                kt_high = ktp_dct[phigh][1][temp_high_idx]
                # Calculate the % difference and see if above threshold
                kt_dif = (abs(kt_low - kt_high) / kt_low) * 100.0
                if kt_dif > tol:
                    is_pdep = True
                    break

    return is_pdep


def filter_ktp_dct(ktp_dct, fit_temps=None, fit_pressures=None):
    """ Removes entries from a ktp_dct if they fall outside specified temp or
        pressure limits

        :param ktp_dct: dictionary of k(T,P) values 
        :type ktp_dct: dct {P1: (temps1, kts1), P2: ...} 
        :param fit_temps: limits of allowed temperatures 
        :type fit_temps: tuple (tmin, tmax)
        :param fit_pressures: limits of allowed pressures
        :type fit_pressures: tuple (pmin, pmax)
    """

    # Filter the temps if fit_temps were given
    if fit_temps is not None:
        filt_ktp_dct = {}
        tmin = fit_temps[0]
        tmax = fit_temps[1]
        for pressure, (temps, kts) in ktp_dct.items():
            filt_temps = []
            filt_kts = []
            for tidx, temp in enumerate(temps):
                if tmin <= temp <= tmax:
                    filt_temps.append(temp)
                    filt_kts.append(kts[tidx])
            filt_ktp_dct[pressure] = (filt_temps, filt_kts)
    # Otherwise, just copy the input ktp_dct
    else:
        filt_ktp_dct = copy.deepcopy(ktp_dct)

    # Remove undesired pressures if fit_pressure were given (ignores 'high')
    if fit_pressures is not None:
        pmin = fit_pressures[0]
        pmax = fit_pressures[1]
        pressures = tuple(pressure for pressure in filt_ktp_dct.keys()
                          if pressure != 'high')
        for pressure in pressures:
            if pmin <= pressure <= pmax:
               pass
            else:
                filt_ktp_dct.pop(pressure)

    return filt_ktp_dct

