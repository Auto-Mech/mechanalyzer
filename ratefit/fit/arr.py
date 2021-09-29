import numpy
import random
from scipy.optimize import leastsq
from mechanalyzer.calculator import rates
from autoreact.params import RxnParams
from phydat import phycon
from ratefit.fit import new_err
# import pyswarms as ps

RC = phycon.RC_CAL  # universal gas constant in cal/mol-K
GUESS_BNDS = ((1, 1e4), (0.1, 20), (1, 100))  # guess bounds for double fitting


def get_params(ktp_dct, dbl_tol=15, dbl_iter=30, t_ref=1.0):

    assert len(ktp_dct) == 1 and ktp_dct.get('high') is not None, (
        "There should only be one pressure, 'high' in the ktp_dct")
    
    # Perform single fit and assess its errors
    (temps, kts) = ktp_dct['high']  # read data
    sing_params = single_arr(temps, kts)
    sing_err_dct = new_err.get_err_dct(ktp_dct, sing_params) 
    sing_max_err = new_err.get_max_err(sing_err_dct)

    # If single fit has acceptable errors, use the single fit
    if sing_max_err < dbl_tol:
        print(f'Single fit error is {sing_max_err:.1f}%, which is less than'
              f' the input limit of {dbl_tol}%. Using single fit.')
        params = sing_params

    # Otherwise, perform a double fit
    else:
        print(f'Single fit error is {sing_max_err:.1f}%, which is more than' 
              f' the input limit of {dbl_tol}%. Attempting double fit...')
        # This call is for the iterative version of double_arr
        doub_params, guess_idx = double_arr_iter(temps, kts, sing_params, 
                                            dbl_tol=dbl_tol, dbl_iter=dbl_iter) 

        # This call is for the original, simple form of double_arr
        #doub_params= double_arr_orig(temps, kts, sing_params, 
        #                        dbl_tol=dbl_tol, dbl_iter=dbl_iter) 
        #guess_idx = 0

        # This call is for the PSO version of double_arr
        #doub_params = double_arr_pso(temps, kts, t_ref)
        #guess_idx = 0

        # Assess errors
        doub_err_dct = new_err.get_err_dct(ktp_dct, doub_params) 
        doub_max_err = new_err.get_max_err(doub_err_dct)
        print(f'Double fit obtained with max error of {doub_max_err:.1f}% after'
              f' {guess_idx + 1} iteration(s).')

        # Use single fit if double fit is worse than single fit
        if doub_max_err > sing_max_err:
            print('Double fit error is worse than single; using single fit.')
            params = sing_params
        else: 
            params = doub_params

    return params


def single_arr(temps, kts, t_ref=1.0):
    """ Fit a set of T-dependent rate constants to a single Arrhenius form

        :param temps: temperatures
        :type temps: numpy.ndarray
        :param kts: rate constants
        :type kts: numpy.ndarray
        :param t_ref: reference temperature (K)
        :type t_ref: float
        :return fit_params: A, n, Ea fitting parameters for function
        :rtype: list(float)
    """

    # Consider several cases depending on the number of valid rate constants
    # no k is positive, so return all zeros
    if kts.size == 0:
        a_fit, n_fit, ea_fit = 0.0, 0.0, 0.0

    # if num(k) is 1: set A = k
    elif kts.size == 1:
        a_fit, n_fit, ea_fit = kts[0], 0.0, 0.0

    # if num(k) > 0 is 2,3: fit A and Ea
    elif kts.size in (2, 3):
        # Build vectors and matrices used for the fitting
        a_vec = numpy.ones(len(temps))
        ea_vec = (-1.0 / RC) * (1.0 / temps)
        coeff_mat = numpy.array([a_vec, ea_vec], dtype=numpy.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = numpy.log(kts)
        # Perform the least-squares fit
        theta = numpy.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = numpy.exp(theta[0]), 0.0, theta[1]

    # if num(k) > 0 is more than 3: fit A, n, and Ea
    elif kts.size > 3:
        # Build vectors and matrices used for the fitting
        a_vec = numpy.ones(len(temps))
        n_vec = numpy.log(temps / t_ref)
        ea_vec = (-1.0 / RC) * (1.0 / temps)
        coeff_mat = numpy.array([a_vec, n_vec, ea_vec], dtype=numpy.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = numpy.log(kts)
        # Perform the least-squares fit
        theta = numpy.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = numpy.exp(theta[0]), theta[1], theta[2]

    # Pack the parameters into an arr_dct and instantiate RxnParams
    arr_dct = {'arr_tuples': [[a_fit, n_fit, ea_fit],]}
    params = RxnParams(arr_dct=arr_dct)    

    return params


def double_arr_iter(temps, kts, sing_params, t_ref=1.0, dbl_tol=15, dbl_iter=30):
    """ This is the version that iterates with random guesses
    """

    def _fit_doub_arr(temps, kts, sing_params, t_ref=1.0):
        """ Performs a double Arrhenius fit
        """

        # Unpack the single Arrhenius params
        sing_a, sing_n, sing_ea = sing_params.arr[0]  # get first (& only) entry

        # Randomly generate variations for the three params within given bnds
        a_bnd = random.uniform(GUESS_BNDS[0][0], GUESS_BNDS[0][1])
        n_bnd = random.uniform(GUESS_BNDS[1][0], GUESS_BNDS[1][1])
        ea_bnd = random.uniform(GUESS_BNDS[2][0], GUESS_BNDS[2][1])
        init_guess = [(sing_a * a_bnd), (sing_n + n_bnd), sing_ea / ea_bnd,
                      (sing_a / a_bnd), (-sing_n - n_bnd), sing_ea * ea_bnd]
        # These were the original initial guesses
        #init_guess = [(sing_a * a_bnd), (sing_n + n_bnd), sing_ea / ea_bnd,
                      #(sing_a / a_bnd), (sing_n - n_bnd), sing_ea * ea_bnd]
    
        # Perform a least-squares fit
        plsq = leastsq(_resid_func, init_guess,
                       args=(temps, kts, t_ref),
                       ftol=1.0E-8, xtol=1.0E-8, maxfev=100000)
        
        # Retrieve the fit params and instantiate RxnParams
        raw_params = list(plsq[0])  # a list of length 6
        arr_dct = {'arr_tuples': [raw_params[:3], raw_params[3:]]}
        params = RxnParams(arr_dct=arr_dct)    

        return params


    # Make a maximum of dbl_iter attempts at a double fit
    max_errs = []
    prev_params = []
    ref_ktp_dct = {'high': (temps, kts)}  # used for err_dct later
    for guess_idx in range(dbl_iter):
        # Perform a double fit
        params = _fit_doub_arr(temps, kts, sing_params, t_ref=t_ref)
        
        # Get an err_dct and the max error
        err_dct = new_err.get_err_dct(ref_ktp_dct, params)
        max_err = new_err.get_max_err(err_dct)
        max_errs.append(max_err)
        prev_params.append(params)

        # Exit the loop if tolerance is satisfied
        if max_err < dbl_tol:
            break
        # If tolerance not satisfied and guesses have been exhausted, grab
        # the result with the lowest error
        elif guess_idx == dbl_iter - 1: 
            min_idx = numpy.argmin(max_errs) 
            params = prev_params[min_idx]

    return params, guess_idx


def double_arr_orig(temps, kts, sing_params, t_ref=1.0, dbl_tol=15, dbl_iter=30):
    """ This is the new version
    """
 
    a_vec = numpy.ones(len(temps))
    n_vec = numpy.log(temps / t_ref)
    ea_vec = (-1.0 / RC) * (1.0 / temps)
    coeff_mat = numpy.array([a_vec, n_vec, ea_vec, a_vec, n_vec, ea_vec,],
                            dtype=numpy.float64)
    coeff_mat = coeff_mat.transpose()
    k_vec = numpy.log(kts)
 
    # Perform the least-squares fit
    theta = numpy.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
    # Set the fitting parameters
    a_fit1, n_fit1, ea_fit1, a_fit2, n_fit2, ea_fit2 = (
        numpy.exp(theta[0]), theta[1], theta[2], numpy.exp(theta[3]), 
        theta[4], theta[5]
    )
   
    # Pack the parameters into an arr_dct and instantiate RxnParams
    arr_dct = {'arr_tuples': [[a_fit1, n_fit1, ea_fit1], [a_fit2, n_fit2, ea_fit2],]}
    params = RxnParams(arr_dct=arr_dct)    

    return params


def double_arr_pso(temps, kts, t_ref):
    """ Performs a double Arrhenius fit using particle swarm optimization
    """
    
    def _resid_func_pso(curr_guess):
        """ Computes the residual between fit and data for double fitter
    
            :param curr_guess: current guess for double Arrhenius params
            :type curr_guess: Numpy array of shape (num_particles, num_dims)
            :return resid: SSE residual between fit and data
            :rtype: Numpy array of shape (num_particles,)
        """

        num_particles = numpy.shape(curr_guess)[0]
        resids = numpy.zeros(num_particles)
        for pidx in range(num_particles):
            k_fit1 = numpy.exp(
                numpy.log(curr_guess[pidx, 0]) +
                curr_guess[pidx, 1] * numpy.log(temps / t_ref) -
                curr_guess[pidx, 2] / (RC * temps)
            )
            k_fit2 = numpy.exp(
                numpy.log(curr_guess[pidx, 3]) +
                curr_guess[pidx, 4] * numpy.log(temps / t_ref) -
                curr_guess[pidx, 5] / (RC * temps)
            )
            k_fit = k_fit1 + k_fit2
            resids[pidx] = sum((numpy.log10(kts) - numpy.log10(k_fit)) ** 2)

        return resids

    a_bnd = 1e100
    n_bnd = 100
    ea_bnd = 150e3

    # Run the optimizer
    options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}  # so-called "hyperparameters"
    bounds = (numpy.array([-a_bnd, -n_bnd, -ea_bnd, -a_bnd, -n_bnd, -ea_bnd]),
              numpy.array([a_bnd, n_bnd, ea_bnd, a_bnd, n_bnd, ea_bnd]))
    optimizer = ps.single.GlobalBestPSO(n_particles=2000, dimensions=6,
                                        options=options, bounds=bounds)
    # Output is a tuple: (cost, final_pos), where final_pos is a list of 
    # length num_dimensions
    output = optimizer.optimize(_resid_func_pso, iters=10000)
    
    
    # Retrieve the fit params and instantiate RxnParams
    arr_params1 = list(output[1][:3]) 
    arr_params2 = list(output[1][3:]) 
    arr_dct = {'arr_tuples': [arr_params1, arr_params2]}
    params = RxnParams(arr_dct=arr_dct)    

    return params


def _resid_func(curr_guess, temps, kts, t_ref):
    """ Computes the residual between fit and data for double fitter

        :param curr_guess: current guess for double Arrhenius params
        :type curr_guess: list [A1, n1, Ea1, A2, n2, Ea2]
        :param temps: temperatures
        :type temps: Numpy.ndarray of shape (num_temps,)
        :param kts: rate constants that are being fit
        :type kts: Numpy.ndarray of shape (num_temps,)
        :param t_ref: reference temperature
        :type t_ref: float
        :return resid: residual between fit and data
        :rtype: Numpy.ndarray of shape (num_temps,)
    """

    # Compute the fitted rate constant
    k_fit1 = numpy.exp(
        numpy.log(curr_guess[0]) +
        curr_guess[1] * numpy.log(temps / t_ref) -
        curr_guess[2] / (RC * temps)
    )
    k_fit2 = numpy.exp(
        numpy.log(curr_guess[3]) +
        curr_guess[4] * numpy.log(temps / t_ref) -
        curr_guess[5] / (RC * temps)
    )
    k_fit = k_fit1 + k_fit2
    resid = numpy.log10(kts) - numpy.log10(k_fit)

    return resid
