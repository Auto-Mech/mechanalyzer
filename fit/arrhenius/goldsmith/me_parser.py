#!/usr/bin/env python
"""
  Library of functions to parse the MESS output file used to
  obtain rate constants and fit them to Arrhenius and PLOG expressions
"""

import sys
import os
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot
import matplotlib.pylab
import matplotlib.gridspec

class paper:
    start_here = 0
    color_list = ['BurlyWood', 'Blue', 'Red', 'Fuchsia',
                  'Chartreuse', 'Black', 'Aqua', 'DarkGray', 'Gold']
    Tmin = 0.0
    Tmax = 1.0E4

class reaction:
    theory = 0.0

class fit_type:
    style = 'single'


def get_temp_pres(data, outfile):
    """ This subroutine reads in the me.out file and pulls the
        temperatures, pressures, and species names
    """
    temperature = []
    pressure = []
    channels = []
    wells = []
    bimolecular = []
    temperature_units = []
    pressure_units = []
    for (i, line) in enumerate(outfile):
        if line.startswith('Wells '):
            j = i + 2
            found_well = True
            while found_well:
                well_line = outfile[j].split()
                if len(well_line) > 0:
                    j = j + 1
                    wells.append(well_line[0])
                else:
                    found_well = False
        elif line.startswith('Bimolecular Products '):
            j = i + 2
            found_bimolecular = True
            while found_bimolecular:
                bimolecular_line = outfile[j].split()
                if len(bimolecular_line) > 0:
                    j = j + 1
                    bimolecular.append(bimolecular_line[0])
                else:
                    found_bimolecular = False
        elif line.startswith('Temperature = '):
            bits = line.split()
            if len(bits) == 8: #we have pdep
                temp = float(bits[2])
                pres = float(bits[6])
                if temp not in temperature:
                    temperature.append(temp)
                if pres not in pressure:
                    pressure.append(pres)
                T_unit = bits[3]
                P_unit = bits[7]
                if T_unit not in temperature_units:
                    temperature_units.append(T_unit)
                if P_unit not in pressure_units:
                    pressure_units.append(P_unit)
        elif line.startswith("From\To"):
            bits = line.split()
            for species in bits[1:]:
                if species not in channels:
                    channels.append(species)
        elif line.startswith('Capture/Escape'):
            start_here = i
            break

    data.wells = wells
    data.bimolecular = bimolecular
    data.temperature = temperature
    data.pressure = pressure
    data.temperature_units = temperature_units
    data.pressure_units = pressure_units
    data.channels = channels
    data.start_here = start_here

    # Now create the reaction class for each reaction.
    print(channels)
    sys.exit()
    for reactant in channels:
        for product in channels:
            if reactant != product:
                rxn = reaction()
                print('SPECIES:')
                print(reactant)
                print(product)
                rxn.reactant = reactant
                rxn.product = product
                rxn.title = reactant + '->' + product
                rxn.theory = np.zeros([len(pressure), len(temperature)], dtype=np.float64)
                rxn.s_fit = []
                rxn.d_fit = []

    # create the PLOG class for each pressure
    for p in range(len(pressure)):
        s_plog = fit_type()
        s_plog.A = [0.0, 0.0]
        s_plog.n = [0.0, 0.0]
        s_plog.Ea = [0.0, 0.0]
        rxn.s_fit.append(s_plog)
        d_plog = fit_type()
        d_plog.A = [0.0, 0.0]
        d_plog.n = [0.0, 0.0]
        d_plog.Ea = [0.0, 0.0]
        rxn.d_fit.append(d_plog)
        data.reactions.append(rxn)

    for reactant in channels:
        # add in the capture rate, too
        rxn = reaction()
        rxn.reactant = reactant
        rxn.product = 'Capture'
        rxn.title = reactant + '=' + 'Capture'
        rxn.theory = np.zeros([len(pressure), len(temperature)], dtype=np.float64)
        rxn.s_fit = []
        rxn.d_fit = []
        # create the PLOG class for each pressure
        for p in range(len(pressure)):
            s_plog = fit_type()
            s_plog.A = [0.0, 0.0]
            s_plog.n = [0.0, 0.0]
            s_plog.Ea = [0.0, 0.0]
            rxn.s_fit.append(s_plog)
            d_plog = fit_type()
            d_plog.A = [0.0, 0.0]
            d_plog.n = [0.0, 0.0]
            d_plog.Ea = [0.0, 0.0]
            rxn.d_fit.append(d_plog)
        data.reactions.append(rxn)

    return


def get_pdep_k(data, outfile):
    """ this subrountine continues to read in the output file and
        parses it into the pressure dependent rate constant for all the reactions
    """

    temperature = data.temperature
    pressure = data.pressure
    channels = data.channels

    for (l, line) in enumerate(outfile[data.start_here:]):
        # use the keyword "Pressure" to begin
        if line.startswith('   Pressure ='):
            bits = line.split()
            pres = float(bits[2])
            p = pressure.index(pres)
            # move forward two lines to read of the reactants and products
            rxn_line = outfile[data.start_here + l + 2].split()[1:]
            print(pres, rxn_line)
            if (len(rxn_line) - len(channels) - 1) != 0:
                print(pres, rxn_line)
            # initialize dummy matrix
            local_k = np.zeros([len(temperature), len(channels)+1], dtype=np.float64)
            for j in range(len(temperature)):
                # move forward to copy the results to the dummy matrix
                bits = outfile[data.start_here + l + 3 + j].split()
                t = temperature.index(float(bits[0]))
                # if PAPER cannot compute the rate constant, it prints "***".
                # replace this with -1.0
                for (b, bit) in enumerate(bits):
                    if bit == '***':
                        bits[b] = -1.0
                local_k[t, :] = bits[1:]
            # now cycle through the reactions
            for (r, rx) in enumerate(rxn_line):
                #Capture has no -> symbol and must be treated differently
                if rx != 'Capture':
                    # determine what the reactants and products are
                    reactant, product = rx.split('->')
                    # check to see if the reaction has bimolecular reactants.
                    # the preferred (default) input in chemkin is moles, not molecules
                    if reactant in data.wells:
                        adjust_units = 1.0
                    elif reactant in data.bimolecular:
                        # this parameter is set at the beginning of the locally run parser script.
                        adjust_units = data.N_avo

                    for rxn in data.reactions:
                        if (rxn.reactant == reactant) and (rxn.product == product): # a match!
                            rxn.theory[p, :] = local_k[:, r] * adjust_units
                        # net rate has no product (e.g. "A->"),
                        # and previously we assigned this to be "A->A"
                        elif (rxn.reactant == reactant) and (rxn.product == reactant) and (product == ''): #we have a match!
                            rxn.theory[p, :] = local_k[:, r] * adjust_units
                #now go back through and add the capture rate
                else:
                    for rxn in data.reactions:
                        if (rxn.reactant == reactant) and (rxn.product == 'Capture'): #we have a match!
                            rxn.theory[p, :] = local_k[:, r] * adjust_units

    return


def fit_pdep(data, nonlin_fit=True, just_reactant=['All'], skip_product=['All']):
    """ the subrountine cycles through all the reactions and
        obtains the Arrhenius fits for the PLOG expressions.
        when appropriate, it will also do the sum of two Arrhenius expressions.
    """

    # cycle through the reactions
    for rxn in data.reactions:
        print("fitting reaction: ", rxn.title)
        # first get single fits for each pressure
        for (p, s_plog) in enumerate(rxn.s_fit):
            k = rxn.theory[p, :]
            get_single_arrhenius(data, k, s_plog)
        # now cycle back through for the double plog for each pressure
        if nonlin_fit:
            if just_reactant[0] == 'All':
                for (p, d_plog) in enumerate(rxn.d_fit):
                    k = rxn.theory[p, :]
                    # if there are fewer than seven rate constants that are positive,
                    # then there is no point to using double plog
                    N_good_data = 0
                    for kk in k:
                        if kk > 0.0:
                            N_good_data = N_good_data + 1
                    # we have at least seven good k's, so call double plog.
                    # The nonlin_fit specifies whether or not to use the nonlinear solver
                    if N_good_data > 6:
                        get_double_arrhenius(data, k, d_plog, nonlin_fit)
                    # We have six or fewer good k's, so copy over the singles data
                    else:
                        local_s_plog = rxn.s_fit[p]
                        d_plog.A = local_s_plog.A
                        d_plog.n = local_s_plog.n
                        d_plog.Ea = local_s_plog.Ea
                        d_plog.fit_range = local_s_plog.fit_range
                        d_plog.MAE = local_s_plog.MAE
            else:
                if (rxn.reactant in just_reactant and rxn.product not in skip_product):
                    for (p, d_plog) in enumerate(rxn.d_fit):
                        k = rxn.theory[p, :]
                        # if there are fewer than seven rate constants that are positive,
                        # then there is no point to using double plog
                        N_good_data = 0
                        for kk in k:
                            if kk > 0.0:
                                N_good_data = N_good_data + 1
                        # we have at least seven good k's, so call double plog.
                        # The nonlin_fit specifies whether or not to use the nonlinear solver
                        if N_good_data > 6:
                            get_double_arrhenius(data, k, d_plog, nonlin_fit)
                        # We have six or fewer good k's, so copy over the singles data
                        else:
                            local_s_plog = rxn.s_fit[p]
                            d_plog.A = local_s_plog.A
                            d_plog.n = local_s_plog.n
                            d_plog.Ea = local_s_plog.Ea
                            d_plog.fit_range = local_s_plog.fit_range
                            d_plog.MAE = local_s_plog.MAE
                else:
                    print("skipping double fit for ", rxn.title)
                    for (p, d_plog) in enumerate(rxn.d_fit):
                        d_plog.fit_range = [0, 0]
                        d_plog.MAE = [0, 0]
        else:
            for (p, d_plog) in enumerate(rxn.d_fit):
                d_plog.fit_range = [0, 0]
                d_plog.MAE = [0, 0]

    return


def get_valid_Tk(data, k):
    """ this subroutine takes in a array of rate constants and
        returns the subset of this array that is positive,
        along with the corresponding Temperature array """

    # start by using only the temperatures at which the rate constant is well defined
    local_T = []
    local_k = []
    for (t, T) in enumerate(data.temperature):
        if (k[t] > 0.0) and (T >= data.Tmin) and (T <= data.Tmax):
            local_T.append(T)
            local_k.append(k[t])

    local_T = np.array(local_T, dtype=np.float64)
    local_k = np.array(local_k, dtype=np.float64)

    return local_T, local_k


def fit_arrhenius(data, k):
    """ this subroutine takes in a vector of rate constants and
        returns the Arrhenius parameters, as well as
        the T-range over which they were fit"""

    # start by using only the temperatures at which the rate constant is well defined
    local_T, local_k = get_valid_Tk(data, k)

    # there are many cases to consider, depending upon the number of valid k's
    # no k is positive, so return all zeros
    if len(local_k) == 0:
        A, n, Ea = 0.0, 0.0, 0.0
        fit_range = [0, 0]
    
    # one k is positive, so set A = k
    elif len(local_k) == 1:
        A, n, Ea = local_k, 0.0, 0.0
        fit_range = [data.temperature.index(min(local_T)), data.temperature.index(max(local_T))]
    
    # 2 or 3 k's are positive, so fit A and Ea
    elif (len(local_k) == 2) or (len(local_k) == 3):
        X = np.array([np.ones(len(local_T)), -1.0 / data.R / local_T], dtype=np.float64)
        X = X.transpose()
        theta = np.linalg.lstsq(X, np.log(local_k))[0]
        A = np.exp(theta[0])
        n = 0.0
        Ea = theta[1]
        fit_range = [data.temperature.index(min(local_T)), data.temperature.index(max(local_T))]

    # more than 3 k's are positive, so fit A, n, Ea
    elif len(local_k) > 3:
        X = np.array([np.ones(len(local_T)), np.log(local_T/data.T0), -1.0 / data.R / local_T],
                     dtype=np.float64)
        X = X.transpose()
        theta = np.linalg.lstsq(X, np.log(local_k))[0]
        A = np.exp(theta[0])
        n = theta[1]
        Ea = theta[2]
        fit_range = [data.temperature.index(min(local_T)), data.temperature.index(max(local_T))]

    return A, n, Ea, fit_range


def get_SSE(data, k, plog):
    """ (1) get the sum of square error (SSE) useful when determining
            which double plog routine will be used to initialize the nonlinear solver
        (2) also get the mean absolute error (MAE), which is written to the plog file
    """

    # start by using only the temperatures at which the rate constant is well defined
    local_T, local_k = get_valid_Tk(data, k)
    if len(local_k) > 2:
        # compute the fitted rate constant for the valid temperature range
        A = plog.A
        n = plog.n
        Ea = plog.Ea
        k_fit = A[0] * (local_T/data.T0)**n[0] * np.exp(-Ea[0]/data.R/local_T)
        k_fit = k_fit + A[1] * (local_T/data.T0)**n[1] * np.exp(-Ea[1]/data.R/local_T)
        # initialize SSE and MAE
        SSE = 0.0
        MAE = []
        for t in range(len(local_k)):
            SSE = SSE + (np.log(local_k[t]) - np.log(k_fit[t]))**2.0
            MAE.append(np.abs((local_k[t] - k_fit[t]) / local_k[t]))
        MAE = np.array(MAE, dtype=np.float64)
        MAE = [np.mean(MAE)*100.0, np.max(MAE)*100.0]
    else:
        SSE = 0.0
        MAE = [0.0, 0.0]
    plog.SSE = SSE
    plog.MAE = MAE

    return


def get_single_arrhenius(data, k, s_plog):
    """ Just a simple wrapper to call the Arrhenius subroutine and
        then get the SSE and MAE from the fits
    """

    s_plog.A[0], s_plog.n[0], s_plog.Ea[0], s_plog.fit_range = fit_arrhenius(data, k)
    get_SSE(data, k, s_plog)

    return


def get_double_arrhenius(data, k, d_plog, nonlin_fit):
    """ the main subroutine for obtaining the sum of two Arrhenius expressions.
        (1) Basically, the subroutine first determines whether the curve is
            bending up or down at high temperature.
        (2) Based upon the curvature, one of two possible cycles are initiated
            to provide a good starting guess for the nonlinear solver."""

    # start by using only the temperatures at which the rate constant is well defined
    local_T, local_k = get_valid_Tk(data, k)
    T_mid = 1.0/np.mean(1.0/local_T)
    fit_range = [data.temperature.index(min(local_T)), data.temperature.index(max(local_T))]

    s_A, s_n, s_Ea, s_fit_range = fit_arrhenius(data, k)

    first_guess = [s_A/2.0, s_n, s_Ea, s_A/2.0, s_n, s_Ea]
    # compute the new fitted rate constant
    local_fit = s_A * (local_T/data.T0)**s_n * np.exp(-s_Ea/data.R/local_T)
    # check to make sure there are no errors
    if min(local_fit) < 0.0:
        print('we have a problem')
    # compute the local sum of square errors
    first_SSE = 0.0
    for t in range(len(local_T)):
        first_SSE = first_SSE + (np.log(local_k[t])-np.log(local_fit[t]))**2.0

    sjk_input = open('arrfit.dat', 'w')
    header = ' 1.0e-10 1.0e-10 \n 1 1 1 1 1 1 \n 1 \n   8.1e-11, -0.01,  2000. \n'
    sjk_input.write(header)
    newline = "%3d  %3d \n"%(len(local_T), len(local_T))
    sjk_input.write(newline)
    for i in range(len(local_T)):
        newline = "%2.1F  %2.4E \n"%(local_T[i], local_k[i])
        sjk_input.write(newline)
    sjk_input.write('\n')
    sjk_input.close()
    command = './dsarrfit.x_cfg'
    os.system(command)

    sjk_input = open('arrfit.out', 'r')
    lines = sjk_input.readlines()
    lines.reverse()
    for line in lines:
        if line.startswith(' params'):
            bits1 = lines[lines.index(line)-1].split()
            bits2 = lines[lines.index(line)-2].split()
            break
    if (len(bits1) == 6) and (len(bits2) == 6):
        job = True
        best_guess = np.zeros(len(bits1), dtype=np.float64)
        for (b, bit) in enumerate(bits1):
            if '+' in bit:
                if bit[bit.index('+') - 1] != 'E':
                    bit = bit.replace('+', 'E+')
                # zador
                if len(bit) >= 9:
                    if (bit[8] == '-') and (bit[7] != 'E'):
                        bit = bit.replace('-', 'X', 1)
                        bit = bit.replace('-', 'E-')
                        bit = bit.replace('X', '-')
                if len(bit) >= 8:
                    if (bit[7] == '-') and (bit[6] != 'E'):
                        bit = bit.replace(bit[7], 'E-')

            best_guess[b] = float(bit)

        best_guess[2] = best_guess[2]*1.987
        best_guess[5] = best_guess[5]*1.987
    else:
        print('FAILED! Falling back to python solver')
        job = False

    if not job == False:
        A1 = s_A / 2.0
        n1 = s_n + 0.1
        # A1 = A1 / (T_mid/data.T0)**n1
        Ea1 = s_Ea
        A2 = s_A / 2.0
        n2 = s_n - 0.1
        # A2 = A2 / (T_mid/data.T0)**n2
        Ea2 = s_Ea

        new_guess = [A1, n1, Ea1, A2, n2, Ea2]
        plsq = leastsq(mod_arr_residuals, new_guess, args=(local_k, local_T, data),
                       ftol=1.0E-9, xtol=1.0E-9, maxfev=100000)

        # Dfun=mod_arr_jacobian_finite_difference,col_deriv=1, # epsfcn=1.0E-6 # factor=0.1,
        [A1, n1, Ea1, A2, n2, Ea2] = plsq[0]

        # compute the new fitted rate constant
        local_fit = A1 * (local_T/data.T0)**n1 * np.exp(-Ea1/data.R/local_T)
        local_fit = local_fit + A2 * (local_T/data.T0)**n2 * np.exp(-Ea2/data.R/local_T)

        # check to make sure there are no errors
        if min(local_fit) < 0.0:
            print('we have a problem')

        # compute the local sum of square errors
        local_SSE = 0.0
        for t in range(len(local_T)):
            local_SSE = local_SSE + (np.log(local_k[t])-np.log(local_fit[t]))**2.0
        # check to see if this run is better than all previous runs
        if local_SSE < first_SSE:
            best_guess = [A1, n1, Ea1, A2, n2, Ea2]
        else:
            print(" nonlinear solver was not any better than single exponential!")
            best_guess = first_guess

    # assign the results
    d_plog.fit_range = fit_range                # fit range
    d_plog.A = [best_guess[0], best_guess[3]]   # A factors
    d_plog.n = [best_guess[1], best_guess[4]]   # n
    d_plog.Ea = [best_guess[2], best_guess[5]]  # Ea

    # finally, get the SSE and MAE for the total fit
    get_SSE(data, k, d_plog)

    print("  Mean A.E %1.2F%%, Max A.E. %1.2F%%,     double_SSE / single_SSE %1.2E"%(d_plog.MAE[0], d_plog.MAE[1], d_plog.SSE/first_SSE)) # best_guess

    return


def get_double_arrhenius_OLD(data, k, d_plog, nonlin_fit):
    """ the main subroutine for obtaining the sum of two Arrhenius expressions.
        Basically, the subroutine first determines whether the curve is bending up or down at high temperature.
        Based upon the curvature, one of two possible cycles are initiated to provide a good starting guess for the nonlinear solver."""
    # start by using only the temperatures at which the rate constant is well defined
    local_T, local_k = get_valid_Tk(data, k)
    fit_range = [data.temperature.index(min(local_T)), data.temperature.index(max(local_T))]

    s_A, s_n, s_Ea, s_fit_range = fit_arrhenius(data, k)


    # call the subroutine to determine the curvature of the rate constant
    d_plog.curvature = get_k_curvature(data, local_T, local_k)
    # start with curves that are concave down, i.e. d^2k/dT^2 <=0
    if d_plog.curvature == 'down':
        factor = 1.0/100.0
        # start by adjusting the initial theory value
        k_new = np.zeros(len(k))

        for i in range(len(k_new)):
            k_new[i] = k[i] * 0.9
            k_new[i] = k_new[i] * (1.0 -  (1.0-factor) * float(i) / float(len(k_new)-1.0))

        A1, n1, Ea1, fit_range1 = fit_arrhenius(data, k_new)
        k_diff = []
        for (t, T) in enumerate(data.temperature):
            k_fit = A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T)
            k_diff.append(k[t]-k_fit)

        k_diff = np.array(k_diff, dtype=np.float64)
        A2, n2, Ea2, fit_range2 = fit_arrhenius(data, k_diff)
        d_plog.fit_range = [max(fit_range1[0], fit_range2[0]), min(fit_range1[1], fit_range2[1])]
        p0 = [A1, n1, Ea1, A2, n2, Ea2]
        if nonlin_fit:
            plsq = leastsq(mod_arr_residuals, best_guess, args=(local_k, local_T, data),
                           Dfun=mod_arr_jacobian_finite_difference, col_deriv=1,
                           factor=0.1, ftol=1.0E-8, xtol=1.0E-8, maxfev=10000)
        else:
            plsq = [p0, 0]

    # now switch to curves that are convex up, i.e. d^2k/dT^2 >=0
    elif d_plog.curvature == 'up':
        N_runs = len(local_T) - 6
        best_guess = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        best_SSE = 1e6
        for run in range(N_runs):
            new_T = local_T[3+run:len(local_T)-1]
            new_k = local_k[3+run:len(local_T)-1]
            X = np.array([np.ones(len(new_T)), np.log(new_T/data.T0), -1.0 / data.R / new_T],
                         dtype=np.float64)
            X = X.transpose()
            theta = np.linalg.lstsq(X, np.log(new_k))[0]
            A1 = np.exp(theta[0])
            n1 = theta[1]
            Ea1 = theta[2]
            local_T_2 = []
            local_k_2 = []
            # now get the difference between the theory and the fitted.
            diff = []
            for (t, T) in enumerate(data.temperature):
                diff.append(k[t] - A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T))
            # fit the difference
            A2, n2, Ea2, fit_range2 = fit_arrhenius(data, diff)
            # compute the new fitted rate constant
            local_fit = A1 * (local_T/data.T0)**n1 * np.exp(-Ea1/data.R/local_T)
            local_fit = local_fit + A2 * (local_T/data.T0)**n2 * np.exp(-Ea2/data.R/local_T)
            # check to make sure there are no errors
            if min(local_fit) < 0.0:
                print('we have a problem')
            # compute the local sum of square errors
            local_SSE = 0.0
            for t in range(len(local_T)):
                local_SSE = local_SSE + (np.log(local_k[t])-np.log(local_fit[t]))**2.0
            # check to see if this run is better than all previous runs
            if local_SSE < best_SSE:
                best_SSE = local_SSE
                best_guess = [A1, n1, Ea1, A2, n2, Ea2]
        # p0 is the initial guess for the nonlinear solver
        p0 = best_guess
        # check to see if the nonlinear solver should be called
        if nonlin_fit:
            plsq = leastsq(mod_arr_residuals, best_guess, args=(local_k, local_T, data),
                           Dfun=mod_arr_jacobian_finite_difference, col_deriv=1,
                           ftol=1.0E-12, xtol=1.0E-12, maxfev=10000)
        else:
            plsq = [p0, 0]
    # assign the results
    d_plog.fit_range = fit_range            # fit range
    d_plog.A = [plsq[0][0], plsq[0][3]]     # A factors
    d_plog.n = [plsq[0][1], plsq[0][4]]     # n
    d_plog.Ea = [plsq[0][2], plsq[0][5]]    # Ea
    # finally, get the SSE and MAE for the total fit
    get_SSE(data, k, d_plog)
    print("%1.2F  %1.2F"%(d_plog.MAE[0], d_plog.MAE[1]))

    return


def get_k_curvature(data, T, k):
    """ this subroutine determines whether
        the rate constant curves up or down.
    """

    # just use the first and last value
    local_T = np.array([T[0], T[len(T)-1]], dtype=np.float64)
    local_k = np.array([k[0], k[len(T)-1]], dtype=np.float64)

    # solve for a standard (non-modified) Arrhenius expression
    X = np.array([np.ones(len(local_T)), -1.0 / data.R / local_T], dtype=np.float64)
    X = X.transpose()
    theta = np.linalg.lstsq(X, np.log(local_k))[0]
    A = np.exp(theta[0])
    Ea = theta[1]
    fit_k = A * np.exp(-Ea/data.R/T)

    # now check to see if the fitted rate constant is consistently above or below the theory
    difference = 0.0
    for t in range(len(T)):
        difference = difference + np.log(fit_k[t]) - np.log(k[t])
    # if the fit is consistently higher, then the curve is convex up
    if difference > 0.0:
        curvature = 'up'
    # otherwise if the fit is consistently below, then the curve is concave down
    else:
        curvature = 'down'
    return curvature


def mod_arr_residuals(p, target, T, data):
    " this subroutine computes the residual that goes into the nonlinear solver"

    # compute the fitted rate constant
    k_fit1_log = np.log(p[0]) + p[1] * np.log((T/data.T0)) - p[2]/data.R/T
    k_fit2_log = np.log(p[3]) + p[4] * np.log((T/data.T0)) - p[5]/data.R/T
    k_fit1 = np.exp(k_fit1_log)
    k_fit2 = np.exp(k_fit2_log)
    k_fit = k_fit1 + k_fit2

    log_yes = True
    if log_yes:
        err = np.log10(target) - np.log10(k_fit)
    elif log_yes == False:
        err = (target - k_fit)/target

    return err


def mod_arr_jacobian_finite_difference(p, target, T, data):
    " this subroutine computes the jacobian that goes into the nonlinear solver"

    # read in the parameters
    A1, n1, Ea1, A2, n2, Ea2 = p
    p0 = p
    # delta
    delta = 1.0E2

    param_step = np.zeros(len(p), dtype=np.float64)
    jac = np.zeros([len(p), len(T)], dtype=np.float64)
    for i in range(len(p)):
        param_step[i] = delta * p[i]
        p[i] = p0[i] + param_step[i]
        fit_plus1_log = np.log(p[0]) + p[1] * np.log((T/data.T0)) - p[2]/data.R/T
        fit_plus2_log = np.log(p[3]) + p[4] * np.log((T/data.T0)) - p[5]/data.R/T

        p[i] = p0[i] - param_step[i]
        fit_minus1_log = np.log(p[0]) + p[1] * np.log((T/data.T0)) - p[2]/data.R/T
        fit_minus2_log = np.log(p[3]) + p[4] * np.log((T/data.T0)) - p[5]/data.R/T

        fit_plus1 = np.exp(fit_plus1_log)
        fit_plus2 = np.exp(fit_plus2_log)
        fit_minus1 = np.exp(fit_minus1_log)
        fit_minus2 = np.exp(fit_minus2_log)

        fit_plus = np.log10(fit_plus1 + fit_plus2)
        fit_minus = np.log10(fit_minus1 + fit_minus2)

        jac[i] = (fit_plus - fit_minus) / (2.0 * param_step[i])
        p = p0
    '''
    # compute the jacobian
    jac_A1  = ( np.log( (A1*(1+delta)) * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T) + A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T)) \
            -   np.log( (A1*(1-delta)) * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T) + A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T))) / (2*delta*A1)

    jac_n1  = ( np.log( A1 * (T/data.T0)**(n1*(1+delta)) * np.exp(-Ea1/data.R/T) + A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T)) \
            -   np.log( A1 * (T/data.T0)**(n1*(1-delta)) * np.exp(-Ea1/data.R/T) + A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T))) / (2*delta*n1)

    jac_Ea1 = ( np.log( A1 * (T/data.T0)**n1 * np.exp(-(Ea1*(1+delta))/data.R/T) + A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T)) \
            -   np.log( A1 * (T/data.T0)**n1 * np.exp(-(Ea1*(1-delta))/data.R/T) + A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T))) / (2*delta*Ea1)

    jac_A2  = ( np.log( A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T) + (A2*(1+delta)) * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T)) \
            -   np.log( A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T) + (A2*(1-delta)) * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T))) / (2*delta*A2)

    jac_n2  = ( np.log( A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T) + A2 * (T/data.T0)**(n2*(1+delta)) * np.exp(-Ea2/data.R/T)) \
            -   np.log( A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T) + A2 * (T/data.T0)**(n2*(1-delta)) * np.exp(-Ea2/data.R/T))) / (2*delta*n2)

    jac_Ea2 = ( np.log( A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T) + A2 * (T/data.T0)**n2 * np.exp(-(Ea2*(1+delta))/data.R/T)) \
            -   np.log( A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T) + A2 * (T/data.T0)**n2 * np.exp(-(Ea2*(1-delta))/data.R/T))) / (2*delta*Ea2)
    jac = [jac_A1, jac_n1, jac_Ea1, jac_A2, jac_n2, jac_Ea2]
    '''
    return jac


def simple_arr_residuals(p, target, T, data, n1, n2):
    """ this subroutine computes the residual that goes into the nonlinear solver
    """
    # read in the parameters
    A1, Ea1, A2, Ea2 = p

    # compute the fitted rate constant
    k_fit1 = A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T)
    k_fit2 = A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T)
    k_fit = k_fit1 + k_fit2

    err = np.zeros(len(T), dtype=np.float64)
    for t in range(len(T)):
        # you have a choice of minimizing the difference in the log of k,
        # or minimizing the relative difference.
        # I prefer the log of k.
        log_yes = True
        if log_yes:
            y0 = np.log(target[t])
            y1 = np.log(k_fit[t])
            err[t] = y0 - y1
        elif log_yes == False:
            y0 = target[t]
            y1 = k_fit[t]
            err[t] = ((y0 - y1)/y0)

    return err


def simplex_residuals(p, target, T, data):
    """ this subroutine computes the residual that goes into the nonlinear solver
    """

    # read in the parameters
    A1, n1, Ea1, A2, n2, Ea2 = p

    # compute the fitted rate constant
    k_fit1 = A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T)
    k_fit2 = A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T)
    k_fit = k_fit1 + k_fit2

    err = np.zeros(len(T), dtype=np.float64)
    for t in range(len(T)):
        # you have a choice of minimizing the difference in the log of k,
        # or minimizing the relative difference.
        # I prefer the log of k.
        log_yes = True
        if log_yes:
            y0 = np.log(target[t])
            y1 = np.log(k_fit[t])
            err[t] = y0 - y1
        elif log_yes == False:
            y0 = target[t]
            y1 = k_fit[t]
            err[t] = ((y0 - y1)/y0)

    return sum(err*err)


def nonlin_wrapper_double(data):
    """ non-linear wrapper function
    """

    def nonlin_double_arrhenius(T, *p):
        """ double Arrhenius
        """

        # read in the parameters
        A1, n1, Ea1, A2, n2, Ea2 = p

        # compute the fitted rate constant
        k_fit1 = A1 * (T/data.T0)**n1 * np.exp(-Ea1/data.R/T)
        k_fit2 = A2 * (T/data.T0)**n2 * np.exp(-Ea2/data.R/T)
        k_fit = k_fit1 + k_fit2

        return np.log(k_fit)

    return nonlin_double_arrhenius


def print_plog(data, me_dot_out):
    """ this subroutine reads in the results from all the fits and
        writes them to a separate file
    """

    pressure = data.pressure
    temperature = data.temperature
    channels = data.channels

    # create new outpu file with same name is the original input prefix, but with .plog appended
    prefix, _ = me_dot_out.split('.')
    plog_file = prefix + '.plog'
    output = open(plog_file, 'w')

    # write a basic header
    header = '!------------------------------------------------------------------------------!\n'
    output.write(header)
    output.write('!\n! RRKM/ME rate constants obtained from PAPER.\n')
    output.write('!\n! MAE = mean/max absolute error, defined as:\n! mean( abs( (theory - fit)/theory))\n!\n')
    output.write(header)

    # start by writing all the single PLOG results.
    output.write('\n\n' + header)
    output.write('!\n! RRKM/ME rate constants fit to single Arrhenius expression\n!\n')
    output.write(header)

    # cycle through each reaction
    for rxn in data.reactions:
        # this corresponds to the net rate constant, so we will replace A=A with A=Net
        if rxn.product == rxn.reactant:
            rxn.title = rxn.reactant + '=Net'

        #start by pulling off the parameters for highest pressure to use as the dummy values
        s_plog = rxn.s_fit[len(pressure)-1]
        rxn_line = "\n" + rxn.title.ljust(40) +"%5.2E    %8.2F    %2.1F\n"%(s_plog.A[0], s_plog.n[0], float("%.4G"%(s_plog.Ea[0])))
        output.write(rxn_line)

        # we don't need pressure dependence for the capture rate constant, so omit it from this stage
        if rxn.product != 'Capture':
            # now cycle through each pressure
            for (p, P) in enumerate(pressure):
                #start by finding the correct pressure, since PLOG must be in atm
                if data.pressure_units[0] == 'atm':
                    P = P
                elif data.pressure_units[0] == 'torr':
                    P = P / 760.0
                elif data.pressure_units[0] == 'bar':
                    P = P / 1.01325
                else:
                    print("valid pressure not recognized")
                    break
                # now begin with the data formatting
                s_plog = rxn.s_fit[p]
                d_plog = rxn.d_fit[p]
                plog_line = "  PLOG/%8.3E    %5.2E    %8.2F    %2.1F/"%(P,s_plog.A[0], s_plog.n[0], float("%.4G"%(s_plog.Ea[0])))
                error_line = "! fit btw. %.F and %.F K with MAE of %.1F%%, %.1F%%\n"%(temperature[s_plog.fit_range[0]], temperature[s_plog.fit_range[1]], s_plog.MAE[0], s_plog.MAE[1])
                local_line = plog_line.ljust(80) + error_line
                output.write(local_line)

    # now repeat everything for the double plog fit
    output.write('\n\n' + header)
    output.write('!\n! RRKM/ME rate constants fit to sum of two Arrhenius expressions\n!\n')
    output.write(header)
    for rxn in data.reactions:
        # this corresponds to the net rate constant, so we will replace A=A with A=Net
        if (rxn.product == rxn.reactant):
            rxn.title = rxn.reactant + '=Net'
        if rxn.product != 'Capture':
            #start by pulling off the parameters for highest pressure to use as the dummy values
            s_plog = rxn.s_fit[len(pressure)-1] #we use the s_plog values for this
            rxn_line = "\n" + rxn.title.ljust(40) +"%5.2E    %8.2F    %2.1F\n"%(s_plog.A[0], s_plog.n[0], float("%.4G"%(s_plog.Ea[0])))
            output.write(rxn_line)
            # we don't need pressure dependence for the capture rate constant, so omit it from this stage
            # now cycle through each pressure
            for (p, P) in enumerate(pressure):
                #start by finding the correct pressure, since PLOG must be in atm
                if data.pressure_units[0] == 'atm':
                    P = P
                elif data.pressure_units[0] == 'torr':
                    P = P / 760.0
                elif data.pressure_units[0] == 'bar':
                    P = P / 1.01325
                else:
                    print("valid pressure not recognized")
                    break
                # now begin with the data formatting
                d_plog = rxn.d_fit[p]
                plog_line_1 = "  PLOG/%8.3E    %5.2E    %8.2F    %.1F/\n"%(P,d_plog.A[0], d_plog.n[0], float("%.4G"%(d_plog.Ea[0])))
                plog_line_2 = "  PLOG/%8.3E    %5.2E    %8.2F    %.1F/"%(P,d_plog.A[1], d_plog.n[1], float("%.4G"%(d_plog.Ea[1])))
                error_line = "! fit btw. %1.F and %1.F K with MAE of %1.1F%%, %1.1F%%\n"%(temperature[d_plog.fit_range[0]], temperature[d_plog.fit_range[1]], d_plog.MAE[0], d_plog.MAE[1])
                local_line = plog_line_1 + plog_line_2.ljust(80) + error_line
                output.write(local_line)
        elif rxn.product == 'Capture':
            d_plog = rxn.d_fit[0] #capture rate constants are not pressure dependent, so just pull the first one
            rxn_line_1 = "\n" + rxn.title.ljust(40) +"%5.2E    %8.2F    %2.1F\n  DUPLICATE"%(d_plog.A[0], d_plog.n[0], float("%.4G"%(d_plog.Ea[0])))
            rxn_line_2 = "\n" + rxn.title.ljust(40) +"%5.2E    %8.2F    %2.1F\n  DUPLICATE\n"%(d_plog.A[1], d_plog.n[1], float("%.4G"%(d_plog.Ea[1])))
            local_line = rxn_line_1 + rxn_line_2
            output.write(local_line)
    output.write('\n\n' + header + '\n\n')
    output.close()

    return

def plot_reactant(data, me_dot_out, show_plot=True, save_plot=True):
    """ this subroutine generates a sequencies of plots.
        There is one figure per species.
        In each figure, there is one pane for each pressure (in order listed in the original input file)
        The lower panes correspond to the errors."""

    temperature = data.temperature
    pressure = data.pressure
    prefix, _ = me_dot_out.split('.')

    # cycle through each species and create a new figure for that species
    for species in data.channels:
        fig = matplotlib.pylab.figure(dpi=120)
        Nplts = len(pressure)
        wdth = np.ones(Nplts)
        hght = [4, 1] # makes the Arrhenius plots four times larger than the relative error plots
        gs = matplotlib.gridspec.GridSpec(2, Nplts, width_ratios=wdth, height_ratios=hght)

        # color-coding for the products
        marker_color = data.color_list

        #first cycle through the pressures to obtain the axis limits.
        #this is purely for visual formatting purposes
        xmin = 1000/min(max(temperature), data.Tmax)
        xmax = 1000/max(min(temperature), data.Tmin)
        ymin = np.inf
        ymax = -np.inf
        for (p, P) in enumerate(pressure):
            for rxn in data.reactions:
                if (rxn.product != 'Capture') and (rxn.reactant != rxn.product) and (rxn.reactant == species):
                    k = rxn.theory[p, :]
                    local_T, local_k = get_valid_Tk(data, k)
                    if max(k) < 0.0:
                        local_min = 0.1
                        local_max = 1.0E-100
                    else:
                        local_min = min(x for x in local_k if x > 0.0)
                        local_max = max(x for x in local_k if x > 0.0)
                    if local_min < ymin:
                        ymin = local_min
                    if local_max > ymax:
                        ymax = local_max

        # finished with the axis formatting.
        # now begin the real process.
        # cycle through each pressure and create a new subplot for that pressure
        for (p, P) in enumerate(pressure):
            ax = matplotlib.pyplot.subplot(gs[p])
            res = matplotlib.pyplot.subplot(gs[Nplts+p])
            # cycle through each reaction at that pressure
            for rxn in data.reactions:
                # check to make sure that only reactants that correspond to this figure are considered.
                # also, exclude the capture rate and the sum of all channels
                if (rxn.reactant == species) and (rxn.product != 'Capture') and (rxn.reactant != rxn.product):
                    # is an index that corresponds to each product. It is used for color-coding purposes only.
                    r = data.channels.index(rxn.product)
                    # right now I only have nine colors coded,
                    # if there are more than nine channels, cycle back through the colors.
                    while r >= len(marker_color):
                        r = r - len(marker_color)
                    # read in the rate constant for this pressure
                    k = rxn.theory[p, :]
                    # start by using only the temperatures at which the rate constant is well defined
                    local_T, local_k = get_valid_Tk(data, k)
                    # plot the raw theory numbers as open circles
                    ax.semilogy(1000.0/local_T, local_k,
                                marker='o', markeredgecolor=marker_color[r],
                                color='w',linestyle='None')

                    # single plog: plot the single plog fit as a dashed line
                    s_plog = rxn.s_fit[p]
                    s_fit = s_plog.A[0] * (local_T/data.T0)**s_plog.n[0] * np.exp(-s_plog.Ea[0]/data.R/local_T)

                    # double plog: plot the double plog fit as a solid line
                    d_plog = rxn.d_fit[p]
                    d_fit1 = d_plog.A[0] * (local_T/data.T0)**d_plog.n[0] * np.exp(-d_plog.Ea[0]/data.R/local_T)
                    d_fit2 = d_plog.A[1] * (local_T/data.T0)**d_plog.n[1] * np.exp(-d_plog.Ea[1]/data.R/local_T)
                    d_fit = d_fit1 + d_fit2
                    ax.semilogy(1000.0/local_T, d_fit, color=marker_color[r],label=rxn.title)

                    # do some axis formatting
                    y_axis_stretch = 6 # give some white space at the top and bottom of each plot
                    ax.set_xlim([xmin/1.25, xmax*1.05])
                    ax.set_ylim([ymin/y_axis_stretch, ymax*y_axis_stretch])
                    ax.tick_params(axis='x', which='major', labelsize=6) # decrease tick size for inverse temperature
                    ax.tick_params(axis='y', which='major', labelsize=8) # decrease tick size for rate constants

                    # now start with the subplots.
                    # for the subplot we plot the ratio of the theory/fit or fit/theory, whichever is greater.
                    # thus, the lowest value will always be zero, which makes it easier to see what the maximum error is
                    s_rat = np.zeros(len(local_k)) # initialize single plog ratio
                    d_rat = np.zeros(len(local_k)) # initialize double plog ratio
                    for t in range(len(local_k)):
                        s_rat[t] = max( local_k[t]/s_fit[t], s_fit[t]/local_k[t]) #single
                        d_rat[t] = max( local_k[t]/d_fit[t], d_fit[t]/local_k[t]) #double
                    res.plot(1000.0/local_T, d_rat, color=marker_color[r]) #double = solid
                    res.set_xlim([xmin/1.25, xmax*1.05])
                    res.tick_params(axis='both', which='major', labelsize=6)

        ax.legend(loc='lower right') # let the code determine where best to put the legend
        fig.set_size_inches(10, 8)
        gs.tight_layout(fig) # tight layout

        # if we want to save the figures, we could do that here
        if save_plot:
            plot_file = prefix + '.' + species.strip() + '.pdf'
            fig.savefig(plot_file)
    if show_plot:
        matplotlib.pyplot.show()

    return
