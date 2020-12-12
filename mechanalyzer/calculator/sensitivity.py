""" Calculates and plots A-factor sensitivities using Cantera
"""

import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib import cm
from PyPDF2 import PdfFileMerger
import os
import multiprocessing
import math

def multiple_sens_arrays(mech_filename, target_spc, T_K, P_atm, X, time, reactor_type='IG const TP', num_top_rxns=10, target_rxns='all', parallel=True):
    """ Runs multiple A-factor sensitivity analyses in Cantera. Returns a list of
        lists of reaction indices and a dictionary of sensitivity array, where each
        value is a different 2-D sensitivity array for a specific simulation.

        Note: all inputs must either a single value or a list of length m. The only exception is the
        time input, which should be a single 1-D array of length n (n does not need to equal m).

        :param mech_filename: name of the .cti mechanism
        :type mech_filename: str
        :param target_spc: target species for sensitivity
        :type target_spc: str
        :param T_K: initial temperature, in Kelvin
        :type T_K: float
        :param P_atm: initial pressure, in atmospheres
        :type P_atm: float
        :param X: initial species and their mole fractions
        :type X: str, example: 'N2O: 0.002, O2: 0.02, AR: 0.996'
        :param time: array of times for evaluating solution, in seconds
        :type time: 1-D Numpy array
        :param reactor_type: type of Cantera reactor
        :type reactor_type: str
        :param num_top_rxns: the number of most-sensitive reactions to be output
        :type num_top_rxns: float
        :param target_rxns: indices of reactions to consider, or 'all'
        :type target_rxns: list or str
        :return: all_rxn_idxs: the indices of the sensitive reactions
        :rtype: list[list]
        :return: all_sens_results: array of the sensitivities for each rxn vs. time
        :rtype: list[2-D Numpy array]
    """
    # Get the max length of the inputs; this is also the number of simulations
    inputs = locals()  # produces a dictionary of the input values
    time = inputs['time']  # store the time for later
    inputs.pop('time')  # remove the time, as it can be a different length
    max_len = 0
    for var_name, input in inputs.items():
        if isinstance(input, list):
            if len(input) > max_len:
                max_len = len(input)
    if max_len == 0:
        max_len = 1

    # Fix the integer inputs and check that the list inputs are all the same length
    for var_name, input in inputs.items():
        if isinstance(input, int) or isinstance(input, str) or isinstance(input, float):
            inputs[var_name] = [input]*(max_len)  # convert to a list of the same integer
        elif isinstance(input, list):
            assert len(input) == max_len, (
                f'The length of {var_name} is {len(input)}, when it should be {max_len}.'
                )
        else:
            raise NotImplementedError(
                f'The variable {var_name} is of the wrong class. It should be an int or a list.'
                )

    all_rxn_idxs = []
    all_sens_arrays = []
    calc_idxs = []
    if parallel:        
        # Set up multiprocessing
        nproc_avail = max(len(os.sched_getaffinity(0)) - 1, 1)
        print('# of processors: ', nproc_avail)
        print('max_len: ', max_len)
        calcs_per_proc = math.ceil( max_len / nproc_avail )
        procs = []
        # Run the sensitivity in parallel
        queue = multiprocessing.Queue()
        for proc_n in range(nproc_avail):
            array_start = proc_n * calcs_per_proc
            if array_start > max_len:
                break
            if proc_n == nproc_avail - 1:                
                array_end = max_len
            else:
                array_end = min((proc_n+1) * calcs_per_proc, max_len)            
            proc = multiprocessing.Process(
                target=_multiple_sens_array,
                args=(queue, inputs, time, array_start, array_end,))

            procs.append(proc)
            proc.start()

        # Collect the results for the multiprocessing queue
        for _ in procs:
             rxn_idxs, sens_arrays, calc_idx_array = queue.get()
             all_rxn_idxs.extend(rxn_idxs)
             all_sens_arrays.extend(sens_arrays)
             calc_idxs.extend(calc_idx_array)

        # Wait for the processors to all finish
        for proc in procs:
            proc.join()

    else:    
        # Run the sensitivity simulations
        calc_idxs = range(max_len)
        for calc in range(max_len):
            rxn_idxs, sens_array = single_sens_array(
                inputs['mech_filename'][calc], inputs['target_spc'][calc],
                inputs['T_K'][calc], inputs['P_atm'][calc], inputs['X'][calc],
                time, inputs['reactor_type'][calc], inputs['num_top_rxns'][calc],
                inputs['target_rxns'][calc]
                )
            all_rxn_idxs.append(rxn_idxs)
            all_sens_arrays.append(sens_array)
    
    # Plot
    plot_sens_arrays(inputs['mech_filename'], inputs['target_spc'], inputs['T_K'],
        inputs['P_atm'], inputs['X'], time, all_rxn_idxs, all_sens_arrays, calc_idxs)

    return all_rxn_idxs, all_sens_arrays


def plot_sens_arrays(mech_filename, target_spc, T_K, P_atm, X, time, all_rxn_idxs, all_sens_arrays, calc_idxs):
    """ Produces plots of one or more sensitivity analyses.

        :param mech_filename: name of the .cti mechanism
        :type mech_filename: str
        :param target_spc: target species for sensitivity
        :type target_spc: str
        :param T_K: initial temperature, in Kelvin
        :type T_K: float
        :param P_atm: initial pressure, in atmospheres
        :type P_atm: float
        :param X: initial species and their mole fractions
        :type X: str, example: 'N2O: 0.002, O2: 0.02, AR: 0.996'
        :param time: array of times for evaluating solution, in seconds
        :type time: 1-D Numpy array
        :param: all_rxn_idxs: the indices of the sensitive reactions
        :type: list[list]
        :param: all_sens_arrays: array of the sensitivities for each rxn vs. time
        :type: list[2-D Numpy array]
    """
    # Loop through each sensitivity calculation
    linestyles = ['-', '--', '-.', ':']  # different line styles
    for idx, calc in enumerate(calc_idxs):#range(len(target_spc)):

        # Get some parameters about the calculations
        num_top_rxns = len(all_rxn_idxs[idx])
        gas = ct.Solution(mech_filename[calc])  # for getting the rxn names

        # Initialize plotting
        fig = plt.figure(figsize=(8,6))  # inches
        ax = fig.add_subplot(111)
        colors = cm.rainbow(np.linspace(0, 1, num_top_rxns))

        # Plot the sensitivity for each of the top rxns
        for colmn in range(num_top_rxns):
            rxn_name = str(gas.reaction(all_rxn_idxs[idx][colmn]))
            myplot, = ax.plot(1000*time, all_sens_arrays[idx][:, colmn], color=colors[colmn], label=rxn_name, linestyle = linestyles[colmn%4])
            title = f'{str(T_K[calc])} K, {str(P_atm[calc]*760)} Torr \n{X[calc]}'

            # Format
            plt.xlabel('Time (ms)', fontsize=12)
            plt.ylabel(target_spc[calc] + ' sensitivity', fontsize=12)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.title(title, fontsize=12)
            plt.tight_layout()
            ax.legend(loc='lower right', fontsize=12)

            # Save
            save_filename = 'dummy' + str(calc) + '.pdf'
            plt.savefig(save_filename)

    # Collate the PDFs
    merger = PdfFileMerger()
    for calc in range(len(target_spc)):
        pdf_filename = 'dummy' + str(calc) + '.pdf'
        merger.append(pdf_filename)
    merger.write('sens_results.pdf')
    merger.close()

    # Delete the dummy PDFs
    for calc in range(len(target_spc)):
        pdf_filename = 'dummy' + str(calc) + '.pdf'
        os.remove(pdf_filename)


def _multiple_sens_array(queue, inputs, time, array_start, array_end):
    """Runs single_sens_array for a group of calculations
    """

    print('Processor {} will run calc {} - {}'.format(os.getpid(), array_start, array_end))
    
    all_rxn_idxs = []
    all_sens_arrays = []
    calc_idxs = []
    for calc in range(array_start, array_end):
        rxn_idxs, sens_array = single_sens_array(
            inputs['mech_filename'][calc], inputs['target_spc'][calc],
            inputs['T_K'][calc], inputs['P_atm'][calc], inputs['X'][calc],
            time, inputs['reactor_type'][calc], inputs['num_top_rxns'][calc],
            inputs['target_rxns'][calc]
            )
        all_rxn_idxs.append(rxn_idxs)
        all_sens_arrays.append(sens_array)
        calc_idxs.append(calc)
    queue.put([all_rxn_idxs, all_sens_arrays, calc_idxs])    

    print('Processor {} is finished'.format(os.getpid()))


def single_sens_array(mech_filename, target_spc, T_K, P_atm, X, time, reactor_type='IG const TP', num_top_rxns=10, target_rxns='all'):
    """ Runs an A-factor sensitivity analysis in Cantera. Returns a list of
        reaction indices and a 2-D array of sensitivity values versus time.

        :param mech_filename: name of the .cti mechanism
        :type mech_filename: str
        :param target_spc: target species for sensitivity
        :type target_spc: str
        :param T_K: initial temperature, in Kelvin
        :type T_K: float
        :param P_atm: initial pressure, in atmospheres
        :type P_atm: float
        :param X: initial species and their mole fractions
        :type X: str, example: 'N2O: 0.002, O2: 0.02, AR: 0.996'
        :param time: array of times for evaluating solution, in seconds
        :type time: 1-D Numpy array
        :param reactor_type: type of Cantera reactor
        :type reactor_type: str
        :param num_top_rxns: the number of most-sensitive reactions to be output
        :type num_top_rxns: float
        :param target_rxns: indices of reactions to consider, or 'all'
        :type target_rxns: list or str
        :return: top_rxn_idxs: the indices of the sensitive reactions
        :rtype: list
        :return: top_sens_results: array of the sensitivities for each rxn vs. time
        :rtype: 2-D Numpy array
    """

    # Set up the mechanism and reactor
    gas = ct.Solution(mech_filename)
    gas.TPX = T_K, 101325*P_atm, X  # pressure is supposed to be in Pa
    if reactor_type == 'IG const TP':
        r = ct.IdealGasConstPressureReactor(gas, energy='off')  # setting energy to 'off' holds T constant
    elif reactor_type == 'IG const PH':
        r = ct.IdealGasConstPressureReactor(gas, energy='on')  # energy is held constant
    else:
        raise NotImplementedError(f'The reactor_type input, {reactor_type},has not been implemented.')

    sim = ct.ReactorNet([r])

    # Enable sensitivity with respect to target reaction(s)
    if target_rxns == 'all':
        target_rxns = range(gas.n_reactions)
    for rxn_idx in target_rxns:
        r.add_sensitivity_reaction(rxn_idx)

    # Set the tolerances for the solution and for the sensitivity coefficients
    sim.rtol = 1.0e-6
    sim.atol = 1.0e-15
    sim.rtol_sensitivity = 1.0e-9
    sim.atol_sensitivity = 1.0e-9

    # Timestep through the solution, evaluating the sensitivity for each rxn
    sens_results = np.zeros((len(time), gas.n_reactions))
    for t_idx, t in enumerate(time):
        sim.advance(t)
        for rxn_idx in target_rxns:
            sens = sim.sensitivity(target_spc, rxn_idx)
            sens_results[t_idx, rxn_idx] = sens

    # Look for the max absolute value for each reaction
    max_sens_vals = np.zeros(sens_results.shape[1])
    for col in range(sens_results.shape[1]):
        max_sens_vals[col] = max(abs(sens_results[:,col]))

    # Get the array of sensitivities for the most-sensitive reactions
    top_rxn_idxs = np.argsort(max_sens_vals)[::-1][:num_top_rxns]  # sorted large to small
    top_sens_results = np.zeros((len(time), num_top_rxns))
    top_rxn_names = []
    for idx in range(num_top_rxns):
        top_sens_results[:, idx] = sens_results[:, top_rxn_idxs[idx]]
        rxn_name = str(gas.reaction(top_rxn_idxs[idx]))
        top_rxn_names.append(rxn_name)

    return top_rxn_idxs, top_sens_results
