from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as plt_pdf
import copy
import numpy as np

# Gets rid of the annoying warning about too many figures being open
plt.rcParams.update({'figure.max_open_warning': 0})

def plot_comparisons(combined_rxn_ktp_dct, combined_rxn_em_dct, input_pressures, mech_names=None):

    assert isinstance(input_pressures, list), (
        f'input_pressures is a {type(input_pressures)}, but should be a list.'
    )

    # Get the number of mechanisms
    vals =combined_rxn_ktp_dct.values()
    val_iter = iter(vals)
    first_val = next(val_iter)
    num_mechs = len(first_val)

    # Check/fix the mech_names input
    if mech_names is None:
        mech_names = []
        for i in range(num_mechs):
            mech_names.append('mech' + str(i+1))
    else:
        assert num_mechs == len(mech_names), (
            f"""The number of mechs is {num_mechs}, while the length of the variable 
                mech_names is {len(mech_names)}."""
            )

    # Create various formatting variables
    if len(input_pressures) <= 6:
        colors = ['r', 'b', 'g', 'm', 'c', 'y']
    else:  # account for the unlikely case where there are more than 6 pressures
        colors = cm.rainbow(np.linspace(0, 1, len(input_pressures)))
    formatting = {}
    for idx, pressure in enumerate(input_pressures):
        formatting[pressure] = (colors[idx], str(pressure) + ' atm')
    formatting['high'] = ('k', 'P-indep')
    linestyles = ['-', '--', '-.', ':']
    k_units_dct = {1: '(s$^{-1}$)', 2: '(cm$^3$ mol$^{-1}$ s$^{-1}$)', 3: '(cm$^6$ mol$^{-2}$ s$^{-1}$)', 4: '(cm$^9$ mol$^{-3}$ s$^{-1}$)'}

    # Loop over each reaction
    figs = []  # yummy
    for rxn_name, ktp_dcts in combined_rxn_ktp_dct.items():

        # Get the units of k
        (rcts, prds) = rxn_name
        em = combined_rxn_em_dct[rxn_name]
        if em:  # if there was '+M' in the reaction name
            em_term = 1
        else:
            em_term = 0
        molecularity = len(rcts) + em_term
        if isinstance(rcts, str):  # correct for the case when there's only one reactant
            molecularity = 1 + em_term
        k_units = k_units_dct[molecularity]
        k_label = 'log$_{10}$$k$ ' + k_units

        # Create figures, axes, etc.
        fig = plt.figure(figsize=(8.5, 11))
        ax1 = fig.add_subplot(211)
        plt.xlabel('1000/$T$ (K$^{-1}$)', fontsize=14)
        plt.ylabel(k_label, fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        ax2 = fig.add_subplot(212)
        plt.xlabel('1000/$T$ (K$^{-1}$)', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        num_rates = len([ktp_dct for ktp_dct in ktp_dcts if ktp_dct is not None])

        # Loop over each ktp_dct in this reaction
        first_rates = True
        for mech_idx, ktp_dct in enumerate(ktp_dcts):
            if ktp_dct is not None:  # skip empty entries
                # Loop over each pressure in this ktp_dct
                for pressure, (temps, ktps) in ktp_dct.items():
                    # Formatting
                    (_color, _label) = formatting[pressure]
                    _label += ', ' + mech_names[mech_idx]

                    # Create the plot of the rate constants
                    plot1 = ax1.plot(1000 / temps, np.log10(ktps), label=_label, color=_color,
                                    linestyle=linestyles[mech_idx])

                    # Calculate and plot ratios. If the first time through, save the current ktp_dct as the reference
                    if first_rates:
                        ref_ktp_dct = copy.copy(ktp_dct)
                        ref_mech_idx = mech_idx
                    # Otherwise, calculate the ratio relative to the ref_vals
                    else:
                        if pressure in ref_ktp_dct.keys():
                            (temps, ref_vals) = ref_ktp_dct[pressure]
                            ratio_vals = ktps/ref_vals
    
                            # Check for the case of numbers very close to one
                            close_to_one = np.isclose(ratio_vals, np.ones(len(ratio_vals)), atol=5e-3)
                            if close_to_one.all():
                                ratio_vals = np.ones(len(ratio_vals))

                            # Check for the case of numbers that are very similar but not close to one
                            close_to_each_other = np.isclose(ratio_vals / ratio_vals[0], np.ones(len(ratio_vals)),
                                atol=5e-3)
                            if close_to_each_other.all():
                                ratio_vals = np.ones(len(ratio_vals)) * ratio_vals[0] 

                            # Plot
                            plot2 = ax2.plot(1000 / temps, ratio_vals, label=_label, color=_color,
                                linestyle=linestyles[mech_idx])
                            plt.yscale('log')                                 

                # After looping through all pressures, mark that the first ktp_dct
                # with actual rate values has been looped over
                first_rates = False

        # Do some formatting
        plt.ylabel(f'$k$ ratio relative to {mech_names[ref_mech_idx]}', fontsize=14)
    
        ax1.legend(fontsize=12, loc='upper right')
        if num_rates > 1:
            ax2.legend(fontsize=12, loc='upper right')
        else:
            plt.annotate('For this reaction, there is only one mechanism with $k$ values', (0.05,0.5))
        rxn_name_formatted = format_rxn_name(rxn_name, em)
        fig.suptitle(rxn_name_formatted, x=0.5, y=0.94, fontsize=20)
        figs.append(fig)

    # Produce a PDF
    print('Producing PDF...')
    pdf = plt_pdf.PdfPages("output.pdf")
    for fig in figs:
        pdf.savefig(fig)
    pdf.close()


#def plot_one_rxn_ktp_dct(rxn_ktp_dct, rxn_param_dct):
#    """ Plots a rxn_ktp_dct
#
#    """
#    # Loop over each reaction in the dictionary    
#    for rxn, ktp_dct in rxn_ktp_dct:
#        
#        # Loop over each pressure in the dictionary
#        for pressure, (temps, ktps) in ktp_dct.items():
#        # Plot
# comment
        
 


#def sort_figures


def format_rxn_name(rxn_key, em):
    """ Receives a rxn_key and the corresponding param_vals
        from a rxn_param_dct and writes it to a string that
        the above functions can handle. Adds +M or if
        applicable.
    """
    rcts = list(rxn_key[0])
    prds = list(rxn_key[1])
    # Convert to list if only one species
    if not isinstance(rcts, list):
        rcts = [rcts]
    if not isinstance(prds, list):
        prds = [prds]

    # Write the strings
    for idx, rct in enumerate(rcts):
        if idx == 0:
            rct_str = rct
        else:
            rct_str += ' + ' + rct
    for idx, prd in enumerate(prds):
        if idx == 0:
            prd_str = prd
        else:
            prd_str += ' + ' + prd

    # Add the +M or (+M) text if it is applicable
    if em:
        rct_str += ' + M'
        prd_str += ' + M'

    rxn_name = rct_str + ' = ' + prd_str
    #print('format_rxn_name: ', rxn_name)
    return rxn_name
