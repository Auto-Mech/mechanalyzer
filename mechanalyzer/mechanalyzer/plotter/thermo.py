from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as plt_pdf
from matplotlib.ticker import FormatStrFormatter
import copy
import numpy as np

# Gets rid of the annoying warning about too many figures being open
plt.rcParams.update({'figure.max_open_warning': 0})

def plot_comparisons(aligned_spc_thermo_dct, mech_names=None, sort_method='difference'):

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

    linestyles = ['-', '--', '-.', ':']

    # Loop over each species
    figs = []  # yummy
    largest_ratios = [0]*num_rxns
    for spc_idx, (spc_name, thermo_dcts) in enumerate(aligned_spc_thermo_dct.items()):

        # Create figures, axes, etc.
        fig = plt.figure(figsize=(8.5, 11))
        ax1 = fig.add_subplot(211)
        plt.xlabel('1000/$T$ (K$^{-1}$)', fontsize=14)
        plt.ylabel(k_label, fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
        ax2 = fig.add_subplot(212)
        plt.xlabel('1000/$T$ (K$^{-1}$)', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))

        # The number of rates that actually exist
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
                            plot2 = ax2.plot(1000 / temps, np.log10(ratio_vals), label=_label, color=_color,
                                linestyle=linestyles[mech_idx])

                            # Get and store the largest ratios 
                            print('max\n', max(ratio_vals))
                            print('1/min\n', 1/min(ratio_vals))
                            if max(ratio_vals) > 1/min(ratio_vals):
                                largest_ratio = max(ratio_vals)
                            else:
                                largest_ratio = 1/min(ratio_vals)
                            if largest_ratio > largest_ratios[rxn_idx]:
                                largest_ratios[rxn_idx] = largest_ratio

                # After looping through all pressures, mark that the first ktp_dct
                # with actual rate values has been looped over
                first_rates = False

        # Do some formatting
        lower_ylabel = 'log$_{10}$ of $k$ ratio relative to' + f' {mech_names[ref_mech_idx]}'
        plt.ylabel(lower_ylabel, fontsize=14)
        ax1.legend(fontsize=12, loc='upper right')
        if num_rates > 1:
            ax2.legend(fontsize=12, loc='upper right')
        else:
            plt.annotate('For this reaction, there is only one mechanism with $k$ values', (0.05,0.5))
        rxn_name_formatted = format_rxn_name(rxn_name, em)
        fig.suptitle(rxn_name_formatted, x=0.5, y=0.94, fontsize=20)

        figs.append(fig)

    # Perform sorting if indicated
    largest_ratios = [i for i in largest_ratios if i != 0]  # remove the entries that are zero
    if sort_method == None:
        sorted_figs = figs
    else:
        if sort_method == 'difference':
            fig_indices = np.argsort(largest_ratios)[::-1]  # the [::-1] flips it to be from largest to smallest
            print('fig_indices\n', fig_indices)
            sorted_figs = []
            for idx in fig_indices:
                sorted_figs.append(figs[idx])
        elif sort_method == 'value': 
            pass  # will add this later...

    # Produce a PDF
    print('Producing PDF...')
    pdf = plt_pdf.PdfPages("output.pdf")
    for fig in sorted_figs:
        pdf.savefig(fig)
    pdf.close()
    
    print('largest_ratios\n', largest_ratios)


