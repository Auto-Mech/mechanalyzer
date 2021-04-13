import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.backends.backend_pdf as plt_pdf
import numpy as np
import pandas as pd
from chemkin_io.writer import _util as writer

K_UNITS_DCT = {1: '(s$^{-1}$)', 2: '(cm$^3$ mol$^{-1}$ s$^{-1}$)',
               3: '(cm$^6$ mol$^{-2}$ s$^{-1}$)', 4: '(cm$^9$ mol$^{-3}$ s$^{-1}$)'}


def build_plots(aligned_rxn_ktp_dct, pressures, filename=None, mech_names=None, sort_method=None):


    # Get the number of mechanisms
    vals = aligned_rxn_ktp_dct.values()
    val_iter = iter(vals)
    first_val = next(val_iter)
    num_mechs = len(first_val)

    # Check/fix the mech_names input
    if mech_names is None:
        mech_names = []
        for i in range(num_mechs):
            mech_names.append('mech' + str(i + 1))
    else:
        assert num_mechs == len(mech_names), (
            f"""The number of mechs is {num_mechs}, while the length of the variable
                   mech_names is {len(mech_names)}."""
        )

    # Get the aligned_rxn_ratio_dct
    aligned_rxn_ratio_dct = get_aligned_rxn_ratio_dct(aligned_rxn_ktp_dct)
    #print('aligned_rxn_ratio_dct in build_plots:\n', aligned_rxn_ratio_dct)

    # Get the format_dct
    format_dct = get_format_dct(pressures)

    # Loop over each rxn and plot
    figs = []
    for rxn, ktp_dcts in aligned_rxn_ktp_dct.items():
        #print(rxn)
        #print(ktp_dcts)
        #print(aligned_rxn_ratio_dct)
        ratio_dcts = aligned_rxn_ratio_dct[rxn]
        # Hacky fix for now
        for ratio_dct in ratio_dcts:
            if ratio_dct is not None:
                for pressure, (temps, ratios) in ratio_dct.items():
                    if max(ratios) > 10 or min(ratios) < 0.1:
                        print(f'Ratio greater than 10x for reaction {rxn}')
        for ktp_dct in ktp_dcts:
            if ktp_dct is not None:
                for pressure, (temps, kts) in ktp_dct.items():
                    if max(kts) > 1e15:
                        print(f'Rate constant greater than 1E+15 for reaction {rxn}')
                        print(temps)
                        print(kts)
    


        molecularity = get_molecularity(rxn)
        fig, axs = build_fig_and_axs(molecularity, ratio_dcts, mech_names)
        fig = plot_single_rxn(rxn, ktp_dcts, ratio_dcts, fig, axs, mech_names, format_dct)
        figs.append(fig)

    # Build the pdf
    if filename is None:
        filename = 'rate_comparison.pdf'
    build_pdf(figs, filename)


def plot_single_rxn(rxn, ktp_dcts, ratio_dcts, fig, axs, mech_names, format_dct):


    linestyles = ['-', '--', '-.', ':']
    for mech_idx, ktp_dct in enumerate(ktp_dcts):
        if ktp_dct is not None:
            for pressure, (temps, kts) in ktp_dct.items():
                (_color, _label) = format_dct[pressure]
                _label += ', ' + mech_names[mech_idx]

                if min(kts) < 0:
                    print(f'Negative rate constant for {rxn} in {mech_names[mech_idx]}')
                    print(f'Pressure (atm): {pressure}')
                    print(temps, kts)
                # Plot the rate constants
                plot1 = axs[0].plot(1000 / temps, np.log10(kts), label=_label, color=_color,
                                 linestyle=linestyles[mech_idx])

                # Plot the ratios if they exist
                ratios_plotted = False
                if ratio_dcts[mech_idx] is not None:  # putting second "if" below prevents errors
                    if ratio_dcts[mech_idx][pressure] is not None:
                        (_, ratios) = ratio_dcts[mech_idx][pressure]
                        ratios_plotted = True
                        plot2 = axs[1].plot(1000 / temps, np.log10(ratios), label=_label,
                                            color=_color, linestyle=linestyles[mech_idx])
    # Add legend
    axs[0].legend(fontsize=12, loc='upper right')
    if ratios_plotted:
        axs[1].legend(fontsize=12, loc='upper right')

    rxn_name_formatted = writer.format_rxn_name(rxn) 
    fig.suptitle(rxn_name_formatted, x=0.5, y=0.94, fontsize=20)

    return fig


def sort_by_rate_constant(aligned_rxn_ktp_dct):

    sorting_vals = {}
    for rxn_idx, (rxn, ktp_dcts) in enumerate(aligned_rxn_ktp_dct.items()):
        #print(rxn_idx)
        #print(rxn)
        max_val = 0
        molecularity = get_molecularity(rxn)
        for ktp_dct in ktp_dcts:
            if ktp_dct is not None:
                for pressure, (_, kts) in ktp_dct.items():
                    if max(kts) > max_val:
                        max_val = max(kts)
        # Store values
        sorting_vals[rxn] = [max_val, molecularity, rxn_idx]

    # Convert to dataframe
    df = pd.DataFrame.from_dict(sorting_vals)
    df = df.transpose()
    #print('df:\n', df)


    # df.sort_values()


def get_aligned_rxn_ratio_dct(aligned_rxn_ktp_dct):


    aligned_rxn_ratio_dct = {}
    for rxn, ktp_dcts in aligned_rxn_ktp_dct.items():
        #print('rxn:\n', rxn)
        #print('ktp_dcts:\n', ktp_dcts)
        ref_ktp_dct = ktp_dcts[0]
        ratio_dcts = []
        for mech_idx, ktp_dct in enumerate(ktp_dcts):
            #print('mech_idx:\n', mech_idx)
            #print('ktp_dct:\n', ktp_dct)
            # If (1) on the first ktp_dct, (2) the ref_ktp_dct is None, or (3) the current_ktp_dct
            # is None, set the ratio_dct to None
            if mech_idx == 0 or ref_ktp_dct is None or ktp_dct is None:
                ratio_dct = None
                #print('first if statement, setting to None')
            # Otherwise, calculate the ratio_dct
            else:
                #print('else statement')
                ratio_dct = {}
                for pressure, (temps, kts) in ktp_dct.items():
                    # If the pressure is defined in the ref ktp_dct, calculate and store the ratio
                    if pressure in ref_ktp_dct.keys():
                        #print(pressure)
                        _, ref_kts = ref_ktp_dct[pressure]
                        ratios = kts / ref_kts
                        ratio_dct[pressure] = (temps, ratios)
                if ratio_dct == {}:  # account for the case when no pressures contain ratios
                    ratio_dct = None
                    #print('second if statement, setting to None')

            # Append the current ratio_dct
            ratio_dcts.append(ratio_dct)
            #print('current ratio dcts:\n', ratio_dcts)
        aligned_rxn_ratio_dct[rxn] = ratio_dcts
        #print('aligned_rxn_ratio_dct:\n', aligned_rxn_ratio_dct)

    return aligned_rxn_ratio_dct


def get_molecularity(rxn):
    """ Get the molecularity of a reaction

    :param rxn:
    :type rxn:
    :return molecularity: the molecularity of the reaction; 1 = unimolec., 2 = bimolec., etc.
    :rtype: int
    """
    rcts, _, third_bods = rxn
    num_rcts = len(rcts)
    third_bod = third_bods[0]
    if third_bod is not None:
        if third_bod[0:2] == '(+':  # if the third body does not count toward molecularity
            molecularity = num_rcts
        else:  # if the third body counts as another species
            molecularity = num_rcts + 1
    # If no third body
    else:
        molecularity = num_rcts

    return molecularity


def get_format_dct(pressures):
    """ Set up the formatting dictionary that describes colors and labels

        :param pressures: pressures at which calculations were performed
        :type pressures: list
        :return format_dct: dct containing color and label for each pressure
        :rtype: dct {pressure1: (color1, label1), pressure2: ...}
    """
    if len(pressures) <= 6:
        colors = ['r', 'b', 'g', 'm', 'c', 'y']
    else:  # account for the unlikely case where there are more than 6 pressures
        colors = cm.rainbow(np.linspace(0, 1, len(pressures)))
    format_dct = {}
    for idx, pressure in enumerate(pressures):
        format_dct[pressure] = (colors[idx], str(pressure) + ' atm')
    format_dct['high'] = ('k', 'P-indep')

    return format_dct


def build_fig_and_axs(molecularity, ratio_dcts, mech_names):



    # Set up the y labels for the rate plot and the ratios plot
    k_units = K_UNITS_DCT[molecularity]
    k_label = 'log$_{10}$$k$ ' + k_units
    ratio_label = ''
    for mech_idx, ratio_dct in enumerate(ratio_dcts):
        if ratio_dct is not None:
            ref_mech_name = mech_names[mech_idx-1]
            ratio_label = 'log$_{10}$ of $k$ ratio relative to' + f' {ref_mech_name}'
            break

    # Set up the figure and axes
    fig = plt.figure(figsize=(8.5, 11))
    axs = []
    axs.append(fig.add_subplot(211))
    plt.xlabel('1000/$T$ (K$^{-1}$)', fontsize=14)
    plt.ylabel(k_label, fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    axs[0].yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    axs.append(fig.add_subplot(212))
    plt.xlabel('1000/$T$ (K$^{-1}$)', fontsize=14)
    plt.ylabel(ratio_label, fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    axs[1].yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))

    return fig, axs


def build_pdf(figs, filename):

    print('Producing PDF...')
    pdf = plt_pdf.PdfPages(filename)
    for fig in figs:
        pdf.savefig(fig)
    pdf.close()

