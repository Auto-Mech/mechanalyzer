"""
Plot the rates from a CHEMKIN mechanism file
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Set plotting options
COLORS = ['k', 'b', 'r', 'g', 'm', 'y', 'c', '#ff9333']
LINESTYLES = ['-', '--', '-.']
MARKERS = ['.', 'o', 's']

FONT = {'size': 14}
matplotlib.rc('font', **FONT)


def build(ktp_dct, temps, dir_prefix='.'):
    """ Generates plots of rate constants for all the reactions
        in two mechanisms.

        :param ktp_dct: k(T,P)s at all temps and pressures
        :type: dict[pressure: temps]
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param dir_prefix: path where the plot directory will be built
        :type dir_prefix: str
        :param names: names of each reaction that serve as titles of their plot
        :type names: list(str)
    """

    # Initialize file string to species and file names
    file_name_str = '{0:40s}{1}\n'.format('Name', 'Filename')

    # Make the directory that holds the plots if it doesn't exist
    plot_dir = '{0}/rate_plots'.format(dir_prefix)
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # Plot the rate constants for each reaction
    reactions = list(ktp_dct.keys())
    nreactions = len(reactions)
    print('numreactions', len(reactions))
    for i in range(0, len(reactions), 4):
        print('idx', i)

        # Determine if plot will have four reactions in it, or fewer
        if i+4 <= nreactions:
            n_plot_reactions = 4
        else:
            n_plot_reactions = nreactions - i
        print('nplots', n_plot_reactions)

        # Create the figure object
        fig, axes = _build_figure()

        for j in range(n_plot_reactions):
            # Set the axes object containing the plotted data for each reaction
            reaction_names = []
            # Determine the reaction dictionaries
            reaction = reactions[i+j]
            reaction_names.append(reaction)
            # Set variables needed for the plotting
            isbimol = _is_bimolecular(reaction)

            # Build the axes objects containing the plotted rate constants
            if j == 0:
                row, col = 0, 0
            elif j == 1:
                row, col = 0, 1
            elif j == 2:
                row, col = 1, 0
            elif j == 3:
                row, col = 1, 1
            bottom = bool(row == 1)
            axes_block = axes[row, col]
            _build_axes(axes_block, reaction, ktp_dct,
                        isbimol, temps, bottom)

        # Set the name of the plot
        file_name = 'r{0}'.format(str(i))
        file_name_str += '{0:40s}{1}\n'.format('reaction', file_name)

        # build and save the figure to a PDF
        fig_name = '{0}/{1}.pdf'.format(plot_dir, file_name)
        fig.savefig(fig_name, dpi=100)
        plt.close(fig)


def _build_figure():
    """ Initialize the size and format of the plot figure.

        :param nreactions: number of reactions on single page of figure
        :type nreactions: int
        :return: figure object for single page of plots
        :rtype: matplotlib.pyplot object
        :return: axes object for single page of plots
        :rtype: matplotlib.pyplot object
    """

    # Initialize plot objects
    # if nreactions == 4:
    fig, axes = plt.subplots(
        nrows=2, ncols=2, figsize=(12, 8))
    #     grid = {'width_ratios': [0.5]}
    #     fig, axes = plt.subplots(
    #         nrows=2, ncols=1, figsize=(12, 8), gridspec_kw=grid)

    # Set various plot options
    fig.tight_layout()
    fig.subplots_adjust(left=0.075,
                        top=0.920, bottom=0.075,
                        wspace=0.2, hspace=0.175)

    return fig, axes


def _build_axes(ax_block, reaction, ktp_dct, isbimol, temps, bottom):
    """ plot the rates for various pressures
        certain checks are made throughout to deal with plotting
        only one reaction on a page

        :param ax_block: column of the axis
        :type ax_block: matplotlib.pyplot axis object
        :param isbimol: signal reaction is bimolecular
        :type isbimol: bool
        :param temps: Temperatures (K)
        :type temps: numpy.ndarray
    """

    # Obtain a list of the pressures and sort from low to high pressure
    reaction_pressures_lst = _get_sorted_pressures(ktp_dct[reaction])
    # print(reaction_mech_dcts)
    # print(reaction_mech_dcts)
    # print('plst', reaction_pressures_lst)
    # print('punion', reaction_pressures_union)

    # Plot the data, setting formatting options for the axes
    _full_plot(ax_block, ktp_dct[reaction], reaction_pressures_lst, temps)
    ax_block.set(**_set_axes_labels(reaction, isbimol, bottom=bottom))


def _full_plot(ax_obj, rxn_ktp_dct, rxn_pressures, temps):
    """ Place data points corresponding to all of the rate constants
        from the two mechanisms on a plot.

        :param ax_obj: axes onject to put points on
        :type ax_obj: matplotlib.pyplot object
        :param mech_ktp_dcts: rate constants for two mechanisms
        :type mech_ktp_dcts: list(dict[rxn: dict[pressure: rate constant]])
        :param mech_pressures: Pressures (atm)
        :type mech_pressures: list(float)
        :param temps: Temperatures (K)
        :type temps: list(float)
    """
    for i, pressure in enumerate(rxn_pressures):
        plab = pressure if pressure != 'high' else 'PIndep'
        if pressure in rxn_ktp_dct:
            ax_obj.plot((1000.0/temps), np.log10(rxn_ktp_dct[pressure]),
                        color=COLORS[i], linestyle=LINESTYLES[0],
                        label=plab)
    ax_obj.legend(loc='upper right')


def _get_sorted_pressures(unsorted_pressures):
    """ get a sorted list of pressures for the reaction
    """
    if unsorted_pressures != ['high']:
        pressures = [pressure for pressure in unsorted_pressures
                     if pressure != 'high']
        pressures.sort()
        if 'high' in unsorted_pressures:
            pressures.append('high')
    else:
        pressures = ['high']

    return pressures


def _is_bimolecular(reaction):
    """ Determines if a reaction is bimolecular
    """
    reactants = reaction[0]
    isbimol = bool(len(reactants) == 2)
    return isbimol


def _set_axes_labels(reaction, isbimol, bottom):
    """ alter the axes dictionary
    """
    if isbimol:
        units = 'cm3/s'
    else:
        units = '1/s'

    axes_dct = {}
    if bottom:
        axes_dct['xlabel'] = '1000/T (1/K)'
    axes_dct['ylabel'] = 'log10 k({0})'.format(units)
    axes_dct['title'] = _set_block_title(reaction)

    return axes_dct


def _set_block_title(reaction):
    """ Update the string for the figure title
    """
    side_lst = []
    for side in reaction:
        side_lst.append('+'.join(side))

    title = '{0:^60s}'.format(
        side_lst[0]+'='+side_lst[1])

    return title
