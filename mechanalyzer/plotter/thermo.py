"""
Plot the thermochemical parameters from a CHEMKIN mechanism file
"""

import os
import subprocess
import matplotlib.pyplot as plt


# Set plotting options
COLORS = ['k', 'b', 'r', 'g', 'm', 'y']
LINESTYLES = ['-', '--', '-.']
MARKERS = ['.', 'o', 's']

# Set various labels for plotting
AXES_DCTS = [
    {'ylabel': 'Enthalpy (kcal/mol)'},
    {'ylabel': 'Entropy (kcal/mol.K)'},
    {'xlabel': 'Temperature (K)',
     'ylabel': 'Gibbs Free Energy (kcal/mol)'},
    {'xlabel': 'Temperature (K)',
     'ylabel': 'Heat Capacity (kcal/K)'}
]


def build(thermo_dct, temps, dir_prefix='.', names=None):
    """ run over the dictionary for plotting
    """

    # Initialize file string to species and file names
    file_name_str = '{0:40s}{1}\n'.format('Name', 'Filename')

    # Set names to dict values if ther aren't anything
    if names is None:
        names = list(thermo_dct.keys())

    # Make the directory that holds the plots if it doesn't exist
    plot_dir = '{0}/thermo_plots'.format(dir_prefix)
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # Plot the thermo data for each species
    for i, (species, _) in enumerate(zip(thermo_dct, names)):

        # Set the name of the plot and update plot name file string
        file_name = 'spc{0}'.format(str(i+1))
        file_name_str += '{0:40s}{1}\n'.format(species, file_name)

        # build and save the figure to a PDF
        fig, axes = _build_figure(species)
        _build_axes(axes, thermo_dct[species], temps)
        fig_name = '{0}/{1}.pdf'.format(plot_dir, file_name)
        fig.savefig(fig_name, dpi=100)
        plt.close(fig)

    # # Collate all of the pdfs together
    # _collate_pdfs(plot_dir)

    # # Write file relating plot.pdf names to species names
    # with open('plot_names.txt', 'w') as plot_name_file:
    #     plot_name_file.write(file_name_str)


def _build_figure(species):
    """ Initialize the size and format of the plot figure.

        :param nreactions: number of reactions on single page of figure
        :type nreactions: int
        :return: figure object for single page of plots
        :rtype: matplotlib.pyplot object
        :return: axes object for single page of plots
        :rtype: matplotlib.pyplot object
    """

    # Initialize plot objects
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    # Set variables
    species_name = species  # Need to convert back to the chemkin name
    fig_title = 'Comparing Data for {0}'.format(species_name)

    # Set various plot options
    fig.suptitle(fig_title)
    fig.tight_layout()
    fig.subplots_adjust(left=0.075,
                        top=0.95, bottom=0.075,
                        wspace=0.2, hspace=0.075)

    return fig, axes


def _build_axes(axes_obj, species, temps):
    """ plot the rates for various pressures
    """

    # Build mech lists for conveniant passing
    [m1_h, m1_cp, m1_s, m1_g] = species['mech1']
    [m2_h, m2_cp, m2_s, m2_g] = species['mech2']
    mech_therm_lists = [[m1_h, m2_h], [m1_cp, m2_cp],
                        [m1_s, m2_s], [m1_g, m2_g]]

    # Create the thermo plots
    for i in range(2):
        for j in range(2):
            _thermo_axes(axes_obj[i, j], mech_therm_lists[i], temps)
            idx = i+j if i == 0 else i+j+1
            axes_obj[i, j].set(**AXES_DCTS[idx])


def _thermo_axes(ax_obj, mech_therm, temps):
    """ plots onto the axes objects with the thermo data for
        each thermochemical property
    """
    for i, vals in enumerate(mech_therm):
        if vals is not None:
            temps, vals = _trim_vals(temps, vals)
            ax_obj.plot(temps, vals,
                        color=COLORS[i],
                        linestyle=LINESTYLES[0],
                        label='Mech {}'.format(str(i+1)))
    ax_obj.legend(loc='lower right')


def _trim_vals(temps, vals):
    """ trim off values that are undefined
    """
    trim_temps, trim_vals = [], []
    for temp, val in zip(temps, vals):
        if val is not None:
            trim_temps.append(temp)
            trim_vals.append(val)

    return trim_temps, trim_vals


def _collate_pdfs(plot_dir):
    """ collate all of the pdfs together
    """
    plots = [directory for directory in os.listdir(plot_dir)
             if 'pdf' in directory]
    plots.sort(key=lambda x: int(x.replace('spc', '').replace('.pdf', '')))
    plots.append('all_thermo.pdf')
    plots = [os.path.join(plot_dir, name) for name in plots]

    command = ['pdfunite'] + plots
    subprocess.call(command)
