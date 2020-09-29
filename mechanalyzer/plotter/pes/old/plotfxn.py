"""
Calls MatPlotLib functionality to create the plot
"""

import numpy as np
from matplotlib import pyplot as plt
from specs import Input
from specs import Molecule
from specs import PlotParameter
from specs import OutFileParameter


def minus_formatter(x):
      return unicode(x).replace('-', u'\u2212') 


def format_axes():
    """ Format the x- and y-axes. Default scaling factors defined in are used, unless user specifies a
        value in the input file.
    """

    plt.axes(frameon=False)
    plt.axvline(0, PlotParameter.y_axis_bot_lim, PlotParameter.y_axis_top_lim, color='k')
    plt.tick_params(which='both', bottom='off', top='off', right='off', labelbottom='off', labelsize=18)
    plt.xlim(0, PlotParameter.x_axis_right_lim)
    plt.ylim(PlotParameter.y_axis_bot_lim, PlotParameter.y_axis_top_lim)
    plt.ylabel(PlotParameter.y_axis_label, fontsize=24)


def generate_xcoords():
    """ Generates a list of the coordinates for the left and right endpoints of each of the horizontal lines
        corresponding to a molecular species.
    """

    for i in range(0, Molecule.ground_species_count):
        tmp1 = (PlotParameter.species_line_spacing * i) + 0.25
        tmp2 = tmp1 + PlotParameter.species_line_length
        Molecule.left_endpt.append(tmp1)
        Molecule.right_endpt.append(tmp2)

    for i in range(0, Molecule.above_species_count):
        tmp1 = Molecule.left_endpt[Molecule.above_state[i] - 1]
        tmp2 = Molecule.right_endpt[Molecule.above_state[i] - 1]
        Molecule.left_endpt.append(tmp1)
        Molecule.right_endpt.append(tmp2)


def plt_spec_lines():
    """ Plot the lines that correspond to the molecular species. """

    for i in range(0, Molecule.species_count):
        mid_line = ( (Molecule.right_endpt[i] + Molecule.left_endpt[i]) / 2 ) - ( (Molecule.right_endpt[i] - Molecule.left_endpt[i]) * 0.1)
        shift1 = Molecule.energy[i] - PlotParameter.energy_vshift
        shift2 = Molecule.energy[i] + PlotParameter.name_vshift

        en = minus_formatter('{0:5.2f}'.format(Molecule.energy[i]))

        plt.plot([Molecule.left_endpt[i], Molecule.right_endpt[i]], [Molecule.energy[i], Molecule.energy[i]],
                 color=PlotParameter.species_line_color, lw=PlotParameter.species_line_width, linestyle='-')
        plt.text(mid_line, shift1, en, weight='bold', horizontalalignment='center',
                 fontsize=PlotParameter.energy_font_size, color='black')
        #plt.text(mid_line, shift2, Molecule.name[i], weight='bold', horizontalalignment='center',
        #         fontsize=PlotParameter.name_font_size, color='black')


        plt.yticks(np.arange(PlotParameter.tick_min, PlotParameter.tick_max, PlotParameter.tick_intvl))


def plt_connecting_lines():
    """ Plot the lines that connect the species lines showing establishing a relationship of two molecular species
        in a reaction mechanism.
    """

    for i in range(0, Molecule.connection_count):
        tmp1 = Molecule.right_endpt[Molecule.left_connection[i] - 1]
        tmp2 = Molecule.left_endpt[Molecule.right_connection[i] - 1]
        tmp3 = Molecule.energy[Molecule.left_connection[i] - 1]
        tmp4 = Molecule.energy[Molecule.right_connection[i] - 1]

        plt.plot([tmp1, tmp2], [tmp3, tmp4],
                 color=PlotParameter.connection_line_color,
                 lw=PlotParameter.connection_line_width,
                 linestyle='--'
        )


def create_pdf():
    """ Creates the outgoing plot in the format requested by the user. """

    fig = plt.gcf()
    fig.set_size_inches(OutFileParameter.width, OutFileParameter.height)
    fig.savefig(OutFileParameter.name + '.' + OutFileParameter.ext, dpi=OutFileParameter.dpi)


