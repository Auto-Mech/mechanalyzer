"""
Calls MatPlotLib functionality to create the plot
"""

import numpy as np
from matplotlib import pyplot as plt
import _par as par


# MAKE THE PLOT
def pes(ene_dct, conn_lst):
    """ Make a plot of the PES
        
        :param ene_dct: relative energies for each species
        :type ene_dct: dict[name: energy]
        :param conn_lst: list of all the connections of the PES
        :type conn_lst: lst(str)
    """

    # Generate the coordinates
    spc_coord_dct = _format_coords(ene_dct)

    # Build the plot object
    _plt_axes()
    _plt_species_lines(spc_coord_dct)
    _plt_connecting_lines(spc_coord_dct, conn_lst)
    _plt_ticks()

    # Create the image
    _create_img()


# FORMAT THE DATA TO PLACE OBJECTS ONTO THE PLOT #
def _format_coords(ene_dct,
                   spc_line_spacing=par.SPC_LINE_SPACING,
                   spc_line_base=par.SPC_LINE_BASE,
                   spc_line_len=par.SPC_LINE_LEN):
    """ Generates the coordinates for each of the species in the list.
        Each species is plotted as a horizontal line with some finite length.
    
        :param ene_dct: relative energies for each species
        :type ene_dct: dict[name: energy]
        :param spc_line_spacing: distance between each line
        :type spc_line_spacing: float
        :param spc_line_base: minimal distance between each horizontal line
        :type spc_line_base: float
        :param spc_line_len: length of the each horizontal line
        :type spc_line_len: float
        :return spc_coord_dct: coordinates for each species
        :rtype: dict[name: ((x1, x2), y)]
    """

    spc_coord_dct = {}
    for idx, (name, ene)  in enumerate(ene_dct.items()):
        spc_coord_dct[name] = (
            ((spc_line_spacing * idx + spc_line_base, spc_line_len), ene)

        )

    return spc_coord_dct


# ADD OBJECTS TO THE PES PLOT #
def _plt_axes(y_axis_bot_lim=par.Y_AXIS_BOT_LIM,
              y_axis_top_lim=par.Y_AXIS_TOP_LIM,
              x_axis_right_lim=par.X_AXIS_RIGHT_LIM,
              y_axis_label=par.Y_AXIS_LABEL):
    """ Format the x- and y-axes. 

        :param y_axis_bot_lim: value to extend the bottom part of y-axis
        :type y_axis_bot_lim: float
        :param y_axis_top_lim: value to extend the top part of y-axis
        :type y_axis_top_lim: float
        :param x_axis_right_lim: value ot extend the right part of the x-axis
        :type x_axis_right_lim: float
        :param y_axis_label: label for the y-axis
        :type y_axis_label: float
    """

    plt.axes(frameon=False)
    plt.axvline(0, y_axis_bot_lim, y_axis_top_lim, color='k')
    plt.tick_params(which='both', bottom='off', top='off', right='off',
                    labelbottom='off', labelsize=18)
    plt.xlim(0, x_axis_right_lim)
    plt.ylim(y_axis_bot_lim, y_axis_top_lim)
    plt.ylabel(y_axis_label, fontsize=24)


def _plt_species_lines(spc_coord_dct,
                       name_vshift=NAME_VSHIFT,
                       species_line_color=SPC_LINE_COLOR,
                       species_line_width=SPC_LINE_WIDTH,
                       energy_font_size=ENE_FONT_SIZE,
                       name_font_size=NAME_FONT_SIZE):
    """ Plot the lines that correspond to the molecular species. Each line
        is a horizontal line.
    
        :param spc_coord_dct: coordinates of each species line
        :type spc_coord_dct: dict[name: ((x1, x2), y)]
        :param name_vshift: vertical distance that text appears over line
        :type name_vshift: float
        :param species_line_color: color of the species line
        :type species_line_color: str
        :param species_line_width: (vertical) width of the species line
        :type species_line_width: float
        :param energy_font_size: size of the font of the energy label
        :type energy_font_size: float
        :param name_font_size: size of the font of the name label
        :type name_font_size: float
    """

    for spc_name, ((xp1, xp2), ene) in spc_coord_dct.items():

        # Plot the line for the species
        plt.plot([xp1, xp2], [ene, ene]
                 color=species_line_color,
                 lw=species_line_width, linestyle='-')

        # Plot the label for the line that has the species name
        label_x, label_y = _position_text(xp2, xp1, ene, name_vshift)
        plt.text(label_x, label_y, spc_name,
                 weight='bold', horizontalalignment='center',
                 fontsize=energy_font_size, color='black')


def _plt_connecting_lines(spc_coord_dct, connections):
                          color=CONNE_LINE_COLOR,
                          lw=CONN_LINE_WIDTH,
    """ Plot the lines that connect each of the horizontal species lines that
        show how the species are connected.
    """

    for spc1, spc2 in connections:
        (x1, _), y1 = spc_coord_dct[spc1]
        (_, x2), y2 = spc_coord_dct[spc2]

        plt.plot([x1, x2], [y1, y2], color=color, lw=lw, linestyle='--')


def _plt_ticks(tick_min=TICK_MIN,
               tick_max=TICK_MAX,
               tick_intvl=TICK_INTVL):
    """ Plot the lines that correspond to the molecular species.

        :param tick_min: minimum value on axis to place tick marks
        :type tick_min: float
        :param tick_max: maximum value on axis to place tick marks
        :type tick_max: float
        :param tick_intvl: interval to place tick marks
        :type tick_intvl: float
    """
    plt.yticks(np.arange(tick_min, tick_max, tick_intvl))


def _create_img(width=PLOT_WIDTH, height=PLOT_HEIGHT,
                name=PLOT_NAME, ext=PLOT_EXT, dpi=PLOT_DPI):
    """ Creates the outgoing plot in the format requested by the user.

        :param width: width of the output file image (inches)
        :type width: float
        :param height: height of the output file image (inches)
        :type height: float
        :param name: base name of the output file image
        :type name: str
        :param name: file extension of the output file image
        :type name: str
        :param dpi: dots-per-inch (resolution) of the output file image
        :type dpi: float
    """

    fig = plt.gcf()
    fig.set_size_inches(width, height)
    fig.savefig(name + '.' + ext, dpi=dpi)


# HELPER FUNCTIONS
def _calc_stuff():
    """ Using input from user, determine parameters for plot formatting. """

    name_vshift = name_vshift_scale_fact * range_energy
    energy_vshift = energy_vshift_scale_fact * range_energy
    y_axis_top_lim = max_energy + y_axis_top_extend
    y_axis_bot_lim = min_energy - y_axis_bot_extend
    x_axis_right_lim = ground_species_count + x_axis_right_extend

    if PlotParameter.name_latex_format == 'on':
        for k in range(0, Molecule.species_count):
            tmp = '$' + name[k] + '$'
            name[k] = tmp


def _position_text(xp2, xp1, ene, name_vshift):
    """ Set the position of the text
    """

    # Set the positions of the text
    ene_x = ene + name_vshift
    ene_y = ((xp2 + xp1) / 2.0) - ((xp2 - xp1) * 0.1)
    
    # Format the energy to get proper minus signs
    # ene_txt= _minus_formatter('{0:5.2f}'.format(ene))

    return label_x, label_y


def _minus_formatter(x):
      return unicode(x).replace('-', u'\u2212') 
