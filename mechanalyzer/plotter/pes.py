"""
Calls MatPlotLib functionality to create the plot
"""

import numpy
import matplotlib
from matplotlib import pyplot as plt


# Set plotting options
COLORS = ['k', 'b', 'r', 'g', 'm', 'y', 'c', '#ff9333']
LINESTYLES = ['-', '--', '-.']
MARKERS = ['.', 'o', 's']

# Set various labels for plotting
FIG_TITLE = 'Potential Energy surface'
FONT = {'size': 14}
matplotlib.rc('font', **FONT)

# Parameters for axes set-up
X_AXIS_RIGHT_EXTEND = 0.75
Y_AXIS_TOP_EXTEND = 2.5
Y_AXIS_BOT_EXTEND = 2.5
Y_AXIS_LABEL = 'Relative Enthalpy at 0 K (kcal/mol)'

# Parameters for name and energy labels
NAME_VSHIFT_SCALEF = 0.05
NAME_FONT_SIZE = 18
NAME_LATEX_FORMAT = 'off'
ENE_VSHIFT_SCALEF = 0.05
ENE_FONT_SIZE = 16

# Parameters to draw lines for each species
SPC_LINE_SPACING = 1.15
SPC_LINE_BASE = 0.25
SPC_LINE_LEN = 0.50
SPC_LINE_WIDTH = 4.5
SPC_LINE_COLOR = 'k'
CONN_LINE_WIDTH = 1.5
CONN_LINE_COLOR = 'k'

# Output File Parameters
PLOT_EXT = 'pdf'
PLOT_NAME = 'surface'
PLOT_WIDTH = 16
PLOT_HEIGHT = 9
PLOT_DPI = 1000


# MAKE THE PLOT
def build(ene_dct, conn_lst):
    """ Make a plot of the PES

        :m ene_dct: relative energies for each species
        :type ene_dct: dict[name: energy]
        :m conn_lst: list of all the connections of the PES
        :type conn_lst: lst(str)
    """

    # Generate the coordinates
    spc_coord_dct = _format_coords(ene_dct)

    # Calculate useful information needed for plot
    max_ene, min_ene, spc_cnt = _ranges(ene_dct)
    name_vshift = _calc_vshifts(max_ene, min_ene)
    _, y_axis_tlim, y_axis_blim = _calc_axis_limits(
        max_ene, min_ene, spc_cnt)

    # Build the plot object
    fig, axes = _build_figure()
    _build_axes(axes, y_axis_tlim, y_axis_blim)
    _plt_species_lines(axes, spc_coord_dct, name_vshift)
    _plt_connecting_lines(axes, spc_coord_dct, conn_lst)

    # Create the image
    _create_img(fig)


# FORMAT THE DATA TO PLACE OBJECTS ONTO THE PLOT #
def _format_coords(ene_dct,
                   spc_line_spacing=SPC_LINE_SPACING,
                   spc_line_base=SPC_LINE_BASE,
                   spc_line_len=SPC_LINE_LEN):
    """ Generates the coordinates for each of the species in the list.
        Each species is plotted as a horizontal line with some finite length.

        :m ene_dct: relative energies for each species
        :type ene_dct: dict[name: energy]
        :m spc_line_spacing: distance between each line
        :type spc_line_spacing: float
        :m spc_line_base: minimal distance between each horizontal line
        :type spc_line_base: float
        :m spc_line_len: length of the each horizontal line
        :type spc_line_len: float
        :return spc_coord_dct: coordinates for each species
        :rtype: dict[name: ((x1, x2), y)]
    """

    spc_coord_dct = {}
    for idx, (name, ene) in enumerate(ene_dct.items()):
        xcoord1 = spc_line_spacing * idx + spc_line_base
        xcoord2 = xcoord1 + spc_line_len
        spc_coord_dct[name] = ((xcoord1, xcoord2), ene)

    return spc_coord_dct


# ADD OBJECTS TO THE PES PLOT #
def _build_figure(width=PLOT_WIDTH, height=PLOT_HEIGHT):
    """ Initialize the size and format of the plot figure.

        :m width: width of the output file image (inches)
        :type width: float
        :m height: height of the output file image (inches)
        :type height: float
        :return: figure object for single page of plots
        :rtype: matplotlib.pyplot object
        :return: axes object for single page of plots
        :rtype: matplotlib.pyplot object
    """

    # Initialize plot objects
    fig, axes = plt.subplots(
        nrows=1, ncols=1, figsize=(width, height))

    # Set various plot options
    fig.tight_layout()
    fig.subplots_adjust(left=0.075,
                        top=0.920, bottom=0.075,
                        wspace=0.2, hspace=0.175)

    return fig, axes


def _build_axes(axes_obj, y_axis_tlim, y_axis_blim, tick_intvl=5.0):
    """ Set various things for the axes
    """

    # Set attributes of the y-axis
    axes_obj.set_ylabel(
        'Relative Enthalpy at 0 K (kcal/mol)',
        fontsize='24'
    )

    y_ticks_range = numpy.arange(
        y_axis_blim, y_axis_tlim+tick_intvl, tick_intvl)
    print('tick range', y_ticks_range)
    axes_obj.set_ylim(bottom=y_axis_blim, top=y_axis_tlim)
    axes_obj.set_yticks(y_ticks_range)
    # axes_obj.yaxis.set_major_formatter(
    #     ticker.FormatStrFormatter('%0.2f'))
    axes_obj.set_yticklabels(
        y_ticks_range, fontsize='18'
    )

    # Set attributes of the x-axis
    axes_obj.xaxis.set_ticks([])

    # Remove all of the axis frames except the left one
    axes_obj.spines['bottom'].set_visible(False)
    axes_obj.spines['top'].set_visible(False)
    axes_obj.spines['right'].set_visible(False)


def _plt_species_lines(axes_obj, spc_coord_dct, name_vshift,
                       name_font_size=NAME_FONT_SIZE,
                       species_line_color=SPC_LINE_COLOR,
                       species_line_width=SPC_LINE_WIDTH):
    """ Plot the lines that correspond to the molecular species. Each line
        is a horizontal line.

        :m spc_coord_dct: coordinates of each species line
        :type spc_coord_dct: dict[name: ((x1, x2), y)]
        :m species_line_color: color of the species line
        :type species_line_color: str
        :m species_line_width: (vertical) width of the species line
        :type species_line_width: float
        :m name_vshift: vertical distance that text appears over line
        :type name_vshift: float
        :m name_font_size: size of the font of the name label
        :type name_font_size: float
    """

    for spc_name, ((xp1, xp2), ene) in spc_coord_dct.items():

        # Plot the line for the species
        axes_obj.plot([xp1, xp2], [ene, ene],
                      color=species_line_color,
                      lw=species_line_width, linestyle='-')

        # Plot the label for the line that has the species name
        label_x, label_y = _position_text(xp2, xp1, ene, name_vshift)
        axes_obj.text(label_x, label_y, spc_name,
                      weight='bold', horizontalalignment='center',
                      fontsize=name_font_size, color='black')


def _plt_connecting_lines(axes_obj, spc_coord_dct, conn_lst,
                          color=CONN_LINE_COLOR,
                          line_width=CONN_LINE_WIDTH):
    """ Plot the lines that connect each of the horizontal species lines that
        show how the species are connected.
    """

    for spc1, spc2 in conn_lst:
        (_, xcoord2), ycoord1 = spc_coord_dct[spc1]
        (xcoord1, _), ycoord2 = spc_coord_dct[spc2]

        axes_obj.plot([xcoord2, xcoord1], [ycoord1, ycoord2],
                      color=color, lw=line_width, linestyle='--')


def _create_img(fig,
                width=PLOT_WIDTH, height=PLOT_HEIGHT,
                name=PLOT_NAME, ext=PLOT_EXT, dpi=PLOT_DPI):
    """ Creates the outgoing plot in the format requested by the user.

        :m width: width of the output file image (inches)
        :type width: float
        :m height: height of the output file image (inches)
        :type height: float
        :m name: base name of the output file image
        :type name: str
        :m name: file extension of the output file image
        :type name: str
        :m dpi: dots-per-inch (resolution) of the output file image
        :type dpi: float
    """

    fig_name = '{0}.{1}'.format(name, ext)
    fig.set_size_inches(width, height)
    fig.savefig(fig_name, dpi=dpi)
    plt.close(fig)


# HELPER FUNCTIONS
def _calc_vshifts(max_ene, min_ene,
                  name_vshift_scalef=NAME_VSHIFT_SCALEF):
    """ Using input from user, determine meters for plot formatting.
    """

    energy_range = max_ene - min_ene

    name_vshift = name_vshift_scalef * energy_range
    # ene_vshift = ene_vshift_scalef * energy_range

    return name_vshift
    # return name_vshift, ene_vshift


def _calc_axis_limits(max_ene, min_ene, spc_cnt,
                      x_axis_right_extend=X_AXIS_RIGHT_EXTEND,
                      y_axis_top_extend=Y_AXIS_TOP_EXTEND,
                      y_axis_bot_extend=Y_AXIS_BOT_EXTEND):
    """ Using input from user, determine meters for plot formatting. """

    x_axis_rlim = spc_cnt + x_axis_right_extend
    y_axis_tlim = max_ene + y_axis_top_extend
    y_axis_blim = min_ene - y_axis_bot_extend

    return x_axis_rlim, y_axis_tlim, y_axis_blim


def _ranges(ene_dct):
    """ calc limits of the energy
    """
    enes = ene_dct.values()
    return max(enes), min(enes), len(enes)


def _position_text(xp2, xp1, ene, name_vshift):
    """ Set the position of the text
    """

    # Set the positions of the text
    label_x = (xp2 + xp1) / 2.0
    # label_x = ((xp2 + xp1) / 2.0) - ((xp2 - xp1) * 0.1)  # old
    label_y = ene + name_vshift

    return label_x, label_y


# if PlotParameter.name_latex_format == 'on':
#     for k in range(0, Molecule.species_count):
#         tmp = '$' + name[k] + '$'
#         name[k] = tmp
