import sys

class Input:
    """ Class to store arguments from command line input. All variable values set by argparse in input_proc.py. """

    # Command line arguments
    input_file_name = ''
    output_file_name = ''

    # Variables for storing the input file stuff.
    input_file_content = []
    input_file_length = 0

    # Run program to not make .pyc
    sys.dont_write_bytecode = True            


class Molecule:
    """ Class of variables which stores info about each species used for plotting. """

    # Info about each species input from user.
    number = []
    name = []
    energy = []
    above_state = []
    max_energy = 0
    min_energy = 0

    # Total number of each of the lower and raised species
    species_count = 0
    ground_species_count = 0
    above_species_count = 0

    # Lists that define how each of the species connect to one another
    connection_type = 'default'
    left_connection = []
    right_connection = []
    connection_count = 0

    # Endpoints for the horizontal lines
    left_endpt = []
    right_endpt = []


class PlotParameter:
    """ Class of variables that help format the axes of the plot. """

    # Parameters for axes set-up
    y_axis_top_lim = 0
    y_axis_bot_lim = 0
    y_axis_top_extend = 2.5
    y_axis_bot_extend = 2.5
    y_axis_label = 'Relative Enthalpy at 0 K (kcal/mol)'
    x_axis_right_lim = 0
    x_axis_right_extend = 0.75
    x_axis_label = ''

    # Parameters for name and energy labels
    name_vshift = 0
    name_vshift_scale_fact = 0.1
    name_font_size = 16
    name_latex_format = 'off'
    energy_vshift = 0
    energy_vshift_scale_fact = 0.045
    energy_font_size = 14

    # Parameters to draw lines for each species
    species_line_spacing = 1.15
    species_line_length = 0.50
    species_line_width = 4.5
    species_line_color = 'k'
    connection_line_width = 1.5
    connection_line_color = 'k'

    # Y-Tick Parameters
    tick_min = -50.0
    tick_max = 30.0
    tick_intvl = 5.0
   

class OutFileParameter:
    """ Variables specify what the format of the output file that is generated. """

    ext = 'pdf'
    name = 'surface'
    width = 16
    height = 9
    dpi = 1000


