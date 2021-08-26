""" Script to read in a PES from a MESS input file and then plot it.

    PLot generated is ReactionCoordinate vs. Relative Energy.
"""

import sys
import os
import argparse
import ioformat
import mess_io.reader
import mechanalyzer.plotter


# Set useful global variables
CWD = os.getcwd()

# Parse the command line
PAR = argparse.ArgumentParser()
PAR.add_argument('-i', '--input', default='mess.inp',
                 help='name of input species rate file (mess.inp)')
# PAR.add_argument('-o', '--output', default='surface.pdf',
#                  help='name of outpt plot file (surface.pdf)')
OPTS = vars(PAR.parse_args())

# Read the MEsS file
INP_PES_STR = ioformat.pathtools.read_file(CWD, OPTS['input'])

if INP_PES_STR is None:
    print('ERROR: Input PES file {} not found'.format(OPTS['input']))
    sys.exit()

# Parse the MESS for information required to plot the PES
ene_dct, conn_lst, _ = mess_io.reader.pes(INP_PES_STR, read_fake=False)

print('Information parsed from MESS input file')
print('Energies of species:')
for name, ene in ene_dct.items():
    print('{}: {} kcal/mol'.format(name, ene))

print('Connections:')
for conn in conn_lst:
    print('{} - {}'.format(conn[0], conn[1]))

# Try and resort to make plot nice
# ord_ene_dct = mechanalyzer.plotter.pes.resort_names(ene_dct, conn_lst)
# print('fin', list(ord_ene_dct.keys()))

# Produce the PES plot
mechanalyzer.plotter.new_pes.make_graph(ene_dct, conn_lst)
# mechanalyzer.plotter.pes.build(ene_dct, conn_lst)

# Exit
print('Plot surface created. Script Execution complete.')
