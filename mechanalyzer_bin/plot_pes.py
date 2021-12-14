""" Script to read in a PES from a MESS input file and then plot it.

    PLot generated is ReactionCoordinate vs. Relative Energy.
"""

import sys
import os
import argparse
import ioformat
import automol
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
    print(f'ERROR: Input PES file {OPTS["input"]} not found')
    sys.exit()

# Parse the MESS for information required to plot the PES
ene_dct, _, conn_lst_dct, pes_lab_dct = mess_io.reader.pes(
    INP_PES_STR, read_fake=False)
pes_lab_dct = automol.util.dict_.invert(pes_lab_dct)

print('Information parsed from MESS input file')
print('Energies of species:')
for name, ene in ene_dct.items():
    print(f'{name}: {ene} kcal/mol')

print('Connections:')
for conn in conn_lst_dct.items():
    print(f'{conn[0]} - {conn[1]}')

# Try and resort to make plot nice
# ord_ene_dct = mechanalyzer.plotter.pes.resort_names(ene_dct, conn_lst)
# print('fin', list(ord_ene_dct.keys()))

# Produce the PES plot
mechanalyzer.plotter.pes.pes_graph(
    ene_dct, conn_lst_dct, label_dct=pes_lab_dct)

# Exit
print('Plot surface created. Script Execution complete.')
