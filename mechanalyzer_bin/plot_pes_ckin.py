""" Script to read in a PES from a MESS input file and then plot it.

    PLot generated is ReactionCoordinate vs. Relative Energy.
"""

import sys
import os
import argparse
import ioformat
import automol.util
import mechanalyzer.plotter
import mechanalyzer.parser


# Set useful global variables
CWD = os.getcwd()

DSTR = 'Generates a graph plot of a PES from Chemkin mechanism'
PAR = argparse.ArgumentParser(description=DSTR)
PAR.add_argument('-s', '--species', default='species.csv',
                 help='species file (species.csv)')
PAR.add_argument('-m', '--mechanism', default='mechanism.dat',
                 help='mechanism file (mechanism.dat)')
PAR.add_argument('-p', '--pes', default='1_1',
                 help='PES_SUBPES to plot (1_1)')
PAR.add_argument('-o', '--output', default='surface.pdf',
                 help='name of output plot file (surface.pdf)')
PAR.add_argument('-l', '--label', default='names',
                 help='labels for each species (names, smiles)')
PAR.add_argument('-g', '--image-size', default=400, type=int,
                 help='set the pixel size N (N,N) for structures')
OPTS = vars(PAR.parse_args())

# Parse the input based on the initial type
INP_SPC_STR = ioformat.pathtools.read_file(CWD, 'species.csv')
INP_MECH_STR = ioformat.pathtools.read_file(CWD, 'mechanism.dat')
if INP_SPC_STR is None and INP_MECH_STR is None:
    print('ERROR: Input species.csv or mechanism.dat not found')
    sys.exit()

# Parse the mechanism file
MECH_SPC_DCT = mechanalyzer.parser.spc.build_spc_dct(
    INP_SPC_STR, 'csv')
pes_dct = mechanalyzer.parser.pes.pes_dictionary(
    INP_MECH_STR, 'chemkin', MECH_SPC_DCT, printlog=False)

# Grab the channels of the PES dct that was requested
PES = OPTS['pes']
[PESNUM, SUBPESNUM] = PES.split('_')
PESNUM, SUBPESNUM = int(PESNUM), int(SUBPESNUM)
print(f'Plotting PES-SubPES: {PESNUM}-{SUBPESNUM}')

CHNLS = ()
for (form, pidx, sidx), chnls in pes_dct.items():
    if pidx == PESNUM-1 and sidx == SUBPESNUM-1:
        CHNLS = chnls
        break

SPC_ICHS, SPC_NAMES = (), ()
CONN_LST = ()
for chnl in CHNLS:

    # Obtain the InChIs and SMILES for reagents from the spc dct
    RCTS, PRDS = chnl[1][0], chnl[1][1]
    RCT_ICHS = tuple(MECH_SPC_DCT[rct]['inchi'] for rct in RCTS)
    PRD_ICHS = tuple(MECH_SPC_DCT[prd]['inchi'] for prd in PRDS)
    RCT_SMIS = tuple(MECH_SPC_DCT[rct]['smiles'] for rct in RCTS)
    PRD_SMIS = tuple(MECH_SPC_DCT[prd]['smiles'] for prd in PRDS)

    # Build master species list with inchis
    SPC_ICHS += RCT_ICHS + PRD_ICHS

    if OPTS['label'] == 'names':
        SPC_NAMES += RCTS + PRDS
        CONN_LST += (('+'.join(RCTS), '+'.join(PRDS)),)
    else:
        SPC_NAMES += RCT_SMIS + PRD_SMIS
        CONN_LST += (('+'.join(RCT_SMIS), '+'.join(PRD_SMIS)),)

# Remove redundant smiles and names
SPC_ICHS = automol.util.remove_duplicates_with_order(SPC_ICHS)
SPC_NAMES = automol.util.remove_duplicates_with_order(SPC_NAMES)

print('Connections:')
for conn in CONN_LST:
    print(f'{conn[0]} - {conn[1]}')

# Produce the 2D image plot
img_path = os.path.join(CWD, f'structs_{OPTS["pes"]}.pdf')
automol.inchi.draw_grid(
    SPC_ICHS,
    sub_img_size=(OPTS['image_size'], OPTS['image_size']),
    legends=SPC_NAMES, save_path=img_path)

# Produce the PES plot
mechanalyzer.plotter.pes.pes_graph(
    CONN_LST, ene_dct=None, label_dct=None,
    file_name=f'surface_{OPTS["pes"]}.pdf')
# mechanalyzer.plotter.pes.pes_graph2(
#     CONN_LST, ene_dct=None, label_dct=None,
#     file_name=f'nx_surface_{OPTS["pes"]}.html')

# Exit
print('\nPlot surface and 2D images created. Script Execution complete.')
