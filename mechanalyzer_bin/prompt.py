"""
Analyze the extent of prompt dissociation
for a single exothermic reaction with successive decomposition
"""

import os
import argparse
from ioformat import pathtools, remove_comment_lines
import autoparse.pattern as app
import mess_io
import chemkin_io
import mechanalyzer
import ratefit


FLDs = ['ME2', 'ME3end'] # folder names - ped, hoten, makes no difference
MODEL = 'equip_phi'
# Parse the command line into an options dictionary
DESC = ('Generate Rate Constants including Prompt Dissociation Effects\n'
        'To Run: python prompt.py')
PAR = argparse.ArgumentParser(description=DESC)
PAR.add_argument('-o', '--output', default='mess.out',
                 help='MESS ouput name (mess.out)')
PAR.add_argument('-omicro', '--outputmicro', default='ke.out',
                 help='MESS microcanonical output name (ke.out)')
PAR.add_argument('-l', '--log', default='mess.log',
                 help='MESS energy log name (mess.log)')
PAR.add_argument('-m', '--model', default='equip_simple',
                 help='Statistical energy distribution model (equip_simple)')
PAR.add_argument('-b', '--bf-threshold', default=0.01,
                 help='Minimum branching fraction to include in Prompt (0.1)')
PAR.add_argument('-f', '--fit_method', default='plog',
                 help='method to fit the rates (plog, chebyshev)')

OPTS = vars(PAR.parse_args())

# Set path to current directory where MESS files exist
CWD = os.getcwd()

list_strs_dct = []
# Read the input and output files for MESS calculation of 1st PES
# Here Product Energy Distributions are calculated
for FLD in FLDs:
    flddct = dict.fromkeys(['inp', 'ktp_out', 'ke_out', 'ped', 'log'])
    FLDPATH = os.path.join(CWD, FLD)
    me_inp = pathtools.read_file(FLDPATH, 'mess.inp')
    me_inp = remove_comment_lines(me_inp, delim_pattern=app.escape('!'))
    me_inp = remove_comment_lines(me_inp, delim_pattern=app.escape('#'))
    flddct['inp'] = me_inp
    flddct['ktp_out'] = pathtools.read_file(FLDPATH, OPTS['output'])
    flddct['ke_out'] = pathtools.read_file(FLDPATH, OPTS['outputmicro'])
    _, pedoutput = mess_io.reader.ped.ped_names(me_inp)
    if pedoutput:
        flddct['ped'] = pathtools.read_file(FLDPATH, pedoutput)
    flddct['log'] = pathtools.read_file(FLDPATH, OPTS['log'])
    list_strs_dct.append(flddct)


# Generate CKIN files
# This includes both the thermal and non-thermal reactions

# Read and calculate the rate constants into a rxn ktp dictionary
rxn_ktp_dct = mechanalyzer.builder.multipes_prompt_dissociation_ktp_dct(
    list_strs_dct,
    MODEL, OPTS['bf_threshold']
)

# Fit
rxn_param_dct, rxn_err_dct = ratefit.fit.fit_rxn_ktp_dct(
    rxn_ktp_dct, 'arr', arrfit_dct={'dbltol': 50.},
    pdep_dct={'temps': (500, 1200, 1800), 'tol': 10.,
              'plow': None, 'phigh': None, 'pval': 1.0}
)

# Get the comments dct and write the Chemkin string
rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
    rxn_err_dct=rxn_err_dct)
ckin_str = chemkin_io.writer.mechanism.write_chemkin_file(
    rxn_param_dct=rxn_param_dct, rxn_cmts_dct=rxn_cmts_dct)

# Write the fitted rate parameter file
pathtools.write_file(ckin_str, CWD, f'rates_{MODEL}.txt')
