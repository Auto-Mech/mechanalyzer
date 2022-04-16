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

PEDFLD = ['',] # equivalent to prompt.py but accepts multiple folder inputs for PED and HOTsp
               #returns a single ktp dct merging everything
HOTFLD = ['',] # '' is here
# Parse the command line into an options dictionary
DESC = ('Generate Rate Constants including Prompt Dissociation Effects\n'
        'To Run: python prompt.py')
PAR = argparse.ArgumentParser(description=DESC)
PAR.add_argument('-iped', '--pedinput', default='me_ktp_ped.inp',
                 help='MESS input name (me_ktp_ped.inp)')
PAR.add_argument('-oped', '--pedoutput', default='me_ktp_ped.out',
                 help='MESS ouput name (me_ktp_ped.out)')
PAR.add_argument('-opedmicro', '--pedoutputmicro', default='ke_ped.out',
                 help='MESS microcanonical output name (ke_ped.out)')
PAR.add_argument('-ihot', '--hotinput', default='me_ktp_hoten.inp',
                 help='MESS hotenergy input name (me_ktp_hoten.inp)')
PAR.add_argument('-ohot', '--hotoutput', default='me_ktp_hoten.out',
                 help='MESS hotenergy output name (me_ktp_hoten.out)')
PAR.add_argument('-lhot', '--hotlog', default='me_ktp_hoten.log',
                 help='MESS hotenergy log name (me_ktp_hoten.log)')
PAR.add_argument('-m', '--model', default='equip_simple',
                 help='Statistical energy distribution model (equip_simple)')
PAR.add_argument('-b', '--bf-threshold', default=0.01,
                 help='Minimum branching fraction to include in Prompt (0.01)')
PAR.add_argument('-f', '--fit_method', default='plog',
                 help='method to fit the rates (plog, chebyshev)')

OPTS = vars(PAR.parse_args())

# Set path to current directory where MESS files exist
CWD = os.getcwd()

list_strs_dct = []
# Read the input and output files for MESS calculation of 1st PES
# Here Product Energy Distributions are calculated
for FLD in PEDFLD:
        flddct = dict.fromkeys(['inp','ktp_out','ke_out','ped','log'])
        FLDPATH = os.path.join(CWD, FLD)
        me_ped_inp = pathtools.read_file(FLDPATH, OPTS['pedinput'])
        me_ped_inp = remove_comment_lines(me_ped_inp, delim_pattern=app.escape('!'))
        me_ped_inp = remove_comment_lines(me_ped_inp, delim_pattern=app.escape('#'))
        flddct['inp']  = me_ped_inp 
        flddct['ktp_out']  = pathtools.read_file(FLDPATH, OPTS['pedoutput'])
        flddct['ke_out']  = pathtools.read_file(FLDPATH, OPTS['pedoutputmicro'])
        _, pedoutput = mess_io.reader.ped.ped_names(me_ped_inp)
        flddct['ped'] = pathtools.read_file(FLDPATH, pedoutput)
        list_strs_dct.append(flddct)

for FLD in HOTFLD:
        flddct = dict.fromkeys(['inp','ktp_out','ke_out','ped','log'])
        FLDPATH = os.path.join(CWD, FLD)
        # Read the input and output files for MESS calculation of 2nd PES
        # Here HotEnergies are calculated
        flddct['inp'] = pathtools.read_file(CWD, OPTS['hotinput'])
        flddct['ktp_out'] = pathtools.read_file(CWD, OPTS['hotoutput'])
        flddct['log'] = pathtools.read_file(CWD, OPTS['hotlog'])
        list_strs_dct.append(flddct)

# Get a list of all the statistical energy distribution models
MODELS = OPTS['model'].split(',')
MODELS = ['fne','equip_simple','equip_phi','beta_phi1a','beta_phi2a','beta_phi3a','rovib_dos']
MODELS = ['thermal']
# Generate CKIN files
# This includes both the thermal and non-thermal reactions

# Read and calculate the rate constants into a rxn ktp dictionary
rxn_ktp_dct_models = mechanalyzer.calculator.multipes_prompt_dissociation_ktp_dct(
    list_strs_dct,
    MODELS, OPTS['bf_threshold']
    )

for modeltype in MODELS:
    rxn_param_dct, rxn_err_dct = ratefit.fit.fit_rxn_ktp_dct(
        rxn_ktp_dct_models[modeltype],
        OPTS['fit_method'],
        arrfit_dct={'dbltol': 300}
    )

    # Get the comments dct and write the Chemkin string
    rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
        rxn_err_dct=rxn_err_dct)
    ckin_str = chemkin_io.writer.mechanism.write_chemkin_file(
        rxn_param_dct=rxn_param_dct, rxn_cmts_dct=rxn_cmts_dct)

    # Write the fitted rate parameter file
    pathtools.write_file(ckin_str, CWD, f'rates_{modeltype}.txt')
