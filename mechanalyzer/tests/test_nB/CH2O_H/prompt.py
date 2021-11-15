"""
Analyze the extent of prompt dissociation
for a single exothermic reaction with successive decomposition
"""
import os
import numpy as np
import json
import argparse
import mechanalyzer
import mess_io
from autofile.io_ import read_file
from ratefit import ktpdct
import autoparse.pattern as app
from ioformat import remove_comment_lines

##################################################################

# WHEN INTEGRATED IN MECHDRIVER: REMOVE, MUST BE AUTOMATIC
# BF THRESHOLD AS INPUT.. WHERE?
# name of the prompt dissociating species NB must be the same in all PESs
# rad_name is the name of the fragment in the first PES and well in the second
# for other cases like R+O2: the product names can be derived afterwards
rad_name = ['HCO']
bf_threshold = 0.1  # minimum 10% of BF to include the species in the products
modeltype_list = ['beta_phi1a','beta_phi2a','beta_phi3a']  # type of model


def _read_json(path, file_name):
    ''' read a json file with dictionary defined
    '''

    json_path = os.path.join(path, file_name)
    if os.path.exists(json_path):
        with open(json_path, 'r') as fobj:
            json_dct = json.load(fobj)
    else:
        json_dct = None

    return json_dct


# INPUT READING. NB THE I/O WILL BE EDITED AUTOMATICALLY UPON INTEGRATION IN MECHDRIVER
CWD = os.getcwd()
# to launch this: python prompt.py &; default stuff will be searched for automatically.
# Parse the command line
PAR = argparse.ArgumentParser()
PAR.add_argument('-iped', '--pedinput', default='me_ktp_ped.inp',
                 help='MESS input name (me_ktp_ped.inp)')
PAR.add_argument('-oped', '--pedoutput', default='me_ktp_ped.out',
                 help='MESS ouput name (me_ktp_ped.out)')
PAR.add_argument('-opedmicro', '--pedoutputmicro', default='ke_ped.out',
                 help='MESS microcanonical ouput name (ke_ped.out)')
PAR.add_argument('-ihot', '--hotinput', default='me_ktp_hoten.inp',
                 help='MESS hoten input name (me_ktp_hoten.inp)')
PAR.add_argument('-ohot', '--hotoutput', default='me_ktp_hoten.log',
                 help='MESS hoten log name (me_ktp_hoten.log)')
PAR.add_argument('-l', '--label', default='label.inp',
                 help='label dct name (label.inp)') 
                 #optional: to rename the reactions - unnecessary if they come from automech

OPTS = vars(PAR.parse_args())  # it's a dictionary

# READ INITIAL FILES AND LABEL DICTIONARY
label_dct = _read_json(CWD, OPTS['label'])

############ DO NOT MODIFY ##################
# SECONDO ME CONVIENE ALLA FINE METTERE TUTTO IN UNA CLASSE
# OPERATIONS
# 0. EXTRACT INPUT INFORMATION
me_ped_inp = read_file(os.path.join(CWD, OPTS['pedinput']))
me_ped_inp = remove_comment_lines(me_ped_inp, delim_pattern=app.escape('!'))
me_ped_inp = remove_comment_lines(me_ped_inp, delim_pattern=app.escape('#'))
me_ped_out = read_file(os.path.join(CWD, OPTS['pedoutput']))
species_blocks_ped = mess_io.reader.get_species(me_ped_inp)
# print(species_blocks_ped)
T_lst, _ = mess_io.reader.rates.temperatures(me_ped_inp, mess_file='inp')
P_lst, _ = mess_io.reader.rates.pressures(me_ped_inp, mess_file='inp')
P_lst = P_lst[:-1]  # drop the last element in the pressure list ('high')
pedspecies, pedoutput = mess_io.reader.ped.ped_names(me_ped_inp)
# print(pedspecies)

# get the energies
energy_dct, _, conn_lst_dct, _ = mess_io.reader.pes(me_ped_inp)
print(energy_dct, conn_lst_dct)


# ktp dictionary and energy barriers
ktp_dct = {}
barriers_dct = {}
dof_dct = {}
# get the rates for all set of pedspecies
for species in pedspecies:
    reacs, prods = species
    label = '->'.join(species)
    ktp_dct[label] = mess_io.reader.rates.ktp_dct(
        me_ped_out, reacs, prods)
    # find the corresponding energy barrier
    print(reacs, prods)
    barrier_label = mess_io.reader.find_barrier(conn_lst_dct, reacs, prods)
    try:
        barriers_dct[label] = energy_dct[barrier_label]
    except KeyError:
        barriers_dct[label] = max(energy_dct[reacs], energy_dct[prods])
    # derive dofs involved
    dof_info = mechanalyzer.calculator.statmodels.get_dof_info(species_blocks_ped[prods], ask_for_ts=True)
    dof_dct[label] = dof_info

print(barriers_dct)
print(ktp_dct,'\n')

print(pedspecies, energy_dct)
# energies_sp, energies_ts = mess_io.reader.rates.energies(me_ped_out)
# print(energies_sp, energies_ts)
#_, ene_bw = mess_io.reader.rates.barriers(energies_ts, energies_sp, reacs, prods)

# 1. READ THE PEDOUTPUT file and reconstruct the energy distribution
pedoutput_str = read_file(os.path.join(CWD, pedoutput))
ped_dct = mess_io.reader.ped.get_ped(pedoutput_str, pedspecies, energy_dct)

exit()
# 1b. READ THE ke_ped.out file and extract the energy density of each fragment
ke_ped_out = read_file(os.path.join(CWD, OPTS['pedoutputmicro']))
dos_df = mess_io.reader.rates.dos_rovib(ke_ped_out)

# 2. READ THE HOTENERGIES OUTPUT
hot_inp = read_file(os.path.join(CWD, OPTS['hotinput']))
hot_out = read_file(os.path.join(CWD, OPTS['hotoutput']))
species_blocks = mess_io.reader.get_species(hot_inp)
T_lst_hot, _ = mess_io.reader.rates.temperatures(hot_inp, mess_file='inp')
P_lst_hot, _ = mess_io.reader.rates.pressures(hot_inp, mess_file='inp')
# drop the last element in the pressure list ('high')
P_lst_hot = P_lst_hot[:-1]
hotspecies = mess_io.reader.hotenergies.get_hot_names(hot_inp)
hoten_dct = mess_io.reader.hotenergies.extract_hot_branching(
    hot_out, hotspecies, list(species_blocks.keys()), T_lst_hot, P_lst_hot)

# print(ktp_dct)

# derive P_E1

for modeltype in modeltype_list:
    ped_df_prod1 = mess_io.reader.ped.ped_prod1(
        ped_df, rad_name[0], modeltype, dos_df=dos_df, dof_info=dof_info, ene_bw=ene_bw)


    # 3. DERIVE T,P DEPENDENT PRODUCT BRANCHING FRACTIONS and decide which species to keep
    bf_tp_df = mess_io.reader.bf.bf_tp_df_full(rad_name, ped_df_prod1, hoten_dct)
    bf_tp_dct = mess_io.reader.bf.bf_tp_dct_filter(
        bf_tp_df, bf_threshold, modeltype, T_all=T_lst)
    print(bf_tp_dct)

    # 4. DO ARRHENIUS FITS FOR THE SELECTED BFs
    rxn_ktp_dct = mess_io.reader.bf.merge_bf_rates(bf_tp_dct, ktp_dct)
    print(rxn_ktp_dct)
exit()

    # rename the ktp dictionary with appropriate reaction names
rxn_ktp_dct = ktpdct.rename_ktp_dct(rxn_ktp_dct, pedspecies, label_dct)
fitted_dct = fit_ktp_dct(rxn_ktp_dct, CWD)
print(fitted_dct)
#

def fit_ktp_dct(ktp_dct, mess_path):
    """ Fit with plog a given ktp dictionary
        returns strings with fits
    """
    chemkin_str = ''
    for reaction in ktp_dct.keys():

        chemkin_str += arrfit.pes(
            ktp_dct[reaction], reaction, mess_path, dbltol=10000.0)

    return chemkin_str
