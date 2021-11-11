"""
Analyze the extent of prompt dissociation
for a single exothermic reaction with successive decomposition
"""
import os
import sys
import json
import argparse
import mechanalyzer
import mess_io
from autofile.io_ import read_file
import autoparse.pattern as app
from ioformat import remove_comment_lines

##################################################################

# edit
bf_threshold = 0.01  # minimum 1% of BF to include the species in the products
modeltype_list = ['beta_phi1a','beta_phi2a','beta_phi3a','equip_phi','equip_simple','rovib_dos']  # type of model

##################################################################

# INPUT READING. NB THE I/O WILL BE EDITED AUTOMATICALLY UPON INTEGRATION IN MECHDRIVER
CWD = os.getcwd()
# to launch this: python prompt.py &; default stuff will be searched for automatically.
# Parse the command line
PAR = argparse.ArgumentParser()
PAR.add_argument('-iped', '--pedinput', default='me_ktp_ped.inp',
                 help='MESS input name (me_ktp_ped.inp)')
PAR.add_argument('-oped', '--pedoutput', default='rate_ped.out',
                 help='MESS ouput name (rate_ped.out)')
PAR.add_argument('-opedmicro', '--pedoutputmicro', default='ke_ped.out',
                 help='MESS microcanonical ouput name (ke_ped.out)')
PAR.add_argument('-ihot', '--hotinput', default='me_ktp_hoten.inp',
                 help='MESS hoten input name (me_ktp_hoten.inp)')
PAR.add_argument('-ohot', '--hotoutput', default='me_ktp_hoten.log',
                 help='MESS hoten log name (me_ktp_hoten.log)')
                 #optional: to rename the reactions - unnecessary if they come from automech

OPTS = vars(PAR.parse_args())  # it's a dictionary

############ DO NOT MODIFY ##################
# SECONDO ME CONVIENE ALLA FINE METTERE TUTTO IN UNA CLASSE
# OPERATIONS
# 0. EXTRACT INPUT INFORMATION from me_ped.inp
me_ped_inp = read_file(os.path.join(CWD, OPTS['pedinput']))
me_ped_inp = remove_comment_lines(me_ped_inp, delim_pattern=app.escape('!'))
me_ped_inp = remove_comment_lines(me_ped_inp, delim_pattern=app.escape('#'))
me_ped_out = read_file(os.path.join(CWD, OPTS['pedoutput']))
species_blocks_ped = mess_io.reader.get_species(me_ped_inp)
T_lst, _ = mess_io.reader.rates.temperatures(me_ped_inp, mess_file='inp')
P_lst, _ = mess_io.reader.rates.pressures(me_ped_inp, mess_file='inp')
P_lst = P_lst[:-1]  # drop the last element in the pressure list ('high')
pedspecies, pedoutput = mess_io.reader.ped.ped_names(me_ped_inp)
energy_dct, _, conn_lst_dct, _ = mess_io.reader.pes(me_ped_inp)

# 1. INFO FROM rate_ped.out and ke_ped.out: rate dct, energy barriers, dofs, fragment names
ktp_dct = {}
E_BW_dct = {}
dof_dct = {}
fragments_dct = {}
# get the rates for all set of pedspecies
for species in pedspecies:
    reacs, prods = species
    label = '->'.join(species)
    ktp_dct[label] = mess_io.reader.rates.ktp_dct(
        me_ped_out, reacs, prods)
    # find the corresponding energy barrier
    barrier_label = mess_io.reader.find_barrier(conn_lst_dct, reacs, prods)
    try:
        E_BW_dct[label] = energy_dct[barrier_label]-energy_dct[prods]
    except KeyError:
        E_BW_dct[label] = energy_dct[reacs]-energy_dct[prods]
    # derive dofs involved
    dof_info = mechanalyzer.calculator.statmodels.get_dof_info(species_blocks_ped[prods], ask_for_ts=True)
    dof_dct[label] = dof_info
    fragments_dct[label] = mess_io.reader.dct_species_fragments(species_blocks_ped)[prods]


# 2. read PED
pedoutput_str = read_file(os.path.join(CWD, pedoutput))
ped_dct = mess_io.reader.ped.get_ped(pedoutput_str, pedspecies, energy_dct)

# 3. READ THE ke_ped.out file and extract the energy density of each fragment
ke_ped_out = read_file(os.path.join(CWD, OPTS['pedoutputmicro']))
dos_df = mess_io.reader.rates.dos_rovib(ke_ped_out)

# 4. READ THE HOTENERGIES OUTPUT
hot_inp = read_file(os.path.join(CWD, OPTS['hotinput']))
hot_out = read_file(os.path.join(CWD, OPTS['hotoutput']))
species_blocks_hoten = mess_io.reader.get_species(hot_inp)
hot_frag_dct = mess_io.reader.dct_species_fragments(species_blocks_hoten)
T_lst_hot, _ = mess_io.reader.rates.temperatures(hot_inp, mess_file='inp')
P_lst_hot, _ = mess_io.reader.rates.pressures(hot_inp, mess_file='inp')
P_lst_hot = P_lst_hot[:-1] #drop last value of pressure
hotspecies = mess_io.reader.hoten.get_hot_names(hot_inp)
hoten_dct = mess_io.reader.hoten.extract_hot_branching(
    hot_out, hotspecies, list(species_blocks_hoten.keys()), T_lst_hot, P_lst_hot)

# DERIVE BF AND RATES
rxns = {}
for species in pedspecies:
    label = '->'.join(species)
    ped_df = ped_dct[label]
    E_BW = E_BW_dct[label]
    # select the fregment of which you want the PED: it is the one in common with hotspecies
    fragments = fragments_dct[label]
    try:
        frag1 = list(set(hotspecies).intersection(fragments))[0]
        fragments.remove(frag1)
        frag2 = fragments[0]
    except IndexError:
        print('no superposition between PED fragments and hot fragments - exiting now \n')
        sys.exit()
    # DERIVE PED OF THE HOT FRAGMENT
    ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
        ped_df, frag1, frag2, modeltype_list, dos_df=dos_df, dof_info=dof_dct[label], E_BW=E_BW)

    # JOIN PED AND HOTEN -> DERIVE PRODUCTS BF
    bf_tp_dct = mechanalyzer.builder.bf.bf_tp_dct(modeltype_list, ped_df_frag1_dct, hoten_dct[frag1], bf_threshold, savefile=True)

    # NEW KTP DICTIONARY
    frag_reacs = mess_io.reader.dct_species_fragments(species_blocks_ped)[species[0]]
    rxn_ktp_dct = mechanalyzer.builder.bf.merge_bf_ktp(bf_tp_dct, ktp_dct[label], frag_reacs, frag1, frag2, hot_frag_dct)
    rxns[label] = rxn_ktp_dct

