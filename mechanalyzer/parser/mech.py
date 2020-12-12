"""
Functions for mechanism reading and sorting
Making script mechanalyzer/bin/mech.py more compact
"""

import os
import sys
import mechanalyzer
import chemkin_io
import pandas as pd
import numpy as np
import copy

def readfiles(SPCFILE,MECHFILE):
    '''
    read the mechanism and the species files provided by the user
    '''
    with open(SPCFILE, 'r') as file_obj:
        SPC_STR = file_obj.read() 

    # Extract species dictionary
    spc_dct = mechanalyzer.parser.spc.build_spc_dct(SPC_STR,'csv')

    # Read input mechanism file
    with open(MECHFILE, 'r') as file_obj:
        MECH_STR = file_obj.read() 

    # extract rxn block and build reaction parameter dictionary
    units = chemkin_io.parser.mechanism.reaction_units(MECH_STR)
    block_str = chemkin_io.parser.mechanism.reaction_block(MECH_STR)
    rxn_param_dct = chemkin_io.parser.reaction.param_dct(block_str,units[0],units[1])
    # extract elements
    el_block = chemkin_io.parser.mechanism.element_block(MECH_STR)
    elem_tuple = chemkin_io.parser.species.names(el_block)

    return spc_dct, rxn_param_dct, elem_tuple


def build_dct(spc_dct,rxn_param_dct):
    '''
    Build required info for mech sorting
    '''
    # extract info from dictionary:
    # reactants and products
    rcts,prds=zip(*rxn_param_dct.keys())
    rct_names_lst = list(rcts)
    prd_names_lst = list(prds)

    # inchis dictionary
    ich_dct = mechanalyzer.parser.ckin_.get_ich_dct(spc_dct)

    # formulas and reaction names
    formula_dct, formula_str, rxn_name = mechanalyzer.parser.ckin_.mech_info(rct_names_lst,prd_names_lst,ich_dct)

    mech_info = [formula_dct,formula_str, rct_names_lst, prd_names_lst, rxn_name]

    return mech_info


def sort_mechanism(mech_info,spc_dct,SORT_STR,ISOLATE_SPECIES):
    '''
    mech_info: formulas, reaction names
    spc_dct: species dictionary
    SORT_STR: list with sorting criteria
    ISOLATE_SPECIES: list of species you want to isolate in the final mechanism

    calls sorting functions in mechanalyzer/pes
    returns the rxn indices associated with the comments about sorting
    '''
    # call the sorting class
    srt_mch = mechanalyzer.parser.pes.SORT_MECH(mech_info,spc_dct)
    # sort according to the desired criteria
    srt_mch.sort(SORT_STR,ISOLATE_SPECIES)
    # returns the sorted indices and the corresponding comments
    sorted_idx,cmts_dct = srt_mch.return_mech_df()
    return sorted_idx,cmts_dct


def reordered_mech(rxn_param_dct,sorted_idx,cmts_dct):
    '''
    rxn_param_dct: non-sorted reactions
    sorted_idx: indices of the rxn_param_dct in the desired order
    cmts_dct: comments related to new_idx containing the rxn class
    '''
    sorted_val = list(map(rxn_param_dct.get,sorted_idx))
    rxn_param_dct_sorted = dict(zip(sorted_idx,sorted_val))

    return rxn_param_dct_sorted


def write_reordered_mech(elem_tuple,spc_dct,rxn_param_dct_sorted,cmts_dct,sortedmech_name):
    '''
    MECH_STR: full mech to extract elements name
    spc_dct: species dictionary
    rxn_param_dct_sorted: reaction parameters dictionary in the desired order
    cmts_dct: comments dictionary associated with the sorted mechanism
    sortedmech_name: name of the final mech
    '''

    # write
    chemkin_io.writer.mechanism.write_chemkin_file(elem_tuple=elem_tuple,spc_dct=spc_dct,rxn_param_dct=rxn_param_dct_sorted,filename=sortedmech_name,comments=cmts_dct)
