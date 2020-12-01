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

    # Read input mechanism file
    with open(MECHFILE, 'r') as file_obj:
        MECH_STR = file_obj.read() 

    return SPC_STR, MECH_STR


def build_dct(SPC_STR,MECH_STR):
    '''
    Build required species dictionary and mechanism dictionary
    This is needed for mechanism sorting
    '''
    # species dictionary
    spc_dct = mechanalyzer.parser.spc.build_spc_dct(SPC_STR,'csv')
    # pes dictionary
    mech_info = mechanalyzer.parser.pes.read_mechanism_file(MECH_STR,'chemkin',spc_dct)
    [formulas_dct,formulas, rct_names, prd_names, rxn_names] = mech_info 
    # extract rxn block and build reaction parameter dictionary
    block_str = chemkin_io.parser.mechanism.reaction_block(MECH_STR)
    rxn_param_dct = chemkin_io.parser.reaction.param_dct(block_str)
    # for consistency: replace the keys with rct and prd names re-ordered
    rxn_param_dct = dict(zip(list(zip(rct_names, prd_names)),rxn_param_dct.values()))

    return spc_dct, mech_info, rxn_param_dct


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


def write_reordered_mech(MECH_STR,spc_dct,rxn_param_dct_sorted,cmts_dct,sortedmech_name):
    '''
    MECH_STR: full mech to extract elements name
    spc_dct: species dictionary
    rxn_param_dct_sorted: reaction parameters dictionary in the desired order
    cmts_dct: comments dictionary associated with the sorted mechanism
    sortedmech_name: name of the final mech
    '''
    # extract element block
    el_block = chemkin_io.parser.mechanism.element_block(MECH_STR)
    elem_tuple = chemkin_io.parser.species.names(el_block)
    # write
    chemkin_io.writer.mechanism.write_mech_file(elem_tuple,spc_dct,rxn_param_dct_sorted,filename=sortedmech_name,comments=cmts_dct)
