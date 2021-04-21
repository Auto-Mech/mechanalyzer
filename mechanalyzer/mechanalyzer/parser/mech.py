"""
Functions for mechanism reading and sorting
Making script mechanalyzer/bin/mech.py more compact
"""

#import os
#import sys
#import copy
import mechanalyzer
import chemkin_io
#import pandas as pd
#import numpy as np


def readfiles(spcfile, mechfile):
    """
    read the mechanism and the species files provided by the user
    :param spcfile: path of csv file for species
    :type spcfile: string
    :param mechfile: path of mech file
    :type mechfile: string
    :return spc_dct: species dictionary
    :rtype spc_dct: dictionary
    :return rxn_param_dct: reaction parameter dictionary
    :rtype rxn_param_dct: dictionary
    :return elem_tuple: elements of the mech
    :rtype elem_tuple: tuple of strings ('el1','el2')
    """
    with open(spcfile, 'r') as file_obj:
        spc_str = file_obj.read()

    # Extract species dictionary
    spc_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, 'csv')

    # Read input mechanism file
    with open(mechfile, 'r') as file_obj:
        mech_str = file_obj.read()

    # extract rxn block and build reaction parameter dictionary
    units = chemkin_io.parser.mechanism.reaction_units(mech_str)
    block_str = chemkin_io.parser.mechanism.reaction_block(mech_str)
    rxn_param_dct = chemkin_io.parser.reaction.param_dct(
        block_str, units[0], units[1])
    # extract elements if present
    try:
        el_block = chemkin_io.parser.mechanism.element_block(mech_str)
        elem_tuple = chemkin_io.parser.species.names(el_block)
    except:
        elem_tuple = None

    return spc_dct, rxn_param_dct, elem_tuple


def build_dct(spc_dct, rxn_dct):
    """
    Build mech_info object for mech sorting
    :param spc_dct: species dictionary
    :type spc_dct: dictionary
    :param rxn_dct: parameter dictionary
    :type rxn_dct_keys: dict
    :return mech_info: objects with mech info
    :rtype: list
    """
    # extract info from dictionary:
    # reactants and products
    rcts, prds, thrdbdy = zip(*rxn_dct.keys())
    rct_names_lst = list(rcts)
    prd_names_lst = list(prds)
    thrdbdy_lst = list(thrdbdy)

    # inchis dictionary
    ich_dct = mechanalyzer.parser.ckin_.get_ich_dct(spc_dct)

    # formulas and reaction names
    formula_dct, formula_str, rxn_name = mechanalyzer.parser.ckin_.mech_info(
        rct_names_lst, prd_names_lst, ich_dct)

    mech_info = [formula_dct, formula_str,
                 rct_names_lst, prd_names_lst, thrdbdy_lst, rxn_name, list(rxn_dct.values())]

    return mech_info

def sorting(mech_info, spc_dct, sort_str, isolate_species):
    """ Build srt_mch object and sort according the desired criteria
    :param mech_info: useful mech info: formulas dct, formulas, reactants and products
                      names, third bodies for each rxn, rxn names, dct parameters
    :type mech_info: list
    :param spc_dct: species dictionary
    :type spc_dct: dictionary
    :param sort_str: sorting criteria
    :type sort_str: list(str)
    :return srt_mch: sorted mechanism
    :rtype: object
    """
    # call the sorting class
    srt_mch = mechanalyzer.parser.sort.SortMech(mech_info, spc_dct)
    # sort according to the desired criteria
    srt_mch.sort(sort_str, isolate_species)

    return srt_mch

def get_sorted_mechanism(srt_mch):
    """ get sorted indexes and comments for a sorted mech object
    :param srt_mch: sorted mechanism
    :type srt_mch: object
    :return sortex_idx, cmts_dct, spc_dct: sorted indexes, dct with comments, species dct
    :rtype: list, dct:str, dct
    """
    # returns the sorted indices and the corresponding comments
    sorted_idx, cmts_dct, spc_dct = srt_mch.return_mech_df()

    return sorted_idx, cmts_dct, spc_dct

def get_sorted_pes_dct(srt_mch):
    """ sort mech info according to the desired criteria and
        get a sorted pes dictionary
    :param srt_mch: sorted mechanism
    :type srt_mch: objectria
    :type sort_str: list(str)
    :return pes_dct: sorted pes dictionary
    :rtype: dct
    """
    # returns the sorted pes dictionary
    pes_dct = srt_mch.return_pes_dct()

    return pes_dct


def reordered_mech(rxn_param_dct, sorted_idx):
    """ Reorder a rxn_param_dct according to the desired sorted indexes

    :param rxn_param_dct: non-sorted reaction parameter dictionary
    :type rxn_param_dct: dct
    :param sorted_idx: indices of the rxn_param_dct in the desired order
    :type sorted_idx: list
    :return rxn_param_dct_sorted: sorted reaction parameter dictionary
    :rtype: dct
    """
    sorted_val = list(map(rxn_param_dct.get, sorted_idx))
    rxn_param_dct_sorted = dict(zip(sorted_idx, sorted_val))

    return rxn_param_dct_sorted


def write_mech(elem_tuple, spc_dct, rxn_param_dct_sorted, sortedmech_name, comments=None):
    """ Write the sorted mechanism - move in writer?

    :param elem_tuple: elements of the mechanism
    :type elem_tuple: tuple
    :param spc_dct: species dictionary
    :type spc_dct: dct
    :param rxn_param_dct_sorted: reaction parameters dictionary in the desired order
    :type rxn_param_dct_sorted: dct
    :param sortedmech_name: name of the sorted mech file
    :type sortedmech_name: string
    :param comments: comments dictionary associated with the sorted mechanism
    :type comments: dct
    """
    # reorder spc_dct before writing to make it nicer
    spc_dct = mechanalyzer.parser.spc.order_species_by_atomcount(spc_dct)
    # write
    chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=elem_tuple, spc_dct=spc_dct, rxn_param_dct=rxn_param_dct_sorted,
        filename=sortedmech_name, comments=comments)
