"""
Read the a CHEMKIN mechanism file
"""

import automol
import chemkin_io
import numpy as np


def parse(mech_str, spc_dct,sort_rxns):
    """ parse a chemkin formatted mechanism file
    """

    # Read the reactions and participaring species from the mech file
    rxn_block_str = chemkin_io.parser.mechanism.reaction_block(mech_str)
    rxn_strs = chemkin_io.parser.reaction.data_strings(rxn_block_str)
    rct_names_lst = list(
        map(chemkin_io.parser.reaction.reactant_names, rxn_strs))
    prd_names_lst = list(
        map(chemkin_io.parser.reaction.product_names, rxn_strs))
    
    # delete duplicate names
    rct_names_lst,prd_names_lst = deldup(rct_names_lst,prd_names_lst)

    # Build the inchi dct
    ich_dct = {}
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich_dct[key] = spc_dct[key]['inchi']

    # order the reactants/products by the heaviest
    rct_names_lst,prd_names_lst = order_rct_prd_bystoich(rct_names_lst,prd_names_lst,ich_dct)
    # Sort reactant and product name lists by formula to facilitate
    # multichannel, multiwell rate evaluations
    formula_str = ''
    rxn_name_lst = []
    formula_str_lst = []
    formula_dct_lst = []

    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        rxn_name_lst.append(rxn_name)
        rct_ichs = list(map(ich_dct.__getitem__, rct_names))
        formula_dct = ''
        for rct_ich in rct_ichs:
            formula_i_dct = automol.inchi.formula(rct_ich)
            formula_dct = automol.formula.join(formula_dct, formula_i_dct)
        formula_str = automol.formula.string2(formula_dct)
        formula_dct_lst.append(formula_dct)
        formula_str_lst.append(formula_str)

    rxn_info_lst = list(	    
        zip(formula_dct_lst, formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst))	

    # Sort the reactions if desired	
    if sort_rxns:	
        rxn_info_lst.sort()	
        formula_dct_lst,formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst = zip(	
            *rxn_info_lst)
    return formula_dct_lst, formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst


def order_rct_prd_bystoich(rct_names_lst,prd_names_lst,ich_dct):
    '''
    reorder reactants and products based on the higher number of atoms
    '''

    for key,val in enumerate(rct_names_lst):
        rct_names= val
        rct_ichs = list(map(ich_dct.__getitem__, rct_names))
        fml_rct = list(map(automol.inchi.formula,rct_ichs))
        atoms_rct = list(map(automol.formula.atom_count,fml_rct))
        if len(rct_names)>1:
            if atoms_rct[1] > atoms_rct[0]:
                # swap places of reactants 1 and 2
                rct_names_lst[key] = (rct_names[1],rct_names[0])

    for key,val in enumerate(prd_names_lst):
        prd_names= val
        prd_ichs = list(map(ich_dct.__getitem__, prd_names))
        fml_prd = list(map(automol.inchi.formula,prd_ichs))
        atoms_prd = list(map(automol.formula.atom_count,fml_prd))
        if len(prd_names)>1:
            if atoms_prd[1] > atoms_prd[0]:
                # swap places of reactants 1 and 2
                prd_names_lst[key] = (prd_names[1],prd_names[0])

    return rct_names_lst,prd_names_lst


def deldup(rct_names_lst,prd_names_lst):
    '''
    delete duplicate entries 
    '''

    rct_names_new = []
    prd_names_new = []
    rxn_name_new = []

    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        flagdup = 0
        for rxn in rxn_name_new:
            if rxn == rxn_name:
                flagdup = 1

        if flagdup == 0:
            rxn_name_new.append(rxn_name)
            rct_names_new.append(rct_names)
            prd_names_new.append(prd_names)

    return rct_names_new,prd_names_new
        


