"""
Read the a CHEMKIN mechanism file
"""

import automol
import chemkin_io


def parse(mech_str, spc_dct, sort_rxns):
    """ parse a chemkin formatted mechanism file
    """

    # Read the reactions and participaring species from the mech file
    rxn_block_str = chemkin_io.parser.mechanism.reaction_block(mech_str)
    rxn_strs = chemkin_io.parser.reaction.data_strings(rxn_block_str)
    rct_names_lst = list(
        map(chemkin_io.parser.reaction.reactant_names, rxn_strs))
    prd_names_lst = list(
        map(chemkin_io.parser.reaction.product_names, rxn_strs))

    # Build the inchi dct
    ich_dct = {}
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich_dct[key] = spc_dct[key]['ich']

    # Sort reactant and product name lists by formula to facilitate
    # multichannel, multiwell rate evaluations
    formula_str = ''
    rxn_name_lst = []
    formula_str_lst = []
    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        rxn_name_lst.append(rxn_name)
        rct_ichs = list(map(ich_dct.__getitem__, rct_names))
        formula_dct = ''
        for rct_ich in rct_ichs:
            formula_i_dct = automol.inchi.formula(rct_ich)
            formula_dct = automol.formula.join(formula_dct, formula_i_dct)
        formula_str = automol.formula.string2(formula_dct)
        formula_str_lst.append(formula_str)

    rxn_info_lst = list(
        zip(formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst))

    # Sort the reactions if desired
    if sort_rxns:
        rxn_info_lst.sort()
        formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst = zip(
            *rxn_info_lst)

    return formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst
