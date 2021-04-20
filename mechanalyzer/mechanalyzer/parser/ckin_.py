"""
Read the a CHEMKIN mechanism file
"""

import automol
import chemkin_io
from mechanalyzer.parser._util import get_ich_dct, get_fml


def parse_new(mech_str):
    """ parse
    """
    
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


def parse(mech_str, spc_dct, sort_rxns):
    """ parse a chemkin formatted mechanism file
    """

    # Read the reactions
    rxn_block_str = chemkin_io.parser.mechanism.reaction_block(mech_str)
    rxn_strs = chemkin_io.parser.reaction.data_strings(rxn_block_str)

    # read reactants and products
    reac_and_prods = list(zip(
        map(chemkin_io.parser.reaction.reactant_names, rxn_strs),
        map(chemkin_io.parser.reaction.product_names, rxn_strs)))

    # delete duplicate names
    reac_and_prods = list(set(reac_and_prods))
    rct_names, prd_names = zip(*reac_and_prods)
    rct_names_lst = list(rct_names)
    prd_names_lst = list(prd_names)

    # Build the inchi dct
    ich_dct = get_ich_dct(spc_dct)

    # extract mech info
    formula_dct_lst, formula_str_lst, rxn_name_lst = mech_info(
        rct_names_lst, prd_names_lst, ich_dct)

    # Sort the reactions if desired

    if sort_rxns:
        rxn_info_lst = list(
            zip(formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst))
        fml_dct = dict(zip(formula_str_lst, formula_dct_lst))
        rxn_info_lst.sort()
        formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst = zip(
            *rxn_info_lst)
        formula_dct_lst = list(map(fml_dct.get, formula_str_lst))

    return formula_dct_lst, formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst


# functions called in parse; may be called separately by other scripts
def mech_info(rct_names_lst, prd_names_lst, ich_dct):
    """
    Derives reactants and products formulas and reaction names

    only called by mech.parse_mechanism

    it is called bu mech.build_Dct bu then augmnted by the rxn dct
    """
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
        formula_dct, formula_str = get_fml(rct_ichs)
        formula_dct_lst.append(formula_dct)
        formula_str_lst.append(formula_str)

    return formula_dct_lst, formula_str_lst, rxn_name_lst
