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

    # delete duplicate names
    rct_names_lst, prd_names_lst = deldup(rct_names_lst, prd_names_lst)

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
def get_ich_dct(spc_dct):
    """Generates inchis dictionary from species dictionary
    """
    ich_dct = {}
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich_dct[key] = spc_dct[key]['inchi']

    return ich_dct


def mech_info(rct_names_lst, prd_names_lst, ich_dct):
    """
    Derives reactants and products formulas and reaction names
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


def get_fml(rxn_ichs):
    '''
    rxn_icn: inchis of the species of one side of a reaction (ich1, ich2, ..)
    returns: formula dictionary, formula string
    '''
    formula_dct = ''
    for rct_ich in rxn_ichs:
        formula_i_dct = automol.inchi.formula(rct_ich)
        formula_dct = automol.formula.join(formula_dct, formula_i_dct)
    formula_str = automol.formula.string2(formula_dct)

    return formula_dct, formula_str


def deldup(rct_names_lst, prd_names_lst):
    '''
    delete duplicate name
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

    return rct_names_new, prd_names_new
