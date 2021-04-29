"""
Read the a CHEMKIN mechanism file
"""

import chemkin_io


def parse(mech_str):
    """ parse
    """

    # Extract rxn block and build reaction parameter dictionary
    units = chemkin_io.parser.mechanism.reaction_units(mech_str)
    block_str = chemkin_io.parser.mechanism.reaction_block(mech_str)
    rxn_param_dct = chemkin_io.parser.reaction.param_dct(
        block_str, units[0], units[1])

    # Extract elements if present
    el_block = chemkin_io.parser.mechanism.element_block(mech_str)
    elem_tuple = chemkin_io.parser.species.names(el_block)

    # Read the names of the reactants and products; delete duplicates
    # rxn_block_str = chemkin_io.parser.mechanism.reaction_block(mech_str)
    # rxn_strs = chemkin_io.parser.reaction.data_strings(rxn_block_str)

    # reac_and_prods = list(zip(
    #     map(chemkin_io.parser.reaction.reactant_names, rxn_strs),
    #     map(chemkin_io.parser.reaction.product_names, rxn_strs)))

    # reac_and_prods = list(set(reac_and_prods))
    # rct_names, prd_names = zip(*reac_and_prods)
    # rct_names_lst = list(rct_names)
    # prd_names_lst = list(prd_names)

    # return rct_names_lst, prd_names_lst, rxn_param_dct, elem_tuple
    return rxn_param_dct, elem_tuple
