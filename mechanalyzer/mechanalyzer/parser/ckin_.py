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

    return rxn_param_dct, elem_tuple
