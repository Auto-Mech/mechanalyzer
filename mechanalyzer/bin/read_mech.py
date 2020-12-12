import chemkin_io.parser as parser
import sys

MECH_FILENAME = sys.argv[1]



def read_mech_file(mech_filename):
    """ Takes a mechanism file in chemkin format and produces a rxn_param_dct,
        an elements tuple, and a spc_nasa7_dct (this last item only if the thermo
        is included).

    """
    # Get the string for the entire file
    with open(mech_filename) as f:
        mech_str = f.read()
    f.close()

    # Get the strings for each block (skipping species for now...)
    elem_block_str = parser.mechanism.element_block(mech_str)
    rxn_block_str = parser.mechanism.reaction_block(mech_str)
    thermo_block_str = parser.mechanism.thermo_block(mech_str)

    # Process the reaction block
    if rxn_block_str:
        rxn_param_dct = parser.reaction.param_dct(rxn_block_str)
    else: 
        print(f'No reaction block detected in the file {mech_filename}')
        rxn_param_dct = None

    # Process the element block
    if elem_block_str:
        elem_tuple = parser.species.names(elem_block_str)
    else:
        print(f'No elements block detected in the file {mech_filename}')
        elem_tuple = None

    # Process the thermo block
    if thermo_block_str:
        spc_nasa7_dct = parser.thermo.create_spc_nasa7_dct(thermo_block_str)
    else:
        print(f'No thermo block detected in the file {mech_filename}')
        spc_nasa7_dct = None

    return rxn_param_dct, elem_tuple, spc_nasa7_dct
