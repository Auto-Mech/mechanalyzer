"""
Read the mechanism file
"""

import mechanalyzer
from lib.amech_io.parser import ptt



def build_pes_dct(job_path, mech_type,
                  spc_dct, run_obj_dct, sort_rxns=True):
    """Build the PES dct
    """

    # Read the string
    mech_str = ptt.read_inp_str(job_path, MECH_INP, remove_comments='!')

    # Build the total PES dct
    mech_info = mechanalyzer.parser.pes.mechanism_file(
        mech_str, mech_type, spc_dct, sort_rxns=sort_rxns)
    pes_dct = mechanalyzer.parser.pes.build_pes_dct(*mech_info)

    # Build an index dct relating idx to formula
    idx_dct, form_dct = build_pes_idx_dct(pes_dct)

    # Reduce the PES dct to only what the user requests
    if run_obj_dct:
        pesnums = [idx_pair[0] for idx_pair in run_obj_dct]
    else:
        pesnums = [idx_pair[0] for idx_pair in run_obj_dct]
    reduced_pes_dct = reduce_pes_dct_to_user_inp(pes_dct, pesnums)

    # Get a dct for all of the connected channels with the PESs to run
    conn_chnls_dct = mechanalyzer.parser.pes.connected_channels_dct(
        reduced_pes_dct)

    # Form the pes dct that has info formatted to run
    # Get the models in here
    run_pes_dct = pes_dct_w_rxn_lsts(
        reduced_pes_dct, idx_dct, form_dct, conn_chnls_dct, run_obj_dct)

    # Print the channels for the whole mechanism file
    print_pes_channels(pes_dct)

    return run_pes_dct
