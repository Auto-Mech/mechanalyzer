"""
Functions for mechanism reading and sorting
"""

import sys
import autoparse.pattern as app
import ioformat.ptt
from mechanalyzer.parser import ckin_ as ckin


def parse_mechanism(mech_str, mech_type):
    """ Get the reactions and species from the mechanism input
    """

    # Parse the info from the chemkin file
    if mech_type == 'chemkin':
        rxn_param_dct = ckin.parse_rxn_param_dct(mech_str)
    else:
        raise NotImplementedError

    return rxn_param_dct

# Parse the auxiliary file used to sort a mechanism
def parse_sort(sort_str):
    """ Parse the string from the sort.dat input file that contains various
        parameters used to sort a mechanism.

        Returns the list of species to isolate from the mechanism as well
        as the criteria used to sort the mechanism.

        :param sort_str: string for the sort.dat file
        :type sort_str: str
        :rtype: (list(str), list(str), int)
    """

    # Read and format information from the isolate_submech block
    sort_block = ioformat.ptt.end_block(sort_str, 'isolate_submech')
    spc_block = ioformat.ptt.paren_blocks(sort_block, key='species')

    if spc_block:
        spc_lst = list(ioformat.ptt.values_from_block(
            spc_block[0][1], val_ptt=app.one_or_more(app.URLSAFE_CHAR)))
    else:
        spc_lst = []

    # Read and format information from the sort_mech block
    isol_block = ioformat.ptt.end_block(sort_str, 'sort_mech')

    crit_block = ioformat.ptt.paren_blocks(
        isol_block, key='criteria')
    head_block = ioformat.ptt.keyword_value_blocks(
        isol_block, key='n_criteria_headers')

    if crit_block:
        crit_tup = ioformat.ptt.values_from_block(
            crit_block[0][1], val_ptt=app.one_or_more(app.URLSAFE_CHAR))
    else:
        crit_tup = ()
    nhead = int(head_block[0][1]) if head_block is not None else 0

    sort_tup = crit_tup + (nhead,)
    sort_lst = list(sort_tup)

    # Print an error message for an isol block
    if isol_block is None:
        print('*ERROR: sort_mech section is not defined')

    return spc_lst, sort_lst
