""" Runs the sorter
"""

from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import spc as sparser
from mechanalyzer.builder import sort_fct


SPC_TYPE = 'csv'
MECH_TYPE = 'chemkin'


# Functions to take the mechanism strings (may want to further simplify)
def sorted_pes_dct(spc_str, mech_str, isolate_spc, sort_lst):
    """ Function that extracts sorted subpes for a mech
    """

    srt_mch, _, _ = _sort_objs(spc_str, mech_str, sort_lst, isolate_spc)

    return srt_mch.return_pes_dct()


def sorted_mech(spc_str, mech_str, isolate_spc, sort_lst):
    """ Function that conducts the sorting process for all of the above tests
    """

    # Build mech information
    srt_mch, rxn_param_dct, spc_dct_ord = _sort_objs(
        spc_str, mech_str, sort_lst, isolate_spc)

    sorted_idx, cmts_dct, spc_dct_ord = srt_mch.return_mech_df()
    rxn_param_dct_sort = reordered_mech(rxn_param_dct, sorted_idx)
    rxn_param_dct_rest = {}

    return rxn_param_dct_sort, rxn_param_dct_rest, spc_dct_ord, cmts_dct


def _sort_objs(spc_str, mech_str, sort_lst, isolate_spc):
    """ Build the sort-mech object
    """

    # Build mech information
    spc_dct = sparser.build_spc_dct(spc_str, SPC_TYPE)
    rxn_param_dct = mparser.parse_mechanism(
        mech_str, MECH_TYPE)

    # Build the sorted mechanism and species objects
    srt_mch = sorting(rxn_param_dct, spc_dct, sort_lst, isolate_spc)
    spc_dct_ord = sparser.reorder_by_atomcount(spc_dct)

    return srt_mch, rxn_param_dct, spc_dct_ord


# Functions that perform the individual sorting process
def sorting(rxn_param_dct, spc_dct, sort_lst, isolate_species):
    """ Uses the SortMech class to sort mechanism info and
        returns the sorted indices and the corresponding comments.

    :param rxn_param_dct: reaction parameter dictionary
    :param spc_dct: species dictionary
    :param sort_lst: list with sorting criteria
    :param isolate_species: species you want to isolate in the final mechanism
    :type isolate_species: list()

    calls sorting functions in mechanalyzer/pes
    returns the rxn indices associated with the comments about sorting
    """

    srt_mch = sort_fct.SortMech(rxn_param_dct, spc_dct)
    srt_mch.sort(sort_lst, isolate_species)

    return srt_mch


def reordered_mech(rxn_param_dct, sorted_idx):
    """ Sort the reaction parameter dcitionary using the indices from
        sort functions.

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
