""" Runs the sorter
"""

from automol.util.dict_ import filter_keys
from chemkin_io.writer.mechanism import write_chemkin_file
from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import spc as sparser
from mechanalyzer.builder import sort_fct

SPC_TYPE = 'csv'
MECH_TYPE = 'chemkin'


def _sort_pes(spc_name, mech_name, isolate_species, sort_str):
    """ Function that extracts sorted subpes for a mech
    """

    # Read the files
    with open(spc_name, 'r') as spc_obj:
        spc_str = spc_obj.read()
    with open(mech_name, 'r') as mech_obj:
        mech_str = mech_obj.read()

    # Build mech information
    spc_dct_full = sparser.build_spc_dct(spc_str, SPC_TYPE)
    _, mech_info, _ = mparser.parse_mechanism(
        mech_str, MECH_TYPE, spc_dct_full)

    # Sorting: sort the mech and build the sorted rxn param dct
    srt_mch = sorting(
        mech_info, spc_dct_full, sort_str, isolate_species)
    pes_dct_sorted = sorted_pes_dct(srt_mch)

    return pes_dct_sorted


def _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str):
    """ Function that conducts the sorting process for all of the above tests
    """

    # Read the files
    with open(spc_name, 'r') as spc_obj:
        spc_str = spc_obj.read()
    with open(mech_name, 'r') as mech_obj:
        mech_str = mech_obj.read()

    # Build mech information
    spc_dct_full = sparser.build_spc_dct(spc_str, SPC_TYPE)
    rxn_param_dct, mech_info, elems = mparser.parse_mechanism(
        mech_str, MECH_TYPE, spc_dct_full)

    # Sorting: sort the mech and build the sorted rxn param dct
    srt_mch = sorting(
        mech_info, spc_dct_full, sort_str, isolate_species)

    sorted_idx, cmts_dct, spc_dct = sorted_mech(srt_mch)
    rxn_param_dct_sorted = reordered_mech(rxn_param_dct, sorted_idx)

    # Write the new mechanism
    spc_dct_ord = sparser.order_species_by_atomcount(spc_dct)
    mech_str = write_chemkin_file(
        elem_tuple=elems, spc_dct=spc_dct_ord,
        rxn_param_dct=rxn_param_dct_sorted,
        comments=cmts_dct)
    with open(sortmech_name, 'w') as mech1_obj:
        mech1_obj.write(mech_str)

    # If isolated species provided, save remaining reactions in another file
    if isolate_species:
        rxn_param_dct_rest = filter_keys(
            rxn_param_dct, rxn_param_dct_sorted)
        mech_rest_str = write_chemkin_file(
            elem_tuple=elems, spc_dct=spc_dct_full,
            rxn_param_dct=rxn_param_dct_rest)
        with open(mech_rest_name, 'w') as mech2_obj:
            mech2_obj.write(mech_rest_str)


def sorting(mech_info, spc_dct, sort_str, isolate_species):
    """ Uses the SortMech class to sort mechanism info and
        returns the sorted indices and the corresponding comments.

    :param mech_info: formulas, reaction names
    :param spc_dct: species dictionary
    :param sort_str: list with sorting criteria
    :param isolate_species: species you want to isolate in the final mechanism
    :type isolate_species: list()

    calls sorting functions in mechanalyzer/pes
    returns the rxn indices associated with the comments about sorting
    """

    srt_mch = sort_fct.SortMech(mech_info, spc_dct)
    srt_mch.sort(sort_str, isolate_species)

    return srt_mch


def sorted_mech(srt_mch):
    """ get sorted indexes and comments for a sorted mech object
    :param srt_mch: sorted mechanism
    :type srt_mch: object
    :return sortex_idx, cmts_dct, spc_dct:
        sorted indexes, dct with comments, species dct
    :rtype: list, dct:str, dct
    """
    sorted_idx, cmts_dct, spc_dct = srt_mch.return_mech_df()
    return sorted_idx, cmts_dct, spc_dct


def sorted_pes_dct(srt_mch):
    """ sort mech info according to the desired criteria and
        get a sorted pes dictionary
    :param srt_mch: sorted mechanism
    :type srt_mch: objectria
    :type sort_str: list(str)
    :return pes_dct: sorted pes dictionary
    :rtype: dct
    """
    pes_dct = srt_mch.return_pes_dct()
    return pes_dct


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
