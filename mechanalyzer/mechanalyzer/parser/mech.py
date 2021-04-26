"""
Functions for mechanism reading and sorting
Making script mechanalyzer/bin/mech.py more compact
"""

import mechanalyzer
from mechanalyzer.parser import ckin_ as ckin
from mechanalyzer.parser._util import get_ich_dct, get_fml


# OPERATIONS ON MECHANISM OBJECTS
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

    srt_mch = mechanalyzer.parser.sort.SortMech(mech_info, spc_dct)
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

# I/O
def parse_mechanism(mech_str, mech_type, spc_dct):
    """ Get the reactions and species from the mechanism input
    """

    # Parse the info from the chemkin file
    if mech_type == 'chemkin':
        rxn_param_dct, elem_tuple = ckin.parse(mech_str)
    else:
        raise NotImplementedError

    # Build mech_info object to pass along to various functions
    mech_info = _mech_info(rxn_param_dct, spc_dct)

    return rxn_param_dct, mech_info, elem_tuple


def _mech_info(rxn_param_dct, spc_dct):
    """ Build mech_info object for mech sorting

        :param spc_dct: species dictionary
        :type spc_dct: dict[?:?]
        :param rxn_dct: parameter dictionary
        :type rxn_dct: dict[?:?]
        :return mech_info: objects with mech info
        :rtype: list
    """

    def _inf(rct_names, prd_names, ich_dct):
        """ Sort reactant and product name lists by formula to facilitate
            multichannel, multiwell rate evaluations
        """
        rxn_name_lst, formula_str_lst, formula_dct_lst = [], [], []
        for _rct_names, _prd_names in zip(rct_names, prd_names):
            rxn_name = '='.join(['+'.join(_rct_names), '+'.join(_prd_names)])
            rxn_name_lst.append(rxn_name)
            rct_ichs = list(map(ich_dct.__getitem__, _rct_names))
            formula_dct, formula_str = get_fml(rct_ichs)
            formula_dct_lst.append(formula_dct)
            formula_str_lst.append(formula_str)

        return formula_dct_lst, formula_str_lst, rxn_name_lst

    # Extract info from dictionary
    rcts, prds, thrdbdy = zip(*rxn_param_dct.keys())
    rct_names, prd_names, thrdbdy_lst = list(rcts), list(prds), list(thrdbdy)

    # formulas and reaction names (repplace with the mech info from ckin
    ich_dct = get_ich_dct(spc_dct)
    formula_dct, formula_str, rxn_name = _inf(rct_names, prd_names, ich_dct)

    return [formula_dct, formula_str,
            rct_names, prd_names, thrdbdy_lst,
            rxn_name, list(rxn_param_dct.values())]
