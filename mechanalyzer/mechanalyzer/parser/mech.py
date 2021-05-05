"""
Functions for mechanism reading and sorting
"""

import mechanalyzer
from mechanalyzer.parser import ckin_ as ckin
from mechanalyzer.parser._util import get_ich_dct, get_fml
from ioformat import ptt
import autoparse.find as apf

# inputs to the sorter


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


def read_sort_section(sort_str):
    """ reads the options for sorting from a file

        :param sort_str: string of sorting file w/o comments
        :type sort_str: str
        :return isolate_species: species to include (if []: all)
        :rtype isolate_species: list
        :return sort_list: sorting criteria
        :rtype sort_list: list
    """

    submech_section = apf.all_captures(
        ptt.end_block_ptt('isolate_submech'), sort_str)

    if submech_section is None:
        # empty section
        isolate_species = []
    else:
        # format the section
        species = apf.first_capture(
            ptt.paren_section('species'), submech_section[0])
        isolate_species = ptt.build_keyword_lst(species)

    sortmech_section = apf.all_captures(
        ptt.end_section('sort_mech'), sort_str)
    # this section is mandatory
    if sortmech_section is None:
        print('*ERROR: sort_mech section is not defined')
    else:
        criteria = apf.first_capture(
            ptt.paren_section('criteria'), sortmech_section[0])
        n_criteria = apf.first_capture(
            ptt.keyword_pattern('n_criteria_headers'), sortmech_section[0])
        sort_list = ptt.build_keyword_lst(criteria+n_criteria)

    return isolate_species, sort_list
