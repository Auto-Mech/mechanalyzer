"""
Functions for mechanism reading and sorting
"""

import sys
import autoparse.pattern as app
import ioformat.ptt
from mechanalyzer.parser import ckin_ as ckin
from mechanalyzer.parser.spc import name_inchi_dct
from mechanalyzer.parser._util import get_fml


def parse_mechanism(mech_str, mech_type, spc_dct):
    """ Get the reactions and species from the mechanism input
    """

    # Parse the info from the chemkin file
    if mech_type == 'chemkin':
        rxn_param_dct = ckin.parse_rxn_param_dct(mech_str)
        elem_tuple = ckin.parse_elem_tuple(mech_str)
    else:
        raise NotImplementedError

    # Build mech_info object to pass along to various functions
    if rxn_param_dct is not None:
        _mech_info = mech_info(rxn_param_dct, spc_dct)
    else:
        _mech_info = None

    return rxn_param_dct, _mech_info, elem_tuple


def mech_info(rxn_param_dct, spc_dct):
    """ Build mech_info object for mech sorting

        :param spc_dct: species dictionary
        :type spc_dct: dict[?:?]
        :param rxn_dct: parameter dictionary
        :type rxn_dct: dict[?:?]
        :return mech_info: objects with mech info
        :rtype: list
    """

    def _check_names(rct_names, prd_names, all_spc_names):
        """ Assess if the reactant and product names provided in the
            rxn_param_dct exist in the spc_dct
        """
        all_mech_names = ()
        for _rct_names, _prd_names in zip(rct_names, prd_names):
            all_mech_names += _rct_names
            all_mech_names += _prd_names
        all_mech_names = set(all_mech_names)

        missing_names = all_mech_names - all_spc_names
        if missing_names:
            print('Names in provided in mechanism, '
                  'but not provided in species list (likely from .csv file):')
            for name in missing_names:
                print('  ', name)
            print('Unable to finish parsing mechanism. Exiting...')
            sys.exit()

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

    # Check if the rxn_param dct may be fully parsed
    _check_names(rct_names, prd_names, set(spc_dct.keys()))

    # formulas and reaction names (repplace with the mech info from ckin
    ich_dct = name_inchi_dct(spc_dct)
    formula_dct, formula_str, rxn_name = _inf(rct_names, prd_names, ich_dct)

    return [formula_dct, formula_str,
            rct_names, prd_names, thrdbdy_lst,
            rxn_name, list(rxn_param_dct.values())]


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
