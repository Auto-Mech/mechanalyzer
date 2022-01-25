""" Combines objects pertaining to mechanisms
"""

import copy
from mechanalyzer.calculator import compare


def comb_mechs(rxn_param_dct1, rxn_param_dct2, spc_nasa7_dct1, spc_nasa7_dct2,
               mech_spc_dct1, mech_spc_dct2):
    """ Receives two sets of objects describing two mechanisms (rates, thermo,
        and species) and outputs a single combined mechanism. In the case of
        duplicate values, the first mechanism takes precedence
    """

    rename_instr = compare.get_rename_instr(mech_spc_dct1, mech_spc_dct2)
    comb_rxn_param_dct = comb_dcts(rxn_param_dct1, rxn_param_dct2,
                                   rename_instr, target_type='rxn')
    comb_spc_nasa7_dct = comb_dcts(spc_nasa7_dct1, spc_nasa7_dct2,
                                   rename_instr, target_type='spc')
    comb_mech_spc_dct = compare.get_comb_spc_dct(mech_spc_dct1, mech_spc_dct2)

    return comb_rxn_param_dct, comb_spc_nasa7_dct, comb_mech_spc_dct


def comb_dcts(dct1, dct2, rename_instr, target_type='rxn'):
    """ Combines two dictionaries; can either be rxn_param_dcts or
        spc_nasa7_dcts.
    """

    # Rename the second dictionary to match the first
    renamed_dct2 = compare.rename_species(dct2, rename_instr, target_type)

    # Remove any instances in dct2 that are in dct1
    for key in dct1.keys():
        if target_type == 'rxn':
            matching_rxn, _ = compare.assess_rxn_match(key, renamed_dct2)
            if matching_rxn:
                renamed_dct2.pop(matching_rxn)
        else:  # 'spc'
            if key in renamed_dct2:
                renamed_dct2.pop(key)

    # Add renamed and cleaned up dct2 to the combined dct
    comb_dct = copy.deepcopy(dct1)  # combined dct is initially just dct1
    for key, value in renamed_dct2.items():
        comb_dct[key] = value

    return comb_dct
