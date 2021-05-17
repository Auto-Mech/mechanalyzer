"""
Useful functions for mechanalyzer.parser
- mechanism reading
- extract useful info from the mechanism (reactions, formulas..)
"""

import copy
import automol
from automol.formula._formula import element_count as n_el


def order_rct_bystoich(rct_names_lst, spc_dct=None):
    """ Reorder reactants and products based on the higher number of atoms
        If no species dictionary is given:
            reorder just according to name length,
        If length or stoichiometry is the same:
            reorder by alphabetical order
    """
    rct_names_lst_ordered = copy.deepcopy(rct_names_lst)
    ich_dct = {}
    if spc_dct:
        for key in spc_dct.keys():
            if 'ts' not in key and 'global' not in key:
                ich_dct[key] = spc_dct[key]['inchi']

        for key, val in enumerate(rct_names_lst_ordered):
            rct_names = val
            rct_ichs = list(map(ich_dct.__getitem__, rct_names))
            fml_rct = list(map(automol.inchi.formula, rct_ichs))
            atoms_rct = list(map(automol.formula.atom_count, fml_rct))
            if len(rct_names) == 2:
                if atoms_rct[1] > atoms_rct[0]:
                    # swap places of reactants 1 and 2
                    rct_names_lst_ordered[key] = (rct_names[1], rct_names[0])
                elif atoms_rct[1] == atoms_rct[0]:
                    rct_names = list(rct_names)
                    rct_names.sort()
                    rct_names_lst_ordered[key] = tuple(rct_names)

    else:
        # rct_names_lst_ordered.sort()
        for key, val in enumerate(rct_names_lst_ordered):
            rct_names = val
            if len(rct_names) == 2:
                if len(rct_names[1]) > len(rct_names[0]):
                    # swap places of reactants 1 and 2
                    rct_names_lst_ordered[key] = (rct_names[1], rct_names[0])
                elif len(rct_names[1]) == len(rct_names[0]):
                    rct_names = list(rct_names)
                    rct_names.sort()
                    rct_names_lst_ordered[key] = tuple(rct_names)

    return rct_names_lst_ordered


def count_atoms(fml_list):
    """ Count C, O, H atoms in formula list

        :param fml_list: stoich chemical formula
        :type fml_list: list of dictionaries [dict[str:int], ]
        :rtype: list [int, ], int = nCnOnH
    """
    fml_num_list = []
    for fml in fml_list:
        fml_num = (1000*n_el(fml, 'C')+100*n_el(fml, 'O')+n_el(fml, 'H'))
        fml_num_list.append(fml_num)

    return fml_num_list


def extract_spc(spc):
    """
    extract species 1 from tuple
    """
    _spc1, _spc2 = [], []
    for _spc in spc:
        _spc1.append(_spc[0])
        if len(_spc) > 1:
            # bimol species
            _spc2.append(_spc[1])
        else:
            # unimol species
            _spc2.append('')

    return _spc1, _spc2


def get_mult(spc_tuple, spc_dct):
    """
    extracts the total multiplicity of the set of species
    spc_tuple = (A,B,C,..)
    spc_dct: species dictionary
    returns integer of the multiplicity
    """
    mult = 1
    if isinstance(spc_tuple, str):
        spc_tuple = tuple([spc_tuple])

    for spc in spc_tuple:
        mult *= spc_dct[spc]['mult']

    return mult


def get_ich_dct(spc_dct):
    """Generates inchis dictionary from species dictionary
    """
    ich_dct = {}
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich_dct[key] = spc_dct[key]['inchi']

    return ich_dct


def get_fml(rxn_ichs):
    '''
    rxn_icn: inchis of the species of one side of a reaction (ich1, ich2, ..)
    returns: formula dictionary, formula string
    '''
    formula_dct = ''
    for rct_ich in rxn_ichs:
        formula_i_dct = automol.inchi.formula(rct_ich)
        formula_dct = automol.formula.join(formula_dct, formula_i_dct)
    formula_str = automol.formula.string2(formula_dct)

    return formula_dct, formula_str
