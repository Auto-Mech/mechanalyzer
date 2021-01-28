"""
Functions for
- submechanism identification
- broad rxn classification
"""

# import os
# import mechanalyzer
# import chemkin_io
import sys
import numpy as np
import pandas as pd
import automol
from mechanalyzer.parser import util


def species_subset(fuel, spc_dct):
    """
    fuel: name of the fuel (one of the possible isomers)
    given the name of a fuel in a spc_dct, it finds all species
    involved in its combustion mech by stoichiometry:
    - fuel isomers
    - fuel radicals
    - main radical additions to fuel: fuel+H/OH/O2/O/HO2, R+O
    - low T mech: RO2, RO4, RO3-H

    returns
    - list of species of the corresponding stoichiometries
    - species dataframe: for each species, identifies its subset: 'fuel',
        'fuel_rad', 'fuel_add_H', ...
    """
    species_list = []
    species_subset_df = pd.Series()
    # extract formulas
    fml_df = extract_fml_df(spc_dct)

    # generate a list of stoichiometries to extract and the corresponding labels
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO']].values

    stoich_dct = {
        'FUEL': stoich_fuel,
        'FUEL_RAD': stoich_fuel+[0, -1, 0],
        'FUEL_ADD_H': stoich_fuel+[0, 1, 0],
        'FUEL_ADD_O': stoich_fuel+[0, 0, 1],
        'FUEL_ADD_OH': stoich_fuel+[0, 1, 1],
        'FUEL_ADD_O2': stoich_fuel+[0, 0, 2],
        'R_O': stoich_fuel+[0, -1, 1],
        'R_O2': stoich_fuel+[0, -1, 2],
        'R_O4': stoich_fuel+[0, -1, 4],
        'R_O3-H': stoich_fuel+[0, -2, 3]}

    # extract desired species and assign labels

    for stoich_type, stoich in stoich_dct.items():
        species = extract_species(stoich, fml_df)
        series = pd.Series(stoich_type, index=species)
        species_subset_df = species_subset_df.append(series)
        species_list = species_list + species

    return species_list, species_subset_df


def extract_fml_df(spc_dct):
    """
    given species dictionary, returns
    fml_df: dataframe with index = species names; columns 'fml' (stoichiometry),
    'nC','nH','nO' (number of C/H/O atmoms in the formula)
    """
    fml_df = pd.DataFrame(index=list(spc_dct.keys()),
                          columns=['fml', 'nC', 'nH', 'nO'])
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich = spc_dct[key]['inchi']
            fml_dct = automol.inchi.formula(ich)
            fml_df['fml'][key] = automol.formula.string2(fml_dct)
            fml_df['nC'][key] = automol.formula.element_count(fml_dct, 'C')
            fml_df['nH'][key] = automol.formula.element_count(fml_dct, 'H')
            fml_df['nO'][key] = automol.formula.element_count(fml_dct, 'O')

    return fml_df


def extract_species(n_cho, fml_df):
    """
    extracts species corresponding to a given stoichiometry n_CHO from the formulas dataframe
    n_CHO: numpy array [x,y,z] where x = n of C atoms, y = n of H atoms, z = number of O atoms
    fml_df: dataframe with index = species names; columns 'fml' (stoichiometry),
    'nC','nH','nO' (number of C/H/O atmoms in the formula)
    returns: list of species with the corresponding n_CHO
    """
    n_carbons = n_cho[0]
    n_hydrogens = n_cho[1]
    n_oxygens = n_cho[2]
    species_set = list(fml_df[(fml_df['nC'] == n_carbons) & (
        fml_df['nH'] == n_hydrogens) & (fml_df['nO'] == n_oxygens)].index)

    return species_set


# list of formulas for products/reactants identified for the sublcasses
fmls_set = np.array(['H1', 'O1', 'H1O1', 'O2', 'H1O2', 'C1H3'])


def classify_unimol(rcts, prds, spc_dct):
    """
    Classifies unimolecular reaction from reactants and products names (tuples)
    - A=B: isomerization
    - A=C+H/O/OH/O2/HO2/CH3: addition-H/O/OH/O2/HO2/CH3  /
        recombination (depends on multiplicity of the products)
        it would be nice to distinguish the type of bond that is being broken
    - A=C+D: decomposition
    - A=C+D+E.. : decomposition(lumped)
    """
    # extract formula dictionary
    fml_df = extract_fml_df(spc_dct)
    mult_rcts = util.get_mult(rcts, spc_dct)
    mult_prds = util.get_mult(prds, spc_dct)

    if len(prds) == 1:
        rxn_class_broad = 'Isomerization'
    elif len(prds) == 2:
        rxn_class_broad = 'Decomposition'
        # derive products multiplicity and formulas
        if mult_rcts == 1 and mult_prds >= 4:
            rxn_class_broad = 'Bond fission'
        elif mult_rcts > 1 and mult_prds == 2:
            rxn_class_broad = 'Beta-scission'

        prds_fmls = np.array(
            [fml_df['fml'][prds[0]], fml_df['fml'][prds[1]]])

        # check product composition
        # give priority to the second product (should be the lightest)
        if any(prds_fmls[1] == fmls_set):
            rxn_class_broad += ' +{}'.format(prds[1])
        elif any(prds_fmls[0] == fmls_set):
            rxn_class_broad += ' +{}'.format(prds[0])

    elif len(prds) > 2:
        rxn_class_broad = 'Decomposition(lumped)'

    return rxn_class_broad


def classify_bimol(rcts, prds, spc_dct):
    """
    Classifies bimolecular reactions from reactants and products names (tuples)
      with 1 product
    - A+H/O/OH/O2/HO2/CH3 = B : addition-H/O/OH/O2/HO2/CH3  /
        recombination (depends on the multiplicity of the reactants)
      with 2 products
    - A+R=B+RH: Habstraction-R (subclass indicates the abstractor)
    - A+R=A+R: isomerization-bim (isomerization aided by a radical.
        ex. CH2+H=CH2(S)+H, C6H6+H=FULV+H)
    - A+B=C+D: addition-decomposition - branch/prop/term
    - A+B=C+D+E..: addition-decomposition(lumped) - branch/prop/term
    """
    # extract formula dictionary
    fml_df = extract_fml_df(spc_dct)
    # extracts reactants and products multiplicity
    mult_rcts = util.get_mult(rcts, spc_dct)
    mult_prds = util.get_mult(prds, spc_dct)

    if mult_rcts < 4:
        rxn_class_broad = 'Addition'
    elif mult_rcts >= 4:
        rxn_class_broad = 'Recombination'

    if len(prds) == 1:

        rcts_fmls = np.array(
            [fml_df['fml'][rcts[0]], fml_df['fml'][rcts[1]]])

        # check reactant composition
        # give priority to the second reactant (should be the lightest)
        if any(rcts_fmls[1] == fmls_set):
            rxn_class_broad += ' {}'.format(rcts[1])
        elif any(rcts_fmls[0] == fmls_set):
            rxn_class_broad += ' {}'.format(rcts[0])

    elif len(prds) == 2:

        rxn_class_broad += '-decomposition'

        # H abstraction
        flag_habs = classify_habs(rcts, prds, fml_df, spc_dct)
        # Bimolecular isomerization
        flag_isom_bim = classify_isom_bim(rcts, prds, fml_df)

        if flag_habs == 1:
            rxn_class_broad = 'H abstraction'
        elif flag_isom_bim == 1:
            rxn_class_broad = 'Bimol Isomerization'

    elif len(prds) > 2:
        rxn_class_broad += '-decomposition(lumped)'

    # add subclass related to branching/propagation/termination for add-deco rxns
    if 'decomposition' in rxn_class_broad:
        # check if branching/propagation/termination
        rxn_class_broad += bran_prop_term(mult_rcts, mult_prds)

    return rxn_class_broad


def classify_isom_bim(rcts, prds, fml_df):
    '''
    Check if an A+B=C+D reaction is a bimolecular isomerization of the kind
    A+R=B+R
    By checking if the reactants and the products both match
        - 1 couple of rct/prd must be the exact same species
        - the other couple of rct/prd must match the stoichiometry
    returns 1 if true
    returns 0 if false
    '''
    rcts = np.array(rcts)
    prds = np.array(prds)
    rcts_fmls = np.array(
        [fml_df['fml'][rcts[0]], fml_df['fml'][rcts[1]]])
    prds_fmls = np.array(
        [fml_df['fml'][prds[0]], fml_df['fml'][prds[1]]])

    flag_species = (any(rcts[0] == prds) or any(rcts[1] == prds))
    flag_stoich = (any(rcts_fmls[0] == prds_fmls)
                   and any(rcts_fmls[1] == prds_fmls))

    if (flag_species == 1 and flag_stoich == 1):
        flag_isom_bim = 1
    else:
        flag_isom_bim = 0

    return flag_isom_bim


def classify_habs(rcts, prds, fml_df, spc_dct):
    '''
    Check if an A+B=C+D reaction is an Habstraction based on the stoichiometries
    and multiplicities of reactants and products
    Returns 0 if it is not an Habs, returns 1 if it s
    '''
    if not isinstance(rcts, tuple) or not isinstance(prds, tuple):
        print('error: reactants and products are not tuples')
        sys.exit()
    elif len(rcts) != 2 or len(prds) != 2:
        print('error: not A+B=C+D reaction')
        sys.exit()

    mult_rct = np.array([util.get_mult(rcts[0], spc_dct),
                         util.get_mult(rcts[1], spc_dct)])
    mult_prd = np.array([util.get_mult(prds[0], spc_dct),
                         util.get_mult(prds[1], spc_dct)])

    if (any(mult_rct == 1) and any(mult_rct > 1) and
            any(mult_prd == 1) and any(mult_prd > 1)):

        species_rct = rcts[np.where(mult_rct == 1)[0][0]]
        species_prd = prds[np.where(mult_prd > 1)[0][0]]
        stoich_add = [0, -1, 0]
        flag_try_prd2 = 0

    # Habs with O2 and O: multiplicities are 1*3=2*2
    elif (any(mult_rct == 1) and any(mult_rct == 3) and
          all(mult_prd == 2)):

        species_rct = rcts[np.where(mult_rct == 3)[0][0]]
        species_prd = prds[0]
        stoich_add = [0, +1, 0]
        # try the second product in case the first is not right
        flag_try_prd2 = 1

    try:
        flag_habs = set_flag_habs(species_rct, species_prd, fml_df, stoich_add)

        if flag_try_prd2 == 1 and flag_habs == 0:
            # try the second product
            species_prd = prds[1]
            flag_habs = set_flag_habs(
                species_rct, species_prd, fml_df, stoich_add)

    except NameError:
        flag_habs = 0

    return flag_habs


def set_flag_habs(species_rct, species_prd, fml_df, stoich_add):
    '''
    given reactant and product single species
    the corresponding formulas
    and the stoichiometry to derive product target
    checks if the rct+stoich_add corresponds to the product one
    returns flag 0 or 1 integer
    '''
    stoich_rct = fml_df.loc[species_rct][['nC', 'nH', 'nO']].values
    stoich_prd = fml_df.loc[species_prd][['nC', 'nH', 'nO']].values
    stoich_prd_target = stoich_rct+stoich_add

    flag_habs = int(all(stoich_prd == stoich_prd_target))

    return flag_habs


def bran_prop_term(mult_rcts, mult_prds):
    '''
    Given reactants and products multiplicity:
    Checks if propagation, termination or branching
    '''

    if mult_rcts == mult_prds:
        add = ' - propagation'
    elif mult_rcts > mult_prds:
        add = ' - termination'
    elif mult_rcts < mult_prds:
        add = ' - branching'
    else:
        add = ''

    return add
