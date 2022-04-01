"""
Functions for
- submechanism identification
- broad rxn classification
"""

import sys
import copy
import numpy
import pandas as pd
import automol
from mechanalyzer.parser._util import get_mult

# list of formulas for products/reactants identified for the sublcasses
FMLS_SET = numpy.array(['H1', 'O1', 'H1O1', 'O2', 'H1O2', 'C1H3'])
# atoms order: C, H, O, N, S, CL
STOICH_DEFAULT = [0, 2, 2, 0, 0, 0]
STOICH_DCT_ADD = {
    'FUEL': [0, 0, 0, 0, 0, 0],
    'FUEL_RAD': [0, -1, 0, 0, 0, 0],
    'FUEL_ADD_H': [0, 1, 0, 0, 0, 0],
    'FUEL_ADD_CH3': [1, 3, 0, 0, 0, 0],
    'FUEL_ADD_O': [0, 0, 1, 0, 0, 0],
    'FUEL_ADD_OH': [0, 1, 1, 0, 0, 0],
    'FUEL_ADD_O2': [0, 0, 2, 0, 0, 0],
    'R_CH3': [1, 2, 0, 0, 0, 0],
    'R_O': [0, -1, 1, 0, 0, 0],
    'R_O2': [0, -1, 2, 0, 0, 0],
    'R_O4': [0, -1, 4, 0, 0, 0],
    'R_O3-H': [0, -2, 3, 0, 0, 0]}

def species_subset_ext(fuel, spc_dct, stoich_limit=STOICH_DEFAULT):
    """ call species_subset but and also extract all species
        below stoich: stoich_fuel+stoich_limit
    """
    # extract species subset
    species_list, species_subset_df = species_subset(fuel, spc_dct)
    # add all stoichiometries below that of interest
    # filter extracted species from spc_dct first
    spc_dct_reduced = copy.deepcopy(spc_dct)
    [spc_dct_reduced.pop(sp) for sp in species_list]
    fml_df_reduced = extract_fml_df(spc_dct_reduced)

    fml_df = extract_fml_df(spc_dct)
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO', 'nN', 'nS', 'nCl']].values
    stoich = stoich_fuel + stoich_limit

    # extract desired species and assign labels
    species = extract_species_sub(stoich, fml_df_reduced)
    series = pd.Series('SUBFUEL', index=species)
    species_subset_df = species_subset_df.append(series)
    species_list += species

    return species_list, species_subset_df

def species_subset(fuel, spc_dct):
    """ Given the name of a fuel in a spc_dct, it finds all species
        involved in its combustion mech by stoichiometry

        - fuel isomers
        - fuel radicals
        - main radical additions to fuel: fuel+H/CH3/OH/O2/O/HO2, R+O
        - low T mech: RO2, RO4, RO3-H

        Returns all species with corresponding stoichiometries with
        subset identifiers: 'fuel', 'fuel_rad', 'fuel_add_H', ...

        :param fuel: name of the fuel (one of the possible isomers)
        :type fuel: str
        :param spc_dct: species dict
        :type spc_dct: dict[]
        :rtype: pandas.dataframe
    """

    species_list = []
    species_subset_df = pd.Series()

    # extract formulas
    fml_df = extract_fml_df(spc_dct)

    # Generate list of stoichiometries to extract and the corresponding labels
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO', 'nN', 'nS', 'nCl']].values

    # extract desired species and assign labels

    for stoich_type, stoich in STOICH_DCT_ADD.items():
        species = extract_species(stoich+stoich_fuel, fml_df)
        series = pd.Series(stoich_type, index=species)
        species_subset_df = species_subset_df.append(series)
        species_list += species

    return species_list, species_subset_df


def extract_fml_df(spc_dct):
    """ Given species dictionary, builds a formula Pandas dataframe:
            index = species names; columns 'fml' (stoichiometry)
            'nC','nH','nO' (number of C/H/O atmoms in the formula)

        :param spc_dct:
        :type spc_dct: dict[]
    """

    fml_df = pd.DataFrame(index=list(spc_dct.keys()),
                          columns=['fml', 'nC', 'nH', 'nO', 'nN', 'nS', 'nCl'])
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich = spc_dct[key]['inchi']
            fml_dct = automol.inchi.formula(ich)
            fml_df['fml'][key] = automol.formula.string2(fml_dct)
            fml_df['nC'][key] = automol.formula.element_count(fml_dct, 'C')
            fml_df['nH'][key] = automol.formula.element_count(fml_dct, 'H')
            fml_df['nO'][key] = automol.formula.element_count(fml_dct, 'O')
            fml_df['nN'][key] = automol.formula.element_count(fml_dct, 'N')
            fml_df['nS'][key] = automol.formula.element_count(fml_dct, 'S')
            fml_df['nCl'][key] = automol.formula.element_count(fml_dct, 'Cl')

    return fml_df


def extract_species(n_at, fml_df):
    """ Extracts species corresponding to a given stoichiometry n_CHONSCl
        from the formulas dataframe.

        :param n_cho: numpy array [x,y,z,...] where
            x = n of C atoms, y = n of H atoms, ...
        :type: n_cho: numpy.ndarray
        :param fml_df: dataframe with
            index = species names; columns 'fml' (stoichiometry),
            'nC','nH','nO' (number of C/H/O atmoms in the formula)
        :returns: list of species with the corresponding n_CHONSCl
    """

    species_set = list(fml_df[(fml_df['nC'] == n_at[0]) & (
        fml_df['nH'] == n_at[1]) & (fml_df['nO'] == n_at[2]) & (
        fml_df['nN'] == n_at[3]) & (fml_df['nS'] == n_at[4]) & (
        fml_df['nCl'] == n_at[5])].index)

    return species_set

def extract_species_sub(n_at, fml_df):
    """ as above but leq
    """

    species_set = list(fml_df[(fml_df['nC'] <= n_at[0]) & (
        fml_df['nH'] <= n_at[1]) & (fml_df['nO'] <= n_at[2]) & (
        fml_df['nN'] <= n_at[3]) & (fml_df['nS'] <= n_at[4]) & (
        fml_df['nCl'] <= n_at[5])].index)

    return species_set

def classify_unimol(rcts, prds, spc_dct):
    """ Classifies unimolecular reaction from reactants and products names:
        - A=B: isomerization
        - A=C+H/O/OH/O2/HO2/CH3: addition
        - A=C+D: decomposition
        - A=C+D+E.. : decomposition(lumped)

        For recombination (depends on multiplicity of the products)
        it would be nice to distinguish type of bond that being broken

        :param rcts: reactant names
        :type rcts: tuple
        :param prds: product names
        :type prds: tuple
        :param spc_dct:
        :type spc_dct: dict[]
    """

    # extract formula dictionary
    fml_df = extract_fml_df(spc_dct)
    mult_rcts = get_mult(rcts, spc_dct)
    mult_prds = get_mult(prds, spc_dct)

    if len(prds) == 1:
        rxn_class_broad = 'Isomerization'
    elif len(prds) == 2:
        rxn_class_broad = 'Decomposition'
        # derive products multiplicity and formulas
        if mult_rcts == 1 and mult_prds >= 4:
            rxn_class_broad = 'Bond fission'
        elif mult_rcts > 1 and mult_prds == 2:
            rxn_class_broad = 'Beta-scission'

        prds_fmls = numpy.array(
            [fml_df['fml'][prds[0]], fml_df['fml'][prds[1]]])

        # check product composition
        # give priority to the second product (should be the lightest)
        if any(prds_fmls[1] == FMLS_SET):
            rxn_class_broad += f' +{prds[1]}'
        elif any(prds_fmls[0] == FMLS_SET):
            rxn_class_broad += f' +{prds[0]}'

    elif len(prds) > 2:
        rxn_class_broad = 'Decomposition(lumped)'

    return rxn_class_broad


def classify_bimol(rcts, prds, spc_dct):
    """ Classifies bimolecular reactions from reactants and products names
        with 1 product:

        - A+H/O/OH/O2/HO2/CH3 = B
        - A+R=B+RH: Habstraction-R (subclass indicates the abstractor)
        - A+R=A+R: isomerization-bim (isomerization aided by a radical.
          ex. CH2+H=CH2(S)+H, C6H6+H=FULV+H)
        - A+B=C+D: addition-decomposition - branch/prop/term
        - A+B=C+D+E..: addition-decomposition(lumped) - branch/prop/term

        For recombination (depends on the multiplicity of the reactants)
        with 2 products.

        :param rcts: reactant names
        :type rcts: tuple
        :param prds: product names
        :type prds: tuple
        :param spc_dct:
        :type spc_dct: dict[]
    """

    # extract formula dictionary
    fml_df = extract_fml_df(spc_dct)

    # extracts reactants and products multiplicity
    mult_rcts = get_mult(rcts, spc_dct)
    mult_prds = get_mult(prds, spc_dct)

    if mult_rcts < 4:
        rxn_class_broad = 'Addition'
    elif mult_rcts >= 4:
        rxn_class_broad = 'Recombination'

    if len(prds) == 1:

        rcts_fmls = numpy.array(
            [fml_df['fml'][rcts[0]], fml_df['fml'][rcts[1]]])

        # check reactant composition
        # give priority to the second reactant (should be the lightest)
        if any(rcts_fmls[1] == FMLS_SET):
            rxn_class_broad += f' {rcts[1]}'
        elif any(rcts_fmls[0] == FMLS_SET):
            rxn_class_broad += f' {rcts[0]}'

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

    # add subclass related to branching/propagation/termination
    # for addition-decomposition rxns
    if 'decomposition' in rxn_class_broad:
        # check if branching/propagation/termination
        rxn_class_broad += bran_prop_term(mult_rcts, mult_prds)

    return rxn_class_broad


def classify_isom_bim(rcts, prds, fml_df):
    """ Check if an A+B=C+D reaction is a bimolecular isomerization
        of the kind A+R=B+R by checking

        if the reactants and the products both match
        - 1 couple of rct/prd must be the exact same species
        - the other couple of rct/prd must match the stoichiometry

        :param rcts: reactant names
        :type rcts: tuple
        :param prds: product names
        :type prds: tuple
        :param fml_df:
        :type fml_df:
        :rtype: bool
    """
    rcts = numpy.array(rcts)
    prds = numpy.array(prds)
    rcts_fmls = numpy.array(
        [fml_df['fml'][rcts[0]], fml_df['fml'][rcts[1]]])
    prds_fmls = numpy.array(
        [fml_df['fml'][prds[0]], fml_df['fml'][prds[1]]])

    flag_species = (any(rcts[0] == prds) or any(rcts[1] == prds))
    flag_stoich = (any(rcts_fmls[0] == prds_fmls)
                   and any(rcts_fmls[1] == prds_fmls))

    return flag_species and flag_stoich


def classify_habs(rcts, prds, fml_df, spc_dct):
    """ Check if an A+B=C+D reaction is an hydrogen abstraction
        based on the stoichiometries and multiplicities of
        reactants and products.

        :rtype: bool
    """

    if not isinstance(rcts, tuple) or not isinstance(prds, tuple):
        print('error: reactants and products are not tuples')
        sys.exit()
    elif len(rcts) != 2 or len(prds) != 2:
        print('error: not A+B=C+D reaction')
        sys.exit()

    mult_rct = numpy.array([get_mult(rcts[0], spc_dct),
                            get_mult(rcts[1], spc_dct)])
    mult_prd = numpy.array([get_mult(prds[0], spc_dct),
                            get_mult(prds[1], spc_dct)])

    if (any(mult_rct == 1) and any(mult_rct > 1) and
            any(mult_prd == 1) and any(mult_prd > 1)):

        species_rct = rcts[numpy.where(mult_rct == 1)[0][0]]
        species_prd = prds[numpy.where(mult_prd > 1)[0][0]]
        stoich_add = [0, -1, 0]
        flag_try_prd2 = False

    # Habs with O2 and O: multiplicities are 1*3=2*2
    elif (any(mult_rct == 1) and any(mult_rct == 3) and
          all(mult_prd == 2)):

        species_rct = rcts[numpy.where(mult_rct == 3)[0][0]]
        species_prd = prds[0]
        stoich_add = [0, +1, 0]
        # try the second product in case the first is not right
        flag_try_prd2 = True

    try:
        flag_habs = set_flag_habs(species_rct, species_prd, fml_df, stoich_add)

        if flag_try_prd2 and not flag_habs:
            # try the second product
            species_prd = prds[1]
            flag_habs = set_flag_habs(
                species_rct, species_prd, fml_df, stoich_add)
    except NameError:
        flag_habs = False

    return flag_habs


def set_flag_habs(spc_rct, spc_prd, fml_df, stoich_add):
    """ Uses the formulae of the reactant and product to derive
        the product target and checks if the rct+stoich_add corresponds
        to the product one.

        :param spc_rct:
        :type spc_rct:
        :param spc_prd:
        :type spc_prd:
        :param fml_df:
        :type fml_df:
        :param stoich_add:
        :type: stoich_add:
    """

    stoich_rct = fml_df.loc[spc_rct][['nC', 'nH', 'nO']].values
    stoich_prd = fml_df.loc[spc_prd][['nC', 'nH', 'nO']].values
    stoich_prd_target = stoich_rct+stoich_add

    return all(stoich_prd == stoich_prd_target)


def bran_prop_term(rct_muls, prd_muls):
    """ Checks if a reaction can be further classified as a
        propagation, termination, or branching reaction using the
        reaction multiplicities.

        :param rct_muls: reactant multiplicities
        :type rct_muls: tuple(int)
        :param prd_muls: product multiplicities
        :type prd_muls: tuple(int)
    """

    if rct_muls == prd_muls:
        add = ' - propagation'
    elif rct_muls > prd_muls:
        add = ' - termination'
    elif rct_muls < prd_muls:
        add = ' - branching'
    else:
        add = ''

    return add
