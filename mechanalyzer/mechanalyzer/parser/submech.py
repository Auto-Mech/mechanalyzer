"""
Functions for
- submechanism identification
- broad rxn classification
"""

import os
import sys
import mechanalyzer
import chemkin_io
import automol
from mechanalyzer.parser import ckin_ as ckin
import pandas as pd
import numpy as np


def species_subset(fuel, spc_dct):
    """
    fuel: name of the fuel (one of the possible isomers)
    given the name of a fuel in a spc_dct, it finds all species involved in its combustion mech by stoichiometry:
    - fuel isomers
    - fuel radicals
    - main radical additions to fuel: fuel+H/OH/O2/O/HO2, R+O
    - low T mech: RO2, RO4, RO3-H

    returns
    - list of species of the corresponding stoichiometries
    - species dataframe: for each species, identifies its subset: 'fuel', 'fuel_rad', 'fuel_add_H', ...
    """
    species_list = []
    species_subset_df = pd.Series()
    # extract formulas
    fml_df = extract_fml_df(spc_dct)

    # generate a list of stoichiometries to extract and the corresponding labels
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO']].values

    stoich_dct = {
        'fuel': stoich_fuel,
        'fuel_rad': stoich_fuel+[0, -1, 0],
        'fuel_add_H': stoich_fuel+[0, 1, 0],
        'fuel_add_O': stoich_fuel+[0, 0, 1],
        'fuel_add_OH': stoich_fuel+[0, 1, 1],
        'fuel_add_O2': stoich_fuel+[0, 0, 2],
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

    print(species_subset_df)
    print(species_list)

    return species_list, species_subset_df


def extract_fml_df(spc_dct):
    """
    given species dictionary, returns
    fml_df: dataframe with index = species names; columns 'fml' (stoichiometry), 'nC','nH','nO' (number of C/H/O atmoms in the formula)
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


def extract_species(n_CHO, fml_df):
    """
    extracts species corresponding to a given stoichiometry n_CHO from the formulas dataframe
    n_CHO: numpy array [x,y,z] where x = n of C atoms, y = n of H atoms, z = number of O atoms
    fml_df: dataframe with index = species names; columns 'fml' (stoichiometry), 'nC','nH','nO' (number of C/H/O atmoms in the formula)
    returns: list of species with the corresponding n_CHO
    """
    nC = n_CHO[0]
    nH = n_CHO[1]
    nO = n_CHO[2]
    species_subset = list(fml_df[(fml_df['nC'] == nC) & (
        fml_df['nH'] == nH) & (fml_df['nO'] == nO)].index)

    return species_subset


def classify_unimol(rct_names, prd_names, spc_dct):
    """
    Classifies unimolecular reaction from reactants and products names (tuples)
    - A=B: isomerization
    - A=C+H/O/OH/O2/HO2/CH3: addition-H/O/OH/O2/HO2/CH3  / recombination (depends on multiplicity of the products)
    - A=C+D: decomposition
    - A=C+D+E.. : decomposition(lumped)
    """
    if len(prd_names) == 1:
        rxn_class_broad = 'isomerization'
    elif len(prd_names) == 2:
        rxn_class_broad = 'decomposition'
        # guarda il reagente 2
    elif len(prd_names) > 2:
        rxn_class_broad = 'decomposition(lumped)'

    return rxn_class_broad


def classify_bimol(rct_names, prd_names, spc_dct):
    """
    Classifies bimolecular reactions from reactants and products names (tuples)
      with 1 product
    - A+H/O/OH/O2/HO2/CH3 = B : addition-H/O/OH/O2/HO2/CH3  / recombination (depends on the multiplicity of the reactants)
      with 2 products
    - A+R=B+RH: Habstraction-R (subclass indicates the abstractor)
    - A+R=A+R: isomerization-bim (isomerization aided by a radical. ex. CH2+H=CH2(S)+H, C6H6+H=FULV+H)
    - A+B=C+D: addition-decomposition
    - A+B=C+D+E..: addition-decomposition(lumped)
    """
    if len(prd_names) == 1:
        mult = spc_dct[rct_names[0]]['mult']*spc_dct[rct_names[1]]['mult']
        if mult == 1:
            rxn_class_broad = 'addition'
        elif mult > 1:
            rxn_class_broad = 'recombination'
    elif len(prd_names) == 2:
        rxn_class_broad = 'bimol-unclassified'
    elif len(prd_names) > 2:
        rxn_class_broad = 'addition-decomposition(lumped)'

    return rxn_class_broad
