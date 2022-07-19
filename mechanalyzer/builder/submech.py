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
from mechanalyzer.calculator import formulas

# atoms order: C, H, O, N, S, CL
STOICH_DEFAULT = [0, 2, 2, 0, 0, 0]
STOICH_DEFAULT_DEL = [1, 0, 2, 0, 0, 0] # DELETE ALL ABOVE THESE STOICH
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
    fml_df_reduced = formulas.extract_fml_df(spc_dct_reduced)

    fml_df = formulas.extract_fml_df(spc_dct)
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO', 'nN', 'nS', 'nCl']].values
    stoich = stoich_fuel + stoich_limit

    # extract desired species and assign labels
    species = formulas.extract_species_sub(stoich, fml_df_reduced)
    series = pd.Series('SUBFUEL', index=species)
    species_subset_df = species_subset_df.append(series)
    species_list += species

    return species_list, species_subset_df

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
    fml_df_reduced = formulas.extract_fml_df(spc_dct_reduced)

    fml_df = formulas.extract_fml_df(spc_dct)
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO', 'nN', 'nS', 'nCl']].values
    stoich = stoich_fuel + stoich_limit

    # extract desired species and assign labels
    species = formulas.extract_species_sub(stoich, fml_df_reduced)
    series = pd.Series('SUBFUEL', index=species)
    species_subset_df = species_subset_df.append(series)
    species_list += species

    return species_list, species_subset_df

def species_subset_del(fuel, spc_dct, stoich_limit=STOICH_DEFAULT_DEL):
    """ call species_subset but and also extract all species
        below stoich: stoich_fuel+stoich_limit
    """
    species_subset_df = pd.Series()
    # add all stoichiometries below that of interest
    # filter extracted species from spc_dct first
    spc_dct_reduced = copy.deepcopy(spc_dct)
    fml_df = formulas.extract_fml_df(spc_dct_reduced)
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO', 'nN', 'nS', 'nCl']].values
    stoich = stoich_fuel + stoich_limit

    # extract desired species and assign labels
    species_del = formulas.extract_species_above(stoich, fml_df)
    [spc_dct_reduced.pop(sp) for sp in species_del]
    species_list = list(spc_dct_reduced.keys())
    series = pd.Series('SUBFUEL', index=species_list)
    species_subset_df = species_subset_df.append(series)

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
    fml_df = formulas.extract_fml_df(spc_dct)

    # Generate list of stoichiometries to extract and the corresponding labels
    stoich_fuel = fml_df.loc[fuel][['nC', 'nH', 'nO', 'nN', 'nS', 'nCl']].values

    # extract desired species and assign labels

    for stoich_type, stoich in STOICH_DCT_ADD.items():
        species = formulas.extract_species(stoich+stoich_fuel, fml_df)
        series = pd.Series(stoich_type, index=species)
        species_subset_df = species_subset_df.append(series)
        species_list += species

    return species_list, species_subset_df

