""" extract formulas from spc dct and put in dataframes
"""

import pandas as pd
import automol

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
            fml_dct = automol.chi.formula(ich)
            fml_df.loc[key, 'fml'] = automol.form.string2(fml_dct)
            for X in ['C', 'H', 'O', 'N', 'S', 'Cl']:
                fml_df.at[key, 'n{}'.format(X)] = automol.form.element_count(fml_dct, X)

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

def extract_species_core(n_at, fml_df):
    """ extract any fml with less C/N/S/Cl atoms than indicated. exclude N containing species
    """
    species_set = list(fml_df[(fml_df['nC'] <= n_at[0] ) & (
        fml_df['nN'] <= n_at[3]) & (fml_df['nS'] <= n_at[4]) & (
        fml_df['nCl'] <= n_at[5])].index)
    
    return species_set

def extract_species_above(n_at, fml_df):
    """ extract any fml above the indicated limits
    """
    species_set = list(fml_df[(fml_df['nC'] >= n_at[0] ) & (
         (fml_df['nO'] >= n_at[2]))].index)
    # also delete anything with N that is not N2
    species_set.extend(list(fml_df[(fml_df['nN'] >= 1 )].index))
    species_set.remove('N2')
    
    return species_set