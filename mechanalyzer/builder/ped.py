""" build the ped of prompt products
"""
from mechanalyzer.calculator import statmodels


def ped_frag1(ped_df, hotfrg, otherfrg, modeltype_list, dos_df=None, dof_info=None, E_BW=None):
    """ call ped_models class in statmodels and compute P(E1)

        :param ped_df: dataframe(columns:P, rows:T) with the Series of energy distrib
        :type ped_df: dataframe(series(float))
        :param hotfrg: selected hot fragment between frag1 and frag2
        :type hotfrg: str
        :param otherfrg: the other fragment
        :type otherfrg: str
        :param modeltype_list: models to be used for P(E1) calculations
        :type modeltype_list: list(str)
        :param dos_df: rovibr dos for each fragment
        :type dos_df: dataframe(index=energy, columns=[frag1, frag2])
        :param dof_info: dataframe with vibrat/rot degrees of freedom and molecular weight
        :type: dataframe(index=species, columns=['vib dof', 'rot dof', 'mw'])
        :param E_BW: backward energy barrier TS-PRODS
        :type E_BW: float
        :param modeltype: type of model to be implemented
        :type modeltype: str

        :return P_E1_prod1: energy distribution of the product prod
        :rtype: dct{modelype: dataframe(series(float, index=energy), index=T, columns=P)}
    """

    # call class
    ped_prod1_fct = statmodels.ped_models(
        ped_df, hotfrg, otherfrg, dos_df=dos_df, dof_info=dof_info, E_BW=E_BW)

    ped_df_frag1_dct = dict.fromkeys(modeltype_list)

    for modeltype in modeltype_list:
        ped_df_frag1_dct[modeltype] = ped_prod1_fct.compute_ped(modeltype)

    return ped_df_frag1_dct
