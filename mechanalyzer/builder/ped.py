""" build the ped of prompt products
"""

from mechanalyzer.calculator import statmodels


def ped_frag1(ped_df, hotfrg, otherfrg, modeltype,
              dos_df=None, dof_info=None):
    """ call ped_models class in statmodels and compute P(E1)

        :param ped_df: dataframe(columns:P, rows:T) energy distrib. series
        :type ped_df: dataframe(series(float))
        :param hotfrg: selected hot fragment between frag1 and frag2
        :type hotfrg: str
        :param otherfrg: the other fragment
        :type otherfrg: str
        :param dos_df: rovibr dos for each fragment
        :type dos_df: dataframe(index=energy, columns=[frag1, frag2])
        :param dof_info: vib-rot degrees of freedom and molecular weight
        :type: dataframe(index=species, columns=['vib dof', 'rot dof', 'mw'])
        :param ene_bw: backward energy barrier TS-PRODS
        :type ene_bw: float
        :param modeltype: type of model to be implemented
        :type modeltype: str

        :return P_E1_prod1: energy distribution of the product prod
        :rtype: dataframe(series(float, index=energy), index=T, columns=P)
    """

    # call class
    ped_prod1_fct = statmodels.PEDModels(
        ped_df, hotfrg, otherfrg,
        dos_df=dos_df, dof_info=dof_info)

    ped_df_frag1 = ped_prod1_fct.compute_ped(modeltype)

    return ped_df_frag1
