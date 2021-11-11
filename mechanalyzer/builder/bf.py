""" build the dictionaries of the product branching fractions
"""

from mechanalyzer.calculator import bf
# probably more useful to turn this into a class


def bf_tp_dct(modeltype_list, ped_dct, hoten_df, bf_threshold, savefile=False):
    """ Build a branching fractions dictionary as a function of temeprature and pressure
        containing the BFs of each product of the PES
        :param modeltype_list: models used for P(E1) calculations
        :type modeltype_list: list(str)
        :param prod_name: name of the product in pedspecies and in hotenergies
        :type prod_name: list(str)
        :param ped_dct: dictionary of modeltypes with corresponding dataframe[P][T] with the Series of energy distrib [en: prob(en)]
        :type ped_dct: dct{str: dataframe(series(float))}
        :param hoten_df: hot branching fractions for hotspecies
        :type hoten_df: df[P][T]:df[allspecies][energies]}
        :return bf_tp_dct: branching fractions at T,P for each product for the selected hotspecies
        :rtype: dct{model: dct{species: {pressure: (array(T), array(BF))}}}
    """
    bf_tp_dct = dict.fromkeys(modeltype_list)
    for modeltype in modeltype_list:
        ped_df = ped_dct[modeltype]
        bf_tp_df = bf.bf_tp_df_full(ped_df, hoten_df)
        bf_tp_dct_species = bf.bf_tp_df_todct(
            bf_tp_df, bf_threshold, model=modeltype, savefile=savefile)

        bf_tp_dct[modeltype] = bf_tp_dct_species

    return bf_tp_dct


def merge_bf_ktp(bf_ktp_dct, ktp_dct, frag_reacs, frag1, frag2, hotsp_dct):
    """ derive k' = bf*k and rename final ktp dictionary appropriately
        :param bf_tp_dct: branching fractions at T,P for each product for the selected hotspecies
        :type: dct{model: dct{species: {pressure: (array(T), array(BF))}}}
        :param ktp_dct: rates of the original reaction A=>B to split in A=>Pi with Pi = species in bf_ktp_dct
        (B decomposes to Pi)
        :type ktp_dct: dictionary {P: (T, k)}
        :param frag_reacs: reactant fragments ['A','B']
        :type frag_reacs: list(str)
        :param frag1: hot dissociating product 
        :type frag1: str
        :param frag2: remaining product
        :type frag2: str
        :param hotsp_dct: dictionary of hotspecies and corresponding fragments (if bimol)
        :type hotsp_dct: {species_unimol: [species_unimol], species_bimol: [frag1, frag2], ...}
        :return rxn_ktp_dct: ktp dct of final rate constants for different channels
        :rtype: {model: {rxn: {P: (T, k)}}}
    """
    modeltype_list = list(bf_ktp_dct.keys())
    rxn_ktp_dct = dict.fromkeys(modeltype_list)
    for modeltype in modeltype_list:
        ktp_dct_model_i = bf.merge_bf_rates(bf_ktp_dct[modeltype], ktp_dct)
        ktp_dct_model_i_new = rename_ktp_dct(
            ktp_dct_model_i, frag_reacs, frag1, frag2, hotsp_dct)
        rxn_ktp_dct[modeltype] = ktp_dct_model_i_new

    return rxn_ktp_dct


def rename_ktp_dct(ktp_dct, frag_reacs, hotprod, otherprod, hotsp_dct):
    """ rename ktp dictionary with appropriate names for prompt dissociation
        ktp_dct.keys(): sp
        renamed_ktp_dct.keys(): rctname=>sp
        if sp is the original product, the reaction is reversible =
        :param rxn_ktp_dct: ktp dct of final rate constants for different channels
        :type rxn_ktp_dct: {sp: {P: (T, k)}}
        :param frag_reacs: reactant fragments ['A','B']
        :type frag_reacs: list(str)
        :param hotprod: hot dissociating product 
        :type hotprod: str
        :param otherprod: remaining product
        :type otherprod: str
        :return rename_ktp_dct: dct with new keys
        :rtype: {rxn_name: {P: (T, k)}}
    """
    rename_ktp_dct = {}
    reacs = '+'.join(frag_reacs)

    for sp in ktp_dct.keys():
        linker = (sp == hotprod)*'=' + \
            (sp != hotprod)*'=>'
        if len(hotsp_dct[sp]) > 1:
            hotprod = '+'.join(hotsp_dct[sp])
        elif len(hotsp_dct[sp]) == 1:
            hotprod = hotsp_dct[sp][0]
        prods = '+'.join([hotprod, otherprod])
        newkey = linker.join([reacs, prods])
        rename_ktp_dct[newkey] = ktp_dct[sp]

    return rename_ktp_dct
