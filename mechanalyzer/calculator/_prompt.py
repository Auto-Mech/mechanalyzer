""" Return the rates from a prompt dissociation process
"""
import sys
import copy
import mess_io
import mechanalyzer


def multipes_prompt_dissociation_ktp_dct(list_strs_dct, models, bf_thresh):
    """ generation of rxn_ktp_dct for multipes with unknown character
        - identifies if ped or hot
        - extracts prompt rates for all
        - provides a final non-redundant ktp dct
        :param list_strs_dct: list of the dictionaries for each rxn
                containing the corresponding file strings
        :type list_strs_dct: [{'inp': str, 'ktp_out':str, ...}, {}, ...]
    """
    # SET UP A LIST OF PED FILES AND HOT FILES DEPENDING ON THE INPUT
    rxn_type_dct = {'ped': [], 'hot': []}
    # base ktp dct: valid for any model. contains rxns on hot PESs (unmodified by prompt)
    # and PED rxns
    rxn_ktp_dct_full = {}
    for strs_dct in list_strs_dct:
        ped_spc, _ = mess_io.reader.ped.ped_names(strs_dct['inp'])
        hot_spc_en = mess_io.reader.hoten.get_hot_species(strs_dct['inp'])
        if ped_spc:
            rxn_type_dct['ped'].append(strs_dct)
        elif hot_spc_en:
            rxn_type_dct['hot'].append(strs_dct)
        else:
            print('Error: no PED/HOT pes, useless to call this fct - exiting')
            sys.exit()
        rxn_ktp_dct_full.update(extract_ktp_dct(strs_dct['ktp_out']))

    rxn_ktp_dct_models_full = dict(zip(models, [{}]*len(models)))

    for ped_strs_dct in rxn_type_dct['ped']:
        rxn_ktp_dct_models_ped = dict(zip(models, [{}]*len(models)))
        for hot_strs_dct in rxn_type_dct['hot']:
            # call the usual function of the prompt dissociation
            rxn_ktp_dct_models = prompt_dissociation_ktp_dct(
                ped_strs_dct['inp'], ped_strs_dct['ktp_out'],
                ped_strs_dct['ped'], ped_strs_dct['ke_out'],
                hot_strs_dct['inp'], hot_strs_dct['log'],
                models, bf_thresh, fullktpout=False)

            # update dct
            for modeltype in models:
                rxn_ktp_dct_models_ped[modeltype].update(
                    rxn_ktp_dct_models[modeltype])

        # print('pedhot dct {}'.format(rxn_ktp_dct_models_ped[models[0]]))
        # merge the ped dct with the full one
        for modeltype in models:
            rxn_ktp_dct_models_full[modeltype] = mechanalyzer.calculator.rates.merge_rxn_ktp_dcts(
                rxn_ktp_dct_models_full[modeltype],
                rxn_ktp_dct_models_ped[modeltype]
            )
        # print('pedhot dct merged {}'.format(rxn_ktp_dct_models_full[models[0]]))

    # print('first model prompt dct {}'.format(rxn_ktp_dct_models_full[models[0]]))
    # Add in the prompt versions of the reactions

    for modeltype in models:
        # pop values in rxn_ktp_dct_full
        ks = rxn_ktp_dct_models_full[modeltype].keys()
        rxn_dct_red = copy.deepcopy(rxn_ktp_dct_full)
        [rxn_dct_red.pop(x) for x in ks if x in rxn_dct_red.keys()]
        rxn_ktp_dct_models_full[modeltype].update(
            rxn_dct_red)
        # alternatively, you can update without popping but you need to update
        # the rxn_dct_red with the other one and then replace the latter

    return rxn_ktp_dct_models_full


def prompt_dissociation_ktp_dct(ped_inp_str, ped_out_str,
                                ped_ped_str, ped_ke_out_str,
                                hot_inp_str, hot_log_str,
                                models, bf_thresh, hot_out_str=None,
                                fullktpout=False):
    """ Parses the MESS files and generates the rate constants
        from a prompt dissociation process
        hot_out_str not strictly necessary - if None/'', the ktp dct
        will only provide ktp dct with original PES and prompt values
        fullktpout: True/False: decide if you want the full set of rates
        (True, final and ready to fit) or if you want just the prompt rxns (False)
    """

    # PED INFO
    ped_spc, spc_blocks_ped, ped_dct, dof_dct, fragments_dct, \
        ene_bw_dct, dos_df, rxn_ktp_dct_ped = prompt_ped_info(
            ped_inp_str, ped_out_str, ped_ped_str, ped_ke_out_str)

    # HOTEN INFO
    hot_frag_dct, hot_spc_en, hoten_dct, fne_bf = \
        prompt_hot_info(hot_inp_str, hot_log_str)

    # OBTAIN ALL OF THE RATE CONSTANTS FROM THE OUTPUT FILES
    # put dictionaries together
    rxn_ktp_dct = {}
    rxn_ktp_dct.update(rxn_ktp_dct_ped)

    # check hot_out_str and if present include it in the dct
    if hot_out_str:
        rxn_ktp_dct_hot = extract_ktp_dct(hot_out_str)
        rxn_ktp_dct.update(rxn_ktp_dct_hot)

    # Derive Branching Fractions, Calculate Prompt Rates
    # Merge Prompt Rates with Thermal Rates
    prompt_rxns = ()
    full_prompt_rxn_ktp_dct = dict(zip(models, [{}]*len(models)))

    for spc in ped_spc:
        reacs, prods = spc
        _reacs = tuple(reacs.split('+'))
        _prods = tuple(prods.split('+'))

        label = ((reacs,), (prods,), (None,))
        relabel = (_reacs, _prods, (None,))

        ped_df = ped_dct[label]
        ene_bw = ene_bw_dct[label]
        # select the fragment of which you want the PED:
        # it is the one in common with hotspecies
        fragments = fragments_dct[label]
        try:
            frag1 = list(set(hot_spc_en).intersection(fragments))[0]
            frag2 = list(set(fragments).difference((frag1,)))[0]
        except IndexError:
            print('no superposition between PED fragments and hot fragments '
                  '- skipping \n')
            continue

        # DERIVE PED OF THE HOT FRAGMENT
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            ped_df, frag1, frag2, models,
            dos_df=dos_df, dof_info=dof_dct[label], ene_bw=ene_bw)

        # JOIN PED AND HOTEN -> DERIVE PRODUCTS BF
        bf_tp_dct = mechanalyzer.builder.bf.bf_tp_dct(
            models, ped_df_frag1_dct, hoten_dct[frag1], bf_thresh,
            savefile=True, reac=reacs, fne=fne_bf[frag1])

        # Calculate Prompt Dissociation Rates
        frag_reacs_dct = mess_io.reader.dct_species_fragments(
            spc_blocks_ped)
        frag_reacs = frag_reacs_dct[spc[0]]

        prompt_rxn_ktp_dct = mechanalyzer.builder.bf.merge_bf_ktp(
            bf_tp_dct, rxn_ktp_dct[relabel],
            frag_reacs, frag1, frag2, hot_frag_dct)

        prompt_rxns += (relabel,)

        for modeltype in models:
            # Merge Prompt Rates with all current; Rates added if rxn is prev. found
            full_prompt_rxn_ktp_dct[modeltype] = mechanalyzer.calculator.rates.merge_rxn_ktp_dcts(
                full_prompt_rxn_ktp_dct[modeltype],
                prompt_rxn_ktp_dct[modeltype]
            )

    if fullktpout:
        # Remove the original reaction
        for rxn in prompt_rxns:
            print('pop reaction?', rxn)
            rxn_ktp_dct.pop(rxn)

        # Add in the prompt versions of the reactions
        rxn_ktp_dct_models = dict(zip(models, [rxn_ktp_dct]*len(models)))
        for modeltype in models:
            rxn_ktp_dct_models[modeltype].update(
                full_prompt_rxn_ktp_dct[modeltype])

        return rxn_ktp_dct_models

    else:
        return full_prompt_rxn_ktp_dct


def prompt_ped_info(ped_inp_str, ped_out_str,
                    ped_ped_str, ped_ke_out_str):
    # Read: ped.inp:
    #  species names, energies
    spc_blocks_ped = mess_io.reader.get_species(ped_inp_str)
    ped_spc, _ = mess_io.reader.ped.ped_names(ped_inp_str)  # can supply
    energy_dct, _, conn_lst_dct, _ = mess_io.reader.pes(ped_inp_str)

    # Read: rate_ped.out and ke_ped.out:
    #  energy barriers, dofs, fragment names
    ene_bw_dct = {}
    dof_dct = {}
    fragments_dct = {}
    for spc in ped_spc:
        reacs, prods = spc
        label = ((reacs,), (prods,), (None,))

        # Find the corresponding energy barrier
        barrier_label = mess_io.reader.find_barrier(conn_lst_dct, reacs, prods)
        try:
            ene_bw_dct[label] = energy_dct[barrier_label]-energy_dct[prods]
        except KeyError:
            ene_bw_dct[label] = energy_dct[reacs]-energy_dct[prods]

        # Derive dofs involved
        dof_info = mechanalyzer.calculator.statmodels.get_dof_info(
            spc_blocks_ped[prods])
        dof_dct[label] = dof_info
        fragments_dct[label] = mess_io.reader.dct_species_fragments(
            spc_blocks_ped)[prods]

        max_ene = mechanalyzer.calculator.statmodels.max_en_auto(
            dof_info['n_atoms']['TS'], ene_bw_dct[label], ref_ene=energy_dct[prods])
        print('test _prompt line 133: max ene to write in ped/hoten is {} kcal/mol \n'.format(max_ene))
        print(
            'for hoten this would be {} kcal/mol \n'.format(max_ene-energy_dct[prods]))
    # Read ped.out file for product energy distributions
    ped_dct = mess_io.reader.ped.get_ped(
        ped_ped_str, ped_spc, energy_dct, sp_labels='inp')

    # Read ke_ped.out file for energy density of each fragment
    dos_df = mess_io.reader.rates.dos_rovib(ped_ke_out_str, sp_labels='inp')

    rxn_ktp_dct_ped = extract_ktp_dct(ped_out_str)

    return ped_spc, spc_blocks_ped, ped_dct, dof_dct, fragments_dct, ene_bw_dct, dos_df, rxn_ktp_dct_ped


def prompt_hot_info(hot_inp_str, hot_log_str):
    """
    Extract required info for hotenergies
    """
    spc_blocks_hoten = mess_io.reader.get_species(hot_inp_str)
    hot_frag_dct = mess_io.reader.dct_species_fragments(spc_blocks_hoten)
    hot_spc_en = mess_io.reader.hoten.get_hot_species(hot_inp_str)

    hoten_dct = mess_io.reader.hoten.extract_hot_branching(
        hot_log_str, hot_spc_en, list(spc_blocks_hoten.keys()), sp_labels='inp')

    fne_bf = mess_io.reader.hoten.extract_fne(hot_log_str)

    return hot_frag_dct, hot_spc_en, hoten_dct, fne_bf


def extract_ktp_dct(out_str):
    """ 
    Extract ktp dct
    """
    # extract rate constants
    rxn_ktp_dct = mess_io.reader.rates.get_rxn_ktp_dct(
        out_str,
        filter_kts=True,
        filter_reaction_types=('fake', 'self',
                                       'loss', 'capture', 'reverse'),
        relabel_reactions=True
    )

    return rxn_ktp_dct
