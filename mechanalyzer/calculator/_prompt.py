""" Return the rates from a prompt dissociation process
"""
from operator import mod
import sys
import copy
import mess_io
import mechanalyzer


def multipes_prompt_dissociation_ktp_dct(list_strs_dct, model, bf_thresh):
    """ generation of rxn_ktp_dct for multipes with unknown character
        - identifies if ped or hot
        - extracts prompt rates for all
        - provides a final non-redundant ktp dct
        :param list_strs_dct: list of the dictionaries for each rxn
                containing the corresponding file strings
        :type list_strs_dct: [{'inp': str, 'ktp_out':str, ...}, {}, ...]
    """
    if isinstance(model, list):
        model = model[0]
    # SET UP A LIST OF PED FILES AND HOT FILES DEPENDING ON THE INPUT
    rxn_type_dct = {'ped': [], 'pedhot': [], 'hot': []}
    # base ktp dct: valid for any model. contains rxns on hot PESs (unmodified by prompt)
    # and PED rxns
    rxn_ktp_dct_full0 = {}
    for strs_dct in list_strs_dct:
        ped_spc, _ = mess_io.reader.ped.ped_names(strs_dct['inp'])
        hot_spc_en = mess_io.reader.hoten.get_hot_species(strs_dct['inp'])
        if ped_spc and hot_spc_en:
            rxn_type_dct['pedhot'].append(strs_dct)
        elif ped_spc and not hot_spc_en:
            rxn_type_dct['ped'].append(strs_dct)
        elif hot_spc_en and not ped_spc:
            rxn_type_dct['hot'].append(strs_dct)

        rxn_ktp_dct_full0.update(extract_ktp_dct(strs_dct['ktp_out']))

    # FOR EACH PED: CHECK IF OVERLAP WITH INTERMEDIATE PESs;
    # THEN FOR THE ONE IDENTIFIED CHECK OVERLAP WITH HOT OR INTERMEDIATE PESs, ITERATIVELY
    # GENERATE LISTS OF INTERMEDIATE PESs CONNECTED TOGETHER - IDK ABOUT "BRANCHES" FROM MULTIPLE HOT PESs
    # WHEN YOU FIND ONE "SPLIT" (MORE THAN 1 COMPATIBILITY)
    # REMEMBER THAT YOU CAN HAVE DIFFERENT OVERLAPS, SO IN THE END YOU CAN HAVE SOMETHING LIKE
    # LIST OF DCTS WITH CHAINS [CHAIN1, CHAIN2, CHAIN3...]
    # {PED:0, PEDHOT:1,2,3, HOT:4}
    # {PED:0, PEDHOT:5, HOT:6} for  now not allowed - only 1 chain
    # {PED:0, HOT:7}
    # {PED:8: HOT:9} THIS OCCURS IF YOU HAVE MULTIPLE HOT PRODUCTS WHICH DO DIFFERENT THINGS
    # RUN EACH SET - THE RXN_KTP_DCT_FULL MUST BE ALWAYS THE SAME - UPDATE IT AT THE END OF THE CHAIN
    # IN THESE CASES AN 'INTERMEDIATE' PES CAN NEVER BE A BEGINNING OR AN END
    rxn_ktp_dct_full = rxn_chains_calc(
        rxn_type_dct['ped'], rxn_type_dct['pedhot'], rxn_type_dct['hot'],
        model, bf_thresh)
    print(rxn_ktp_dct_full.keys(), rxn_ktp_dct_full0.keys())
    # Add in the prompt versions of the reactions
    # pop values in rxn_ktp_dct_full0
    ks = rxn_ktp_dct_full.keys()
    [rxn_ktp_dct_full0.pop(x) for x in ks if x in rxn_ktp_dct_full0.keys()]
    rxn_ktp_dct_full.update(
        rxn_ktp_dct_full0)
    # alternatively, you can update without popping but you need to update
    # the rxn_dct_red with the other one and then replace the latter

    return rxn_ktp_dct_full


def rxn_chains_calc(PES_0, PES_i, PES_end, model, bf_thresh):
    """ generate reaction chains of successive prompt/hot PESs
    """
    def listconnected(PES_str0, PESs_str):
        retlist = []
        [retlist.append(PES_str)
         for PES_str in PESs_str if isconnected(PES_str0, PES_str)]
        return retlist

    def isconnected(ped_spc_str, hot_en_spc_str):
        """ return true if ped species is hot species of the following PES """
        connected = []
        ped_spc_all, _ = mess_io.reader.ped.ped_names(ped_spc_str['inp'])
        hot_spc_en = mess_io.reader.hoten.get_hot_species(
            hot_en_spc_str['inp'])
        for ped_spc in ped_spc_all:
            frags = []
            [frags.extend(fragset.split('+')) for fragset in ped_spc]
            connected.extend(
                list(set(list(hot_spc_en.keys())).intersection(frags)))

        return len(connected) > 0

    rxn_ktp_dct_full = {}
    # SIMPLE PED/HOT PESs: COUPLE ANY N OF PED AND HOT
    # INTERMEDIATE PEDs: JUST ALLOW 1 REACTION CHAIN FOR NOW
    for pes0 in PES_0:
        rxn_ktp_dct_ped = {}  # ped dct

        hotend = listconnected(pes0, PES_end)
        if len(hotend) > 0:
            for hot_strs_dct in hotend:
                # call the usual function of the prompt dissociation
                rxn_ktp_dct, _, _ = prompt_dissociation_ktp_dct(
                    pes0['inp'], pes0['ktp_out'],
                    pes0['ped'], pes0['ke_out'],
                    hot_strs_dct['inp'], hot_strs_dct['log'],
                    model, bf_thresh)

                # update dct
                rxn_ktp_dct_ped.update(
                    rxn_ktp_dct)

        pedhot = listconnected(pes0, PES_i)

        if len(pedhot) == 1:
            # create chain
            pedhot_all = copy.deepcopy(PES_i)
            pedhot_i = pedhot[0]
            pedhot_all.remove(pedhot_i)

            # PED: PES0, HOT: PEDHOT_I - EXTRACT RATES AND NEW ENERGY DISTRIBUTION
            # call the usual function of the prompt dissociation
            rxn_ktp_dct, pednew, ene_bw = prompt_dissociation_ktp_dct(
                pes0['inp'], pes0['ktp_out'],
                pes0['ped'], pes0['ke_out'],
                pedhot_i['inp'], pedhot_i['log'],
                model, bf_thresh, hot_ped_str=pedhot_i['ped'], hot_ke_out_str=pedhot_i['ke_out'])
            # you need the reaction ktp dictionary as a new input
            # you need the starting energy distribution of the products of pedhot_i
            prompt_ktp_dct = copy.deepcopy(rxn_ktp_dct)
            # rxn_ktp_dct contains all reactions in the chain,
            # prompt_ktp_dct contains only the rxns to be "expanded"

            while len(hotend) == 0:  # until you find a termination to the chain

                pedhot_new = listconnected(pedhot_i, pedhot_all)

                if len(pedhot_new) == 1:
                    pedhot_i = copy.deepcopy(pedhot_new[0])
                    pedhot_all.remove(pedhot_i)

                    prompt_ktp_dct, pednew, ene_bw = prompt_chain_ktp_dct(
                        prompt_ktp_dct, pednew,
                        pedhot_i['inp'], pedhot_i['log'],
                        model, bf_thresh, hot_ped_str=pedhot_i['ped'],
                        hot_ke_out_str=pedhot_i['ke_out'], ene_bw=ene_bw)

                    # replace rxns every time
                    rxn_ktp_dct.update(prompt_ktp_dct)

                elif len(pedhot_new) > 1:
                    print('*Error: more than 1 hot chan starting from {}. Exiting')
                    sys.exit()

                # search for end of connection
                hotend = listconnected(pedhot_i, PES_end)
                if len(hotend) > 0:
                    hot_strs_dct = hotend[0]
                    print('CALL PROMPTFUNCTION HERE')
                    prompt_ktp_dct, _, _ = prompt_chain_ktp_dct(
                        prompt_ktp_dct, pednew,
                        hot_strs_dct['inp'], hot_strs_dct['log'],
                        model, bf_thresh)

                    rxn_ktp_dct.update(prompt_ktp_dct)

            # update dct
            rxn_ktp_dct_ped.update(
                rxn_ktp_dct)

        elif len(pedhot) > 1:
            print('*Error: only chains with 1 connected pedhot species PES available')
            sys.exit()

        # merge the ped dct with the full one
        rxn_ktp_dct_full = mechanalyzer.calculator.rates.merge_rxn_ktp_dcts(
            rxn_ktp_dct_full, rxn_ktp_dct_ped
        )

    return rxn_ktp_dct_full


def prompt_chain_ktp_dct(rxn_ktp_dct, ene_start_df_dct,
                         pedhot_inp_str, pedhot_log_str,
                         model, bf_thresh, hot_ped_str=None,
                         hot_ke_out_str=None, ene_bw=0):
    """ Starts from rxn_ktp_dct
        rxn_ktp_dct: ktp dct to update with prompt fractions
        ene_start_df_dct: starting energy distribution for hotspecies (dct keys)

    """
    # HOTEN INFO
    hot_frag_dct, hot_spc_en, hoten_dct, fne_bf = \
        hot_info(pedhot_inp_str, pedhot_log_str)

    # Derive Branching Fractions, Calculate Prompt Rates
    # Merge Prompt Rates with Thermal Rates
    full_prompt_rxn_ktp_dct = {}
    print(ene_start_df_dct.keys(), hot_spc_en.keys())
    for spc in ene_start_df_dct.keys():

        # if fragment in hotspecies, proceed - might not be
        if spc not in hot_spc_en.keys():
            continue

        # CALCULATE PROMPT DISSOCIATION RATES K*BF
        # find rxn in the ktp dct to replace
        for rxn_key in rxn_ktp_dct.keys():
            if spc in rxn_key[1]:
                rxn = '+'.join(rxn_key[0])+'='+'+'.join(rxn_key[1])
                ktp_label = rxn_key
                frag_reacs = rxn_key[0]
                frag2 = tuple([sp for sp in rxn_key[1] if sp != spc])

        # JOIN PED AND HOTEN -> DERIVE PRODUCTS BF
        
        full_prompt_rxn_ktp_dct = calc_bf_ktp(full_prompt_rxn_ktp_dct, model, bf_thresh, 
                ene_start_df_dct[spc], fne_bf[spc], frag2, frag_reacs, 
                hoten_dct[spc], rxn, rxn_ktp_dct[ktp_label], hot_frag_dct)
        
        # if hot options: it means hoten has also ped - extract info
        pedhot_df_dct = {}
        if hot_ke_out_str and hot_ped_str:
            # frag1 è hotspecies - devi prendere quella distribuzione dalle hot
            pedhot_df_dct, ene_bw = build_pedhot_df_dct(pedhot_inp_str, hot_ped_str, hot_ke_out_str,
                                                        spc, ene_start_df_dct[spc], ene_bw, model)

    return full_prompt_rxn_ktp_dct, pedhot_df_dct, ene_bw

def prompt_dissociation_ktp_dct(ped_inp_str, ped_out_str,
                                ped_ped_str, ped_ke_out_str,
                                hot_inp_str, hot_log_str,
                                model, bf_thresh, hot_ped_str=None, hot_ke_out_str=None):
    """ Parses the MESS files and generates the rate constants
        from a prompt dissociation process
    """

    # PED INFO
    ped_spc, ped_dct, dof_dct, \
        dos_df, energy_dct = ped_info(
            ped_inp_str, ped_ped_str, ped_ke_out_str)

    # HOTEN INFO
    hot_frag_dct, hot_spc_en, hoten_dct, fne_bf = \
        hot_info(hot_inp_str, hot_log_str)

    # OBTAIN ALL OF THE RATE CONSTANTS FROM THE OUTPUT FILES
    # put dictionaries together
    rxn_ktp_dct = extract_ktp_dct(ped_out_str)

    # Derive Branching Fractions, Calculate Prompt Rates
    # Merge Prompt Rates with Thermal Rates
    full_prompt_rxn_ktp_dct = {}

    for spc in ped_spc:
        rxn = '='.join(spc)
        reacs, prods = spc
        frag_reacs = tuple(reacs.split('+'))
        frag_prods = tuple(prods.split('+'))
        label = ((reacs,), (prods,), (None,))
        ktp_label = (frag_reacs, frag_prods, (None,))
        # select the fragment of which you want the PED:
        # it is the one in common with hotspecies

        try:
            frag1 = list(set(hot_spc_en).intersection(frag_prods))[0]
            frag2 = list(set(frag_prods).difference((frag1,)))[0]
        except IndexError:
            print('no superposition between PED fragments and hot fragments '
                  '- skipping \n')
            continue

        ped_df = ped_dct[label]
        ene_bw = energy_dct[reacs] - energy_dct[prods]

        if ene_bw < 0:
            print(
                'Warning: endothermic reaction: whatever the model, only fne will be used')
            model = 'fne'

        # DERIVE PED OF THE HOT FRAGMENT
        if model != 'fne':
            ped_df_frag1 = mechanalyzer.builder.ped.ped_frag1(
                ped_df, frag1, frag2, model,
                dos_df=dos_df, dof_info=dof_dct[prods])
        else:
            ped_df_frag1 = None

        full_prompt_rxn_ktp_dct = calc_bf_ktp(full_prompt_rxn_ktp_dct, model, bf_thresh, 
                        ped_df_frag1, fne_bf[frag1], (frag2,), frag_reacs, 
                        hoten_dct[frag1], rxn, rxn_ktp_dct[ktp_label], hot_frag_dct)
        
        # if hot options: it means hoten has also ped - extract info
        pedhot_df_dct = {}
        if hot_ke_out_str and hot_ped_str:
            # frag1 è hotspecies - devi prendere quella distribuzione dalle hot
            pedhot_df_dct, ene_bw = build_pedhot_df_dct(hot_inp_str, hot_ped_str, hot_ke_out_str,
                                                        frag1, ped_df_frag1, ene_bw, model)

    return full_prompt_rxn_ktp_dct, pedhot_df_dct, ene_bw


def calc_bf_ktp(full_prompt_rxn_ktp_dct, model, bf_thresh, 
                ped_df_frag1, fne_bf_frag1, frag2, frag_reacs, 
                hoten_dct_frag, rxn, ktp_dct, hot_frag_dct):
    """ generate full bf and ktp dct from ped and hoten
        NB frag2 must be tuple
    """
    # JOIN PED AND HOTEN -> DERIVE PRODUCTS BF
    bf_tp_dct = mechanalyzer.builder.bf.bf_tp_dct(
        model, ped_df_frag1, hoten_dct_frag, bf_thresh,
        savefile=True, rxn=rxn, fne=fne_bf_frag1)

    # CALCULATE PROMPT DISSOCIATION RATES K*BF

    prompt_rxn_ktp_dct = mechanalyzer.builder.bf.merge_bf_ktp(
        bf_tp_dct, ktp_dct,
        frag_reacs, frag2, hot_frag_dct)

    # Merge Prompt Rates with all current; Rates added if rxn is prev. found
    # print(modeltype, full_prompt_rxn_ktp_dct[modeltype], '\n', prompt_rxn_ktp_dct[modeltype], '\n','\n','stop')
    full_prompt_rxn_ktp_dct = mechanalyzer.calculator.rates.merge_rxn_ktp_dcts(
        full_prompt_rxn_ktp_dct,
        prompt_rxn_ktp_dct
    )

    return full_prompt_rxn_ktp_dct

def build_pedhot_df_dct(hot_inp_str, hot_ped_str, hot_ke_out_str,
                        starthotfrag, starthotfrag_df, ene_bw, model):

    pedhot_df_dct = {}

    hot_ped_spc, hot_ped_dct, hot_dof_dct, \
        hot_dos_df, hot_energy_dct = ped_info(
            hot_inp_str, hot_ped_str, hot_ke_out_str)
    # starting energy distribution
    starthot_df = starthotfrag_df

    # just 1 ped set in this case
    _, prods = hot_ped_spc[0]
    label = ((starthotfrag,), (prods,), (None,))

    # calculate new PED (beta)
    ped_df_fromhot = hot_ped_dct[label]

    ped_df_rescaled = mechanalyzer.calculator.bf.ped_df_rescale(
        starthot_df, ped_df_fromhot)
    ene_bw += (hot_energy_dct[starthotfrag] - hot_energy_dct[prods])

    if ene_bw < 0:
        print(
            'Warning: endothermic reaction: whatever the model, only fne will be used')
        model = 'fne'

    frag1, frag2 = prods.split('+')
    print('\n', ped_df_rescaled, '\n')
    # DERIVE PED OF THE HOT FRAGMENT - BOTH , CHECK FRAG2 NOT ATOM
    pedhot_df_dct[frag1] = mechanalyzer.builder.ped.ped_frag1(
        ped_df_rescaled, frag1, frag2, model,
        dos_df=hot_dos_df, dof_info=hot_dof_dct[prods])
    if hot_dof_dct[prods]['n_atoms'][frag2] > 1:
        pedhot_df_dct[frag2] = mechanalyzer.builder.ped.ped_frag1(
            ped_df_rescaled, frag2, frag1, model,
            dos_df=hot_dos_df, dof_info=hot_dof_dct[prods])

    return pedhot_df_dct, ene_bw


def ped_info(ped_inp_str, ped_ped_str, ped_ke_out_str):
    """ file strings
        wellreac: str if you want PEDs of well->bimol 
    """
    # Read: ped.inp:
    #  species names, energies
    spc_blocks_ped = mess_io.reader.get_species(ped_inp_str)
    ped_spc, _ = mess_io.reader.ped.ped_names(ped_inp_str)  # can supply
    energy_dct, _, _, _ = mess_io.reader.pes(ped_inp_str)
    # if wellreac:
    # Read: rate_ped.out and ke_ped.out:
    #  energy barriers, dofs, fragment names
    dof_dct = {}
    for spc in ped_spc:
        _, prods = spc

        # Derive dofs involved
        dof_dct[prods] = mechanalyzer.calculator.statmodels.get_dof_info(
            spc_blocks_ped[prods])

    # Read ped.out file for product energy distributions
    ped_dct = mess_io.reader.ped.get_ped(
        ped_ped_str, ped_spc, energy_dct, sp_labels='auto')

    # Read ke_ped.out file for energy density of each fragment
    dos_df = mess_io.reader.rates.dos_rovib(ped_ke_out_str, sp_labels='auto')

    return ped_spc, ped_dct, dof_dct, dos_df, energy_dct


def hot_info(hot_inp_str, hot_log_str):
    """
    Extract required info for hotenergies
    """
    spc_blocks_hoten = mess_io.reader.get_species(hot_inp_str)
    hot_frag_dct = mess_io.reader.dct_species_fragments(spc_blocks_hoten)
    hot_spc_en = mess_io.reader.hoten.get_hot_species(hot_inp_str)

    hoten_dct = mess_io.reader.hoten.extract_hot_branching(
        hot_log_str, hot_spc_en, list(spc_blocks_hoten.keys()), sp_labels='auto')

    fne_bf = mess_io.reader.hoten.extract_fne(hot_log_str, sp_labels='auto')

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
