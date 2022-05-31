""" Return the rates from a prompt dissociation process
"""
from operator import mod
import sys
import copy
import mess_io
import mechanalyzer


def prompt_dissociation_ktp_dct(ped_inp_str, ped_out_str,
                                ped_ped_str, ped_ke_out_str,
                                hot_inp_str, hot_log_str,
                                model, bf_thresh, hot_ped_str=None, hot_ke_out_str=None):
    """ Parses the MESS files and generates the rate constants
        from a prompt dissociation process
    """

    # PED INFO
    _, ped_dct, dof_dct, \
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
    ene_bw_dct = {}
    pedhot_df_dct = {}

    for label, ped_df in ped_dct.items():
        # start prompt chains from wells/bimol with ground energies:
        # A+B=>C*+D*, (A+B=>)E=>C*+D*
        frag_reacs = label[0]
        frag_prods = label[1]
        reacs = '+'.join(frag_reacs)
        prods = '+'.join(frag_prods)
        rxn = '='.join([reacs, prods])

        ############################ various checks and derivation of frag names ###################
        # skip PEDs going to other wells or PEDs coming from hotenergy for some reason
        # - which means that the ped_df[T][p][energy] is a Series and not a float

        if len(frag_prods) < 2 or not isinstance(ped_df.iloc[0].iloc[0].iloc[0], float):
            continue
        # select the fragment of which you want the PED:
        try:
            frag1 = list(set(hot_spc_en).intersection(frag_prods))[0]
            frag2 = list(set(frag_prods).difference((frag1,)))[0]
        except IndexError:
            print('no superposition between PED fragments {} and hot fragments '
                  '- skipping \n'.format(prods))
            continue

        ene_bw_dct[label] = energy_dct[reacs] - energy_dct[prods]
        if ene_bw_dct[label] < 0 and len(frag_reacs) > 1:
            print(
                'Warning: endothermic reaction {} with DH of {:.2f} kcal/mol: \
                    whatever the model, only fne will be used'.format(label, -ene_bw_dct[label]))
            model = 'fne'
        elif ene_bw_dct[label] < 0 and len(frag_reacs) == 1:
            print('Endothermic reaction {} with DH of {:.2f} kcal/mol excluded from loop'.format(
                label, -ene_bw_dct[label]))
            continue
        #################################################################################################

        # DERIVE PED OF THE HOT FRAGMENT
        print(
            'Processing reaction {} with DH of {:.2f} kcal/mol'.format(label, -ene_bw_dct[label]))
        if model != 'fne':
            ped_df_frag1 = mechanalyzer.builder.ped.ped_frag1(
                ped_df, frag1, frag2, model,
                dos_df=dos_df, dof_info=dof_dct[prods])
        else:
            ped_df_frag1 = None

        # calculate prompt fraction
        full_prompt_rxn_ktp_dct = calc_bf_ktp(full_prompt_rxn_ktp_dct, model, bf_thresh,
                                              ped_df_frag1, fne_bf[frag1], label,
                                              hoten_dct[frag1], rxn, rxn_ktp_dct[label], hot_frag_dct)
        # print(full_prompt_rxn_ktp_dct.keys(), '\n')
        # if hot options: it means hoten has also ped - extract info

        if hot_ke_out_str and hot_ped_str:
            # frag1 è hotspecies - devi prendere quella distribuzione dalle hot
            pedhot_df_dct_spc, ene_bw_spc = build_pedhot_df_dct(hot_inp_str, hot_ped_str, hot_ke_out_str,
                                                                frag1, ped_df_frag1, ene_bw_dct[label], model, label)
            pedhot_df_dct.update(pedhot_df_dct_spc)
            ene_bw_dct.update(ene_bw_spc)

    return full_prompt_rxn_ktp_dct, pedhot_df_dct, ene_bw_dct


def prompt_chain_ktp_dct(rxn_ktp_dct, ene_start_df_dct,
                         pedhot_inp_str, pedhot_log_str,
                         model, bf_thresh, ene_bw_dct, hot_ped_str=None,
                         hot_ke_out_str=None):
    """ Starts from rxn_ktp_dct
        rxn_ktp_dct: ktp dct to update with prompt fractions
        ene_start_df_dct: starting energy distribution for given reaction chain
                     {(reacs,),(prods,),(None,): {frag1: energy_distr, frag2: energy_distr}

    """
    # HOTEN INFO
    hot_frag_dct, hot_spc_en, hoten_dct, fne_bf = \
        hot_info(pedhot_inp_str, pedhot_log_str)

    # Derive Branching Fractions, Calculate Prompt Rates
    # Merge Prompt Rates with Thermal Rates
    full_prompt_rxn_ktp_dct = {}
    # print(ene_start_df_dct.keys(), hot_spc_en.keys())
    for label, ene_start_df in ene_start_df_dct.items():

        # find dissociating fragment if available
        frag10, frag20 = ene_start_df.keys()
        frag = frag10*(frag10 in hot_spc_en.keys()) + \
            frag20*(frag20 in hot_spc_en.keys())

        # continue if no overlap
        if not frag:
            continue

        rxn = '+'.join(label[0])+'='+'+'.join(label[1])
        # model checks
        if ene_bw_dct[label] < 0:
            print(
                'Warning: endothermic reaction {} with DH of {:.2f} kcal/mol: \
                    whatever the model, only fne will be used'.format(label, -ene_bw_dct[label]))
            model = 'fne'
        
        full_prompt_rxn_ktp_dct = calc_bf_ktp(full_prompt_rxn_ktp_dct, model, bf_thresh,
                                              ene_start_df[frag], fne_bf[frag], label,
                                              hoten_dct[frag], rxn, rxn_ktp_dct[label], hot_frag_dct)

        # if hot options: it means hoten has also ped - extract info
        pedhot_df_dct = {}
        if hot_ke_out_str and hot_ped_str:
            # frag1 è hotspecies - devi prendere quella distribuzione dalle hot
            pedhot_df_dct_spc, ene_bw_spc = build_pedhot_df_dct(pedhot_inp_str, hot_ped_str, hot_ke_out_str,
                                                                frag, ene_start_df[frag], ene_bw_dct[label], model, label)
            pedhot_df_dct.update(pedhot_df_dct_spc)
            ene_bw_dct.update(ene_bw_spc)

    return full_prompt_rxn_ktp_dct, pedhot_df_dct, ene_bw_dct


def build_pedhot_df_dct(hot_inp_str, hot_ped_str, hot_ke_out_str,
                        starthotfrag, starthotfrag_df, ene_bw, model, label_start):
    """ label_start: label of the starting reaction (reaction generating the radical)
    """
    
    pedhot_df_dct_tot = {}
    pedhot_df_dct = {}
    ene_bw_dct = {}

    hot_ped_spc, hot_ped_dct, hot_dof_dct, \
        hot_dos_df, hot_energy_dct = ped_info(
            hot_inp_str, hot_ped_str, hot_ke_out_str)
    # starting energy distribution
    starthot_df = starthotfrag_df

    # just 1 ped set in this case
    _, prods = hot_ped_spc[0]
    frag1, frag2 = prods.split('+')
    label = ((starthotfrag,), tuple(prods.split('+')), (None,))

    newprods = tuple([pr for pr in label_start[1] if pr !=
                      starthotfrag]) + tuple(prods.split('+'))
    
    newlabel = (label_start[0], newprods, label_start[2])
    ene_bw += (hot_energy_dct[starthotfrag] - hot_energy_dct[prods])
    
    print('deriving hot fragment distribution for {} \n'.format(newlabel))
    if ene_bw < 0:
        print(
            'Warning: endothermic reaction {} with DH of {:.2f} kcal/mol: \
                whatever the model, only fne will be used'.format(newlabel, -ene_bw))
        model = 'fne'
        pedhot_df_dct[frag1] = {}
        pedhot_df_dct[frag2] = {}
        
    # calculate new PED (beta)
    ped_df_fromhot = hot_ped_dct[label]

    ped_df_rescaled = mechanalyzer.calculator.bf.ped_df_rescale(
        starthot_df, ped_df_fromhot)

    ########################################################################################
    # stupid test to delete later: try to just have the starthot_df but rescale the energy
    ped_df_rescaled = mechanalyzer.calculator.bf.ped_df_rescale_test(
        starthot_df, hot_energy_dct[prods]-hot_energy_dct[starthotfrag])

    #########################################################################################

    # assign name to the global reaction
    if model != 'fne':

        # DERIVE PED OF THE HOT FRAGMENT - BOTH , CHECK FRAG2 NOT ATOM
        pedhot_df_dct[frag1] = mechanalyzer.builder.ped.ped_frag1(
            ped_df_rescaled, frag1, frag2, model,
            dos_df=hot_dos_df, dof_info=hot_dof_dct[prods])

        if hot_dof_dct[prods]['n_atoms'][frag2] > 1:
            pedhot_df_dct[frag2] = mechanalyzer.builder.ped.ped_frag1(
                ped_df_rescaled, frag2, frag1, model,
                dos_df=hot_dos_df, dof_info=hot_dof_dct[prods])

    pedhot_df_dct_tot[newlabel] = copy.deepcopy(pedhot_df_dct)
    ene_bw_dct[newlabel] = ene_bw

    return pedhot_df_dct_tot, ene_bw_dct


def calc_bf_ktp(full_prompt_rxn_ktp_dct, model, bf_thresh,
                ped_df_frag1, fne_bf_frag1, label,
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
        label, hot_frag_dct)

    # Merge Prompt Rates with all current; Rates added if rxn is prev. found
    # print(modeltype, full_prompt_rxn_ktp_dct[modeltype], '\n', prompt_rxn_ktp_dct[modeltype], '\n','\n','stop')
    full_prompt_rxn_ktp_dct = mechanalyzer.calculator.rates.merge_rxn_ktp_dcts(
        full_prompt_rxn_ktp_dct,
        prompt_rxn_ktp_dct
    )

    return full_prompt_rxn_ktp_dct


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
        ped_ped_str, energy_dct, sp_labels='auto')

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
        out_str, filter_kts=True,
        filter_reaction_types=('fake', 'self',
                                       'loss', 'capture', 'reverse'),
        relabel_reactions=True
    )
    
    return rxn_ktp_dct

def resort_ktp_labels(ktp_dct):
    """
    takes ktp dct
    resort labels of products in alphabetical order
    - only way to guarantee consistency in dct keys when applying prompt
    """
    
    new_ktp_dct = {}
    for key, val in ktp_dct.items():
        prods = list(key[1])
        prods.sort()
        new_key = (key[0], tuple(prods), key[2])
        new_ktp_dct[new_key] = copy.deepcopy(val)
        
    return new_ktp_dct