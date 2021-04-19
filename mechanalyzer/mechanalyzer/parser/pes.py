"""
Extract PES and SUBPESs from a given mechanism
"""

import pandas as pd
import numpy as np
from mechanalyzer.parser import util

# functions working with dictionaries
# FUNTIONS FOR THE PES DICT OBJECTS CONTAINING INFO FOR THE REACTIONS ON PES


def build_pes_dct(formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst):
    """ Build a dictionary of the PESs
    """

    pes_dct = {}
    current_formula = ''
    for fidx, formula in enumerate(formula_str_lst):
        if current_formula == formula:
            pes_dct[formula]['rct_names_lst'].append(rct_names_lst[fidx])
            pes_dct[formula]['prd_names_lst'].append(prd_names_lst[fidx])
            pes_dct[formula]['rxn_name_lst'].append(rxn_name_lst[fidx])
        else:
            current_formula = formula
            pes_dct[formula] = {}
            pes_dct[formula]['rct_names_lst'] = [rct_names_lst[fidx]]
            pes_dct[formula]['prd_names_lst'] = [prd_names_lst[fidx]]
            pes_dct[formula]['rxn_name_lst'] = [rxn_name_lst[fidx]]

    return pes_dct


def connected_channels_dct(pes_dct):
    """ Determine all the connected reaction channels for each PES
        Build a dictionary for each PES with lists of connected channels:
            dct[PES_FORMULA] = [ [SUB_PES_1], [SUB_PES_2], ... , [SUB_PES_N] ]
            where each SUB_PES = [n1, n2, ... , nN],
            where n1 to nN correspond to ixds for channels that are
            connected to each other
        For efficiency we only determine channels for PESs we wish to run.
    """
    conn_chn_dct = {}
    for _, formula in enumerate(pes_dct):
        # Set the names lists for the rxns and species needed below
        pes_rct_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_lst = pes_dct[formula]['prd_names_lst']
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']

        connchnls = find_conn_chnls(
            pes_rct_lst, pes_prd_lst, pes_rxn_name_lst)

        # Add connected channels list to the dictionary
        conn_chn_dct[formula] = connchnls

    return conn_chn_dct


def find_conn_chnls(pes_rct_lst, pes_prd_lst, pes_rxn_name_lst):
    '''
    Given rxn names, reactants, products belonging to 1 PES:
    generate SUB PESs dictionaries
    conndct = {0: [['S1','S2'],['S3','S4'],[S5']], 1:[['S6','S7'],['S8','S9']]}
    connchnls = {0: [0,1,2] , 1:[3,4,5]}...
    corresponding to each subpes
    '''
    # preprocessing:
    # order (bimol) reactants and products in the same fashion
    # example if you have A+B and B+A they will be ordered in the same way
    pes_rct_lst = util.order_rct_bystoich(pes_rct_lst)
    pes_prd_lst = util.order_rct_bystoich(pes_prd_lst)
    # put everything in a dataframe. indices are the numbers given by enumerate
    len_rct_prd = np.array(list(map(len, pes_rct_lst))) + \
        np.array(list(map(len, pes_prd_lst)))
    pes_df = pd.DataFrame(np.array([pes_rct_lst, pes_prd_lst, len_rct_prd], dtype=object).T,
                          index=np.arange(0, len(pes_rxn_name_lst)),
                          columns=['rcts', 'prds', 'N_rcts_prds'])
    # order according to the total number of species (N of reactants + N of products)
    pes_df = pes_df.sort_values(by='N_rcts_prds')
    # Split up channels into a connected sub-pes within a formula
    subpes_idx = 0
    conndct = {}
    connchnls = {}

    for chnl_idx in pes_df.index:
        connected_to = []
        chnl_species = [list(pes_df['rcts'][chnl_idx]),
                        list(pes_df['prds'][chnl_idx])]

        for conn_chnls_idx in conndct:
            for spc_pair in chnl_species:
                if len(spc_pair) == 1 and spc_pair in conndct[conn_chnls_idx]:
                    # this works for unimol species; need also to verify bimol wellskipping channels
                    if conn_chnls_idx not in connected_to:
                        connected_to.append(conn_chnls_idx)
                elif len(spc_pair) == 1 and spc_pair[::-1] in conndct[conn_chnls_idx]:
                    if conn_chnls_idx not in connected_to:
                        connected_to.append(conn_chnls_idx)

            if len(chnl_species[0]) == 2 and len(chnl_species[1]) == 2:
                # bimol bimol reactions
                if ((chnl_species[0] in conndct[conn_chnls_idx]) and
                        (chnl_species[1] in conndct[conn_chnls_idx]) and
                        (conn_chnls_idx not in connected_to)):
                    connected_to.append(conn_chnls_idx)

        if not connected_to:
            conndct[subpes_idx] = chnl_species
            connchnls[subpes_idx] = [chnl_idx]
            subpes_idx += 1
        else:
            conndct[connected_to[0]].extend(chnl_species)
            connchnls[connected_to[0]].append(chnl_idx)
            if len(connected_to) > 1:
                for cidx, cval in enumerate(connected_to):
                    if cidx > 0:
                        conn_specs = conndct.pop(cval, None)
                        conn_chnls = connchnls.pop(cval, None)
                        conndct[connected_to[0]].extend(conn_specs)
                        connchnls[connected_to[0]].extend(conn_chnls)
            for cidx in conndct:
                conndct[cidx].sort()
                conndct[cidx] = [
                    conndct[cidx][i] for i in
                    range(len(conndct[cidx])) if i == 0 or
                    conndct[cidx][i] != conndct[cidx][i-1]]

    return connchnls

# ORIGINALLY IN MECHDRIVER
# FUNTIONS FOR THE PES DICT OBJECTS CONTAINING INFO FOR THE REACTIONS ON PES
def build_pes_idx_dct(pes_dct):
    """ build a dct relating index to formulas
    """
    idx_dct = {}
    form_dct = {}
    for pes_idx, formula in enumerate(pes_dct):
        idx_dct[pes_idx+1] = formula
        form_dct[formula] = pes_idx+1

    return idx_dct, form_dct


def reduce_pes_dct_to_user_inp(pes_dct, pesnums):
    """ get a pes dictionary containing only the PESs the user is running
    """
    run_pes_dct = {}
    for pes_idx, formula in enumerate(pes_dct):
        if pes_idx+1 in pesnums:
            run_pes_dct[formula] = pes_dct[formula]
    return run_pes_dct


def print_pes_channels(pes_dct):
    """ Print the PES
    """

    print('\n  Sorted Mechanism read from file:')
    for pes_idx, formula in enumerate(pes_dct):
        print('! PES:', pes_idx+1, formula)
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        for chn_idx, _ in enumerate(pes_rxn_name_lst):
            print('  {} = {}   1.0 0.0 0.0'.format(
                ' + '.join(pes_rct_names_lst[chn_idx]),
                ' + '.join(pes_prd_names_lst[chn_idx])))


def pes_dct_w_rxn_lsts(pes_dct, idx_dct, form_dct,
                       conn_chnls_dct, pes_idxs):
    """ Form a new PES dictionary with the rxn_lst formatted to work
        with the drivers currently
    """
    run_pes_dct = {}
    for formula in pes_dct:

        # Set correct pes index based on the formula
        pes_idx = form_dct[formula]

        # Build the names list
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']

        # Get a list of the idxs corresponding to which channels to run
        run_chnls = pes_idxs[pes_idx]

        # Select names from the names list corresponding to chnls to run
        for sub_pes_idx, sub_chnl_idxs in conn_chnls_dct[formula].items():
            rct_names_lst = []
            prd_names_lst = []
            rxn_name_lst = []
            rxn_chn_idxs = []
            for chn_idx in run_chnls:
                if chn_idx-1 in sub_chnl_idxs:
                    rct_names_lst.append(pes_rct_names_lst[chn_idx-1])
                    prd_names_lst.append(pes_prd_names_lst[chn_idx-1])
                    rxn_name_lst.append(pes_rxn_name_lst[chn_idx-1])
                    # rxn_model_lst.append(run_obj_dct[(pes_idx, chn_idx)])  # need dif way to get kin,spc models
                    rxn_chn_idxs.append(chn_idx)

            # Form reaction list (is empty if no chnls requested on sub pes)
            rxn_lst = format_run_rxn_lst(
                rct_names_lst, prd_names_lst, rxn_chn_idxs)

            # Add the rxn lst to the pes dictionary if there is anythin
            if rxn_lst:
                run_pes_dct[(formula, pes_idx, sub_pes_idx+1)] = rxn_lst

    return run_pes_dct


def format_run_rxn_lst(rct_names_lst, prd_names_lst, rxn_chn_idxs):
    """ Get the lst of reactions to be run
    """

    # Get a list of all the species in the pes
    spc_queue = []
    for idx, _ in enumerate(rct_names_lst):
        rxn_spc = list(rct_names_lst[idx])
        rxn_spc.extend(list(prd_names_lst[idx]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)

    # Now loop over all the reactions to build rxn_lst
    run_lst = []
    for idx, _ in enumerate(rct_names_lst):
        spc_queue = []
        rxn_spc = list(rct_names_lst[idx])
        rxn_spc.extend(list(prd_names_lst[idx]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)
        run_lst.append(
            {'species': spc_queue,
             'reacs': list(rct_names_lst[idx]),
             'prods': list(prd_names_lst[idx]),
             # 'model': rxn_model_lst[idx],
             'chn_idx': rxn_chn_idxs[idx],
             'dummy': []}
        )

    return run_lst
