"""
Extract PES and SUBPESs from a given mechanism
"""

import pandas as pd
import numpy
from mechanalyzer.parser._util import order_rct_bystoich


# functions working with dictionaries
# FUNTIONS FOR THE PES DICT OBJECTS CONTAINING INFO FOR THE REACTIONS ON PES
def connected_channels_dct(pes_dct):
    """ Determine all the connected reaction channels for each PES
        Build a dictionary for each PES with lists of connected channels:
        dct[PES_FORMULA] = [ [SUB_PES_1], [SUB_PES_2], ... , [SUB_PES_N] ]
        where each SUB_PES = [n1, n2, ... , nN],
        where n1 to nN correspond to ixds for channels that are
        connected to each other.

        For efficiency we only determine channels for PESs we wish to run.

        :param pes_dct: Dictionary of PEss
        :type pes_dct:
        :rtype: dict[]
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
    """ Determine all of the connected reaction channels on a PES.
        Information compiled in SUB-PES dictionaries:

        conndct = {0: [['S1','S2'],['S3','S4'],[S5']],
                   1:[['S6','S7'],['S8','S9']]}
        connchnls = {0: [0,1,2] , 1:[3,4,5]}...

        :param pes_rct_lst:
    """

    # preprocessing:
    # order (bimol) reactants and products in the same fashion
    # example if you have A+B and B+A they will be ordered in the same way
    pes_rct_lst = order_rct_bystoich(pes_rct_lst)
    pes_prd_lst = order_rct_bystoich(pes_prd_lst)

    # put everything in a dataframe. indices are the numbers given by enumerate
    len_rct_prd = numpy.array(list(map(len, pes_rct_lst))) + \
        numpy.array(list(map(len, pes_prd_lst)))
    pes_df = pd.DataFrame(
        numpy.array([pes_rct_lst, pes_prd_lst, len_rct_prd], dtype=object).T,
        index=numpy.arange(0, len(pes_rxn_name_lst)),
        columns=['rcts', 'prds', 'N_rcts_prds'])
    # order by total number of species (N of reactants + N of products)
    pes_df = pes_df.sort_values(by=['N_rcts_prds', 'rcts', 'prds'])
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
                    # This works for unimol species
                    # Need also to verify bimol wellskipping channels
                    if conn_chnls_idx not in connected_to:
                        connected_to.append(conn_chnls_idx)
                elif (len(spc_pair) == 1 and
                      spc_pair[::-1] in conndct[conn_chnls_idx]):
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
def print_pes_channels(pes_dct):
    """ Print the PES
    """

    for (form, pidx, sidx), chnls in pes_dct.items():
        pes_str = 'PES: {},'.format(pidx+1)
        sub_pes_str = 'SUB-PES: {},'.format(sidx+1)
        chnl_str = 'Channels: ' + ','.join(str(chnl[0]+1) for chnl in chnls)
        print('!', form, pes_str, sub_pes_str, chnl_str)
        for chnl in chnls:
            _, rxn = chnl
            print('  {} = {}'.format(
                ' + '.join(rxn[0]),
                ' + '.join(rxn[1])))
