"""
Read the mechanism file
"""

from mechanalyzer.parser import ckin_ as ckin


MECH_INP = 'inp/mechanism.dat'


# Read the mechanism file
def mechanism_file(mech_str, mech_type, spc_dct, sort_rxns=False):
    """ Get the reactions and species from the mechanism input
    """

    # Parse the info from the chemkin file
    if mech_type == 'chemkin':
        formulas, rct_names, prd_names, rxn_names = ckin.parse(
            mech_str, spc_dct, sort_rxns)
    else:
        raise NotImplementedError

    return [formulas, rct_names, prd_names, rxn_names]


# Write a new mechanism file
def write_mechanism_file(pes_dct, path, outname):
    """ Write the mechanism file from a mech dct
    """

    mech_str = ''

    for pes_idx, formula in enumerate(pes_dct):
        print('! PES:', pes_idx+1, formula)
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        for chn_idx, _ in enumerate(pes_rxn_name_lst):
            mech_str += (
                '  {} = {}   1.0 0.0 0.0'.format(
                ' + '.join(pes_rct_names_lst[chn_idx]),
                ' + '.join(pes_prd_names_lst[chn_idx]))
            )


    # Write the file
    mech_file = os.path.join(path, outname)
    with open(mech_file, 'w') as file_obj:
        file_obj.write(spc_str)


# FUNTIONS FOR THE PES DICT OBJECTS CONTAINING INFO FOR THE REACTIONS ON PES
def build_sub_pes_dct_w_rxn_lsts(pes_dct, idx_dct, form_dct,
                                 conn_chnls_dct, run_obj_dct):
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
        run_chnls = []
        print('run dct', run_obj_dct)
        for pes_chn_pair in run_obj_dct:
            pes_num, chn_num = pes_chn_pair
            if idx_dct[pes_num] == formula:
                run_chnls.append(chn_num)

        # Select names from the names list corresponding to chnls to run
        # conn_chnls_dct[formula] = {sub_pes_idx: [channel_idxs]}
        for sub_pes_idx, sub_chnl_idxs in conn_chnls_dct[formula].items():
            rct_names_lst = []
            prd_names_lst = []
            rxn_name_lst = []
            rxn_model_lst = []
            rxn_chn_idxs = []
            for chn_idx in run_chnls:
                if chn_idx-1 in sub_chnl_idxs:
                    rct_names_lst.append(pes_rct_names_lst[chn_idx-1])
                    prd_names_lst.append(pes_prd_names_lst[chn_idx-1])
                    rxn_name_lst.append(pes_rxn_name_lst[chn_idx-1])
                    rxn_model_lst.append(run_obj_dct[(pes_idx, chn_idx)])
                    rxn_chn_idxs.append(chn_idx)

            # Form reaction list (is empty if no chnls requested on sub pes)
            rxn_lst = format_run_rxn_lst(
                rct_names_lst, prd_names_lst, rxn_model_lst, rxn_chn_idxs)

            # Add the rxn lst to the pes dictionary if there is anythin
            if rxn_lst:
                run_pes_dct[(formula, pes_idx, sub_pes_idx+1)] = rxn_lst

    return run_pes_dct


def format_run_rxn_lst(rct_names_lst, prd_names_lst,
                       rxn_model_lst, rxn_chn_idxs):
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
             'model': rxn_model_lst[idx],
             'chn_idx': rxn_chn_idxs[idx],
             'dummy': []}
        )

    return run_lst


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


# FUNCTIONS FOR THE CHANNELS DICT OBJECTS
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
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']

        # Split up channels into a connected sub-pes within a formula
        subpes_idx = 0
        conndct = {}
        connchnls = {}
        for chnl_idx, _ in enumerate(pes_rxn_name_lst):
            connected_to = []
            chnl_species = [list(pes_rct_names_lst[chnl_idx]),
                            list(pes_prd_names_lst[chnl_idx])]
            for conn_chnls_idx in conndct:
                for spc_pair in chnl_species:
                    if len(spc_pair) == 1:
                        if spc_pair in conndct[conn_chnls_idx]:
                            if conn_chnls_idx not in connected_to:
                                connected_to.append(conn_chnls_idx)
                        elif spc_pair[::-1] in conndct[conn_chnls_idx]:
                            if conn_chnls_idx not in connected_to:
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

        # Add connected channels list to the dictionary
        conn_chn_dct[formula] = connchnls

    return conn_chn_dct


# OLD FUNCTION FOR TERMOL RXNS
def conv_termol_to_bimol(rct_zmas, prd_zmas):
    """ Convert termolecular reaction to a bimolecular reaction
    """
    # Force a trimolecular reaction to behave like a bimolecular.
    # Termolecular spc often arise from the direct decomp of some init product.
    # Need to be able to find the TS for channel preceding direct decomp.
    rct_tors_names = []
    if len(rct_zmas) > 2:
        ret = automol.zmatrix.ts.addition(rct_zmas[1:-1], [prd_zmas[-1]])
        new_zma, dist_name, rct_tors_names = ret
        new_zma = automol.zmatrix.standard_form(new_zma)
        babs2 = automol.zmatrix.get_babs2(new_zma, dist_name)
        new_zma = automol.zmatrix.set_values(
            new_zma, {dist_name: 2.2, babs2: 180. * phycon.DEG2RAD})
        rct_zmas = [rct_zmas[0], new_zma]
    elif len(prd_zmas) > 2:
        ret = automol.zmatrix.ts.addition(prd_zmas[1:-1], [prd_zmas[-1]])
        new_zma, dist_name, rct_tors_names = ret
        new_zma = automol.zmatrix.standard_form(new_zma)
        babs1 = automol.zmatrix.get_babs1(new_zma, dist_name)
        aabs1 = babs1.replace('D', 'A')
        new_zma = automol.zmatrix.set_values(
            new_zma, {dist_name: 2.2, aabs1: 170. * phycon.DEG2RAD})
        prd_zmas = [prd_zmas[0], new_zma]

    return rct_zmas, prd_zmas, rct_tors_names
