""" Determine how various PESs may connect

    Let PES1 and PES2 be two independent potential energy surfaces.

    Let W is a well on PES1, and P1 is a bimol reaction on PES2

    If W is a part of the reaction P1, then PES1 and PES2 can be
    kinetically connected via non-thermal processeses.
"""

def find_connections(well_lst_dct, spc_lst_dct):
    """ find places where rxns can connect

        dct = {idx: (name1, name2, ..., name3)}
    """

    # Just build a list of wells
    wells = []
    for well_lst in wells_lst_dct.values():
        wells.extend(well_lst)

    wells_in_pes = {}
    for well in wells:
        idx_lst = []
        for pes_idx, spc_lst in spc_lst_dct.items():
            if well in spc_lst:
                idx_lst.append(pes_idx)
        wells_in_pes.update({well: idx_lst})

    # OLD
    for idx1, wells1 in wells_dct.items():
        conns = {} 
        conn_lst = []
        for idx2, spc_lst in spc_dct.items():
            # See if the well in pes1 in spclst of pes2
            conn_wells = wells1 & spc_lst
            if conn_wells:
                conn_lst.append(idx2)

        if conn_lst:
            conns.update({idx1: (conn_wells, conn_lst)})

    return conns


def find_species_and_wells(pes_dct):
    """ Find all the wells on each PES

        :param pes_dct: dictionary of pess
        :type pes_dct: dict[?]
        :rtype: dict[int: tuple]
    """

    # Build idx_name dct {idx: [names]}}
    wells = {}
    spc = {}
    for pes_idx, pes_rxns in enumerate(ped_dct):
        pes_wells = []
        pes_species = []
        # Only look at multichannel PESs for the wells?
        for rxn in pes_rxns:
            rcts, prds = rxn[0], rxn[1]

            # Find the wells
            if len(rcts) == 1:
                pes_wells.append(rcts[0])
            if len(prds) == 1:
                pes_wells.append(prds[0])

            # Get the species for the reaction
            pes_spc.extend(rcts + prds)

        wells.update({ped_idx: set(pes_wells)}) 
        spc.update({ped_idx: set(pes_species)}) 

    # Build name_idx dct for wells {name: idx}
    name_idx_dct = {}
    for idx, names in wells.items():
        for name in names:
            name_idx_dct[name] = idx

    return wells, spc
