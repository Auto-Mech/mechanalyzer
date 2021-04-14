""" Determine how various PESs may connect

    Let PES1 and PES2 be two independent potential energy surfaces.

    Let W is a well on PES1, and P1 is a bimol reaction on PES2

    If W is a part of the reaction P1, then PES1 and PES2 can be
    kinetically connected via non-thermal processeses.
"""


def conn_pes(pes_dct):
    """ Determine how PESs are connected

        Fails to link multiple PESs with dif wells

        Linked:
        R+O2 = QOOH
        QOOH + O2 = OOQOOH

        Not Linked to Above I think (need new code)
        OOQOOH = HOOQOOH

        above will come back in two conn dcts

    """

    well_lst_dct, bimols_lst_dct, spc_lst_dct = _find_spc(pes_dct)
    wells_in_pes = _well_locations(well_lst_dct, spc_lst_dct)

    conn_lst = _sort_connected_pes(pes_dct, wells_in_pes)

    return conn_lst


# Determine where well species and all species exist on PESs
def _find_spc(pes_dct):
    """ Find all the wells on each PES.

        Currently gets wells from all PESs, maybe should have
        a check for getting the multichannel.

        :param pes_dct: dictionary of pess
        :type pes_dct: dict[?]
        :rtype: dict[int: tuple]
    """

    # Build idx_name dct {idx: [names]}}
    wells, bimols, spc = {}, {}, {}
    for pes_idx, pes_rxns in enumerate(pes_dct):
        pes_wells, pes_species, pes_bimols = [], [], []

        for rxn in pes_rxns:
            rcts, prds = rxn[0], rxn[1]

            # Find the wells (could get the bimol)
            if len(rcts) == 1:
                pes_wells.append(rcts[0])
            else:
                pes_bimols.append(rcts)
            if len(prds) == 1:
                pes_wells.append(prds[0])
            else:
                pes_bimols.append(prds)

            # Get the species for the reaction
            pes_species.extend(rcts + prds)

        wells.update({pes_idx: set(pes_wells)})
        bimols.update({pes_idx: set(pes_bimols)})
        spc.update({pes_idx: set(pes_species)})

    # Get names where the wells are bimols

    return wells, bimols, spc


def _well_locations(well_lst_dct, spc_lst_dct):
    """ Taking a list of all the wells (from multichannel PESs)

        We determine what PESs the species that comprise these
        wells exist in.

        first get wells_in_pes = {well: [idx_lst]}

        dct = {name: (idx1, idx2, ..., idx3)}
    """

    # Just build a list of wells
    wells = []
    for well_lst in well_lst_dct.values():
        wells.extend(well_lst)

    # Determine what PESs a well exists inside of
    wells_in_pes = {}
    for well in wells:
        idx_lst = []
        for pes_idx, spc_lst in spc_lst_dct.items():
            if well in spc_lst:
                idx_lst.append(pes_idx)
        wells_in_pes.update({well: idx_lst})

    return wells_in_pes


# Parse, modify, and sort well in PES list to write connections
def _sort_connected_pes(pes_dct, wells_in_pes):
    """ figures out what PESs are connected to do a chain of PESs
    """

    # Build new dictionary of well species that appear in multiple PESs
    msurf_wells = {well: idxs
                   for (well, idxs) in wells_in_pes.items()
                   if len(idxs) > 1}

    # Now figure out how the idxs should be lists
    # For a given well:
    # (1) find pes idx1 where it is a well
    # (2) find pes idx2 where it is in a bimol (where idx1 != idx2)
    # (3) order indices to that idx1 comes first


    return conn_lst
