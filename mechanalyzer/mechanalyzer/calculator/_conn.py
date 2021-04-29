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
    msurf_wells = _well_locations(well_lst_dct, spc_lst_dct)
    print('msurf_wells', msurf_wells)

    conn_lst = _sort_connected_pes(msurf_wells, well_lst_dct, bimols_lst_dct)

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
    for pes_form, pes_rxns in pes_dct.items():
        pes_wells, pes_species, pes_bimols = (), (), ()

        for rxn in pes_rxns:
            rcts, prds = rxn[0], rxn[1]

            # Find the wells (could get the bimol)
            if len(rcts) == 1:
                pes_wells += (rcts[0],)
            else:
                pes_bimols += (rcts,)
            if len(prds) == 1:
                pes_wells += (prds[0],)
            else:
                pes_bimols += (prds,)

            # Get the species for the reaction
            pes_species += (rcts + prds)

        wells.update({pes_form: set(pes_wells)})
        bimols.update({pes_form: set(pes_bimols)})
        spc.update({pes_form: set(pes_species)})

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
        form_lst = []
        for pes_form, spc_lst in spc_lst_dct.items():
            if well in spc_lst:
                form_lst.append(pes_form)
        wells_in_pes.update({well: form_lst})

    # Trim to wells on multiple surfaces
    msurf_wells = {well: forms
                   for (well, forms) in wells_in_pes.items()
                   if len(forms) > 1}

    return msurf_wells


# Parse, modify, and sort well in PES list to write connections
def _sort_connected_pes(msurf_wells, wells, bimols):
    """ figures out what PESs are connected to do a chain of PESs

        start with the PES where the well-species is an actual well
        then have all the bimolecular species
    """

    msurf_well_wtyp = {}
    for well, surfs in msurf_wells.items():
        for surf in surfs:
            surf_wells = wells[surf]
            if well in surf_wells:
                loc = 'well'
            else:
                surf_bimols = bimols[surf]
                if well in surf_bimols:
                    loc = 'bimol'
            locs.append((well,loc))
    
    # Now figure out how the idxs should be lists
    # For a given well:
    # (1) find pes idx1 where it is a well
    # (2) find pes idx2 where it is in a bimol (where idx1 != idx2)
    # (3) order indices to that idx1 comes first

    return conn_lst


def _check_conns(msurf_wells):
    """ For wells in the dct see if they are connected
    """
    return NotImplementedError
