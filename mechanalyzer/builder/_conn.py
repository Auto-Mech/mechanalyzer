""" Determine how various PESs may connect

    Let PES1 and PES2 be two independent potential energy surfaces.

    Let W is a well on PES1, and P1 is a bimol reaction on PES2

    If W is a part of the reaction P1, then PES1 and PES2 can be
    kinetically connected via non-thermal processeses.

    Fails to link multiple PESs with dif wells

    Linked:
    R+O2 = QOOH
    QOOH + O2 = OOQOOH

    Not Linked to Above I think (need new code)
    OOQOOH = HOOQOOH

    above will come back in two conn dcts
"""


def connected_surfaces(pes_dct, excl_spc=()):
    """ Determine how PESs are connected.
    """
    return _find_spc(pes_dct, excl_spc=excl_spc)


def _find_spc(pes_dct, excl_spc=()):
    """ New version of connector code
    """

    # Initial dict of spc[name] = ((pes_id, rxn1), (pes_id2, rxn2), ...)
    spcdct = {}
    for pes_id, pes_rxns in pes_dct.items():
        for rxn in pes_rxns:
            rcts, prds = rxn[0], rxn[1]
            for rct in rcts:
                if rct in spcdct:
                    spcdct[rct] += ((pes_id, rxn),)
                else:
                    spcdct[rct] = ((pes_id, rxn),)
            for prd in prds:
                if prd in spcdct:
                    spcdct[prd] += ((pes_id, rxn),)
                else:
                    spcdct[prd] = ((pes_id, rxn),)

    # Build list of PESs that each species exists in
    # Only enumerate species not in excl list and appear in multiple PESs
    spcdct2 = {}
    for name, pes_rxn_lst in spcdct.items():
        if name not in excl_spc:
            _lst = set(tuple(pes_id for pes_id, _ in pes_rxn_lst))
            if len(_lst) > 1:
                spcdct2[name] = _lst

    return spcdct2
