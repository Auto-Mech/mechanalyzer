def stoich_gra(gra):
    """ Stoichiometry dictionary from a graph
    """

    atms = automol.graph.atoms(gra)
    stoich_dct = {}
    hcount = 0
    for atm in atms:
        hcount += np.floor(atms[atm][1])
        if atms[atm][0] in stoich_dct:
            stoich_dct[atms[atm][0]] += 1
        else:
            stoich_dct[atms[atm][0]] = 1
    if 'H' not in stoich_dct:
        stoich_dct['H'] = hcount
    else:
        stoich_dct['H'] += hcount

    return stoich_dct


def stoich(ich):
    """
    Finds the stoichiometry of a molecule
    INPUT:
    ich  -- STR inchi
    OUTPUT:
    stoich -- dictionary with key = STR atomsymbol,
                val = INT number of atomsymbol in molecule

    replace with automol.inchi.formula
    """

    stoich_dct = {'H': 0}
    gra = automol.inchi.graph(ich)
    atms = automol.graph.atoms(gra)
    for atm in atms:
        stoich_dct['H'] += atms[atm][1]
        if atms[atm][0] in stoich_dct:
            stoich_dct[atms[atm][0]] += 1
        else:
            stoich_dct[atms[atm][0]] = 1
    return stoich_dct
