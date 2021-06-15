"""
  Toy code for mechanism building
"""

import itertools
# import functools
import automol


# GENERATE ALL OF THE REACTIONS
def build_mechanism(ini_spc_ichs, rxn_series='gen'):
    """ Use the lst of reactions to build objects to describe mechanism
    """

    # Initialize the master spc lst and gra dcts using the seed spc lst
    master_spc_lst = _initialize_spc(ini_spc_ichs)

    # Initialize empty master rxn lst since we assume no rxns given by user
    master_rxn_lst = _initialize_rxn()

    # Set the reactions series and oop over the series to
    # generate reactions to add to the mechanism
    rseries = _set_reaction_series(rxn_series=rxn_series)
    for series in rseries:

        # Unpack series
        spc_set, rxn_set, _ = series
        print('sets', spc_set, rxn_set)

        # Build all of the gras needed for the new reaction generation
        spc_lst = _set_spc_lst(master_spc_lst, spc_set)
        print('spc lst', spc_lst)

        # Generate reactions for desired reactants and classes
        new_rxns, new_spc = _generate_rxns(spc_lst, rxn_set)

        # Update the master spc and rxn lsts
        master_spc_lst, master_rxn_lst = _update_mech(
            master_spc_lst, master_rxn_lst, new_rxns, new_spc)

    return master_spc_lst, master_rxn_lst


# Build the reactions for some master lst of inchis
def _initialize_spc(ini_spc_ichs):
    """ Build the first master species lst
    """

    master_spc_lst = tuple()
    for spc_ich in ini_spc_ichs:
        print(spc_ich)
        master_spc_lst += ((spc_ich, automol.inchi.graph(spc_ich)),)

    return master_spc_lst


def _initialize_rxn():
    """ Build an initial set of reactions
    """
    return ()


def _set_spc_lst(master_spc_lst, spc_set):
    """ Determine a list of species
    """

    # Build the lst of spc ichs to generate reactions for
    if spc_set == 'all':
        spc_ichs = master_spc_lst
    elif spc_set == 'radicals':
        spc_ichs = radicals(master_spc_lst)

    return spc_ichs


def _generate_rxns(rct_ichs, rxnclasses):
    """ Loop over all of the inchis
    """

    rxns = tuple()

    for ich in rct_ichs:
        found_rxns = tuple()
        for rclass in rxnclasses:
            found_rxns += RXN_DCT[rclass](ich)
        if found_rxns:
            rxns += found_rxns

    return rxns


def _bimol(func, rct1_gra, rct2_gras):
    """ Generate the inchis for the reacs and prods for
        hydrogen abstractions for a given species
    """

    rct_ich = automol.graph.inchi(rct1_gra)
    rad_ich = ''  # wrong

    new_rxns = tuple()
    new_spc = tuple()

    for rct2_gra in rct2_gras:
        for gras in func(rct1_gra, rct2_gra):
            prod_ichs = tuple(map(automol.graph.inchi, gras))

            new_rxns += (((rct_ich, rad_ich), prod_ichs),)
            new_spc += tuple(zip(prod_ichs, gras))

    return new_rxns, new_spc


def _unimol(func, rct_gra):
    """ Generate the inchis for the reacs and prods for
        hydrogen abstractions for a given species
    """

    rct_ich = automol.graph.inchi(rct_gra)

    new_rxns = tuple()
    new_spc = tuple()

    for gras in func(rct_gra):

        prod_ichs = tuple(map(automol.graph.inchi, gras))

        new_rxns += (((rct_ich,), prod_ichs),)
        new_spc += tuple(zip(prod_ichs, gras))

    return new_rxns, new_spc


def _unimol_to_unimol(func, rct_gra):
    """ Generate the inchis for the reacs and prods for
        hydrogen abstractions for a given species
    """

    rct_ich = automol.graph.inchi(rct_gra)

    new_rxns = tuple()
    new_spc = tuple()

    for gra in func(rct_gra):

        prd_ich = automol.graph.inchi(gra)

        new_rxns += ((rct_ich,), (prd_ich,))
        new_spc += ((prd_ich, gra))

    return new_rxns, new_spc


RXN_DCT = {
    # 'hydrogen_abstractions': functools.partial(
    #     _bimol, automol.reac.prod_hydrogen_abstraction),
    # 'homolytic_scissions': functools.partial(
    #     _unimol, automol.reac.prod_homolytic_scission),
    # 'beta_scissions': functools.partial(
    #     _unimol, automol.reac.prod_beta_scission),
    # 'hydrogen_migrations': functools.partial(
    #     _unimol_to_unimol, automol.reac.prod_hydrogen_migration)
}


def _set_reaction_series(rxn_series='gen'):
    """ Based on a desired mechanism, set the series of steps to build reactions

        Now, we assume just reactions, no seed species; those given by user
        Impose restrictions on bimol rcts?
    """

    if rxn_series == 'gen':
        rseries = (
            ('all', ('hydrogen_abstractions', 'homo_scissions'), {}),
            ('radicals', ('hydrogen_migrations', 'beta_scissions'), {})
        )

    return rseries


def _update_mech(master_spc_tup, master_rxn_tup, new_rxn_tup, new_spc_tup):
    """ Update the mechanism master lsts

        master_spc_lst: ( (ich1, gra1), (ich2, gra2), ... )
        master_rxn_lst:
           ( ('FORM', (rct_ich1, rct_ich2), (prd_ich1, prd_ich2)), ... )
        new_rxn_lst: ( ((rct_ich1, rct_ich2), (prd_ich1, prd_ich2)), ... )
    """

    # maybe just work with inchis and get around the graphs
    # Grab all of the species from new rxns and add ones not in master tup
    new_rxns_spc_tup = tuple(itertools.chain(*new_spc_tup))  # wrong
    new_spc_tup = tuple(set(new_rxns_spc_tup).difference(master_spc_tup))
    master_spc_tup += new_spc_tup

    # Do same for reactions
    new_spc_tup = tuple(set(new_rxn_tup).difference(master_rxn_tup))
    master_rxn_tup += new_spc_tup

    return master_spc_tup, master_rxn_tup


# Various helper functions for building reactions
def _ich_to_gra(spc_ichs, gra_dct):
    """ convert lst of ichs to gras
    """
    spc_gras = tuple(gra_dct[ich] for ich in spc_ichs)
    for ich in spc_ichs:
        spc_gras += (gra_dct[ich],)

    return spc_gras


def unique_ichs_in_rxns(rxns):
    """ Determine the inchis that are unique in lst of reactions
    """

    uni_rxn_ichs = set()
    for rxn in rxns:
        uni_rxn_ichs = uni_rxn_ichs.union(set(itertools.chain(*rxn)))

    return tuple(uni_rxn_ichs)


def radicals(ich_lst):
    """ Determine the radicals
    """

    rad_ich_lst = tuple()
    for ich, gra in ich_lst:
        if automol.graph.radical_species(gra):
            rad_ich_lst += ((ich, gra),)

    return rad_ich_lst


def _combine_bimol(rct1_ichs, rct2_ichs):
    return tuple(itertools.product(rct1_ichs, rct2_ichs))


# Handle building a chemkin mechanism ands species file objets
def mechanism_strs(rxns, ich_name_dct):
    """ Build strings
    """

    # Build the objects for the mechanism and species files
    ich_name_dct = build_spc_dct(rxns)
    mechanism = build_mech_dat(rxns, ich_name_dct)

    mech_str = build_mech_str(mechanism)
    spc_str = build_spc_str(ich_name_dct)

    return mech_str, spc_str


def build_spc_dct(rxns):
    """ Check if the ichs currently exist in the species dict, if not, add them
    """

    ich_name_dct = {}

    # Determine all of the unique inchi strings from the reactions
    uni_rxn_ichs = unique_ichs_in_rxns(rxns)

    # Add the ichs to a dct indexed by formula
    fml_dct = {}
    for ich in uni_rxn_ichs:
        fml = automol.formula.string2(automol.inchi.formula(ich))
        if fml in fml_dct:
            fml_dct[fml].append(ich)
        else:
            fml_dct[fml] = [ich]

    # Now build final dct where the mechanism names use the formula
    ich_name_dct = {}
    for fml, ich_lst in fml_dct.items():
        for idx, ich in enumerate(ich_lst):
            name = '{}({})'.format(fml, str(idx+1))
            ich_name_dct[ich] = name

    return ich_name_dct


def build_mech_dat(rxns, ich_name_dct):
    """ build mechanism file
    """

    mechanism = tuple()
    for rxn in rxns:
        rcts, prds = rxn
        rct_names = tuple(ich_name_dct[rct] for rct in rcts)
        prd_names = tuple(ich_name_dct[prd] for prd in prds)
        mechanism += ((rct_names, prd_names),)

    return mechanism


def build_spc_str(spc_dct):
    """ str
    """

    spc_str = ''
    for ich, name in spc_dct.items():
        smi = automol.inchi.smiles(ich)
        spc_str += '{},{},{}\n'.format(name, smi, ich)

    return spc_str


def build_mech_str(mech):
    """ str
    """

    mech_str = ''
    for rxn in mech:
        rcts, prds = rxn
        rct_str = '+'.join(rcts)
        prd_str = '+'.join(prds)
        mech_str += '{}={}\n'.format(rct_str, prd_str)

    return mech_str
