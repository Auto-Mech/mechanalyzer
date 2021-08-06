"""
  Code to expand the mechanism using stereochemistry
"""

import automol
import mechanalyzer


# MAIN: Reaction Stereochem
def expand_mech_stereo(rxn_dct, spc_dct):
    """ Build list of stereochemistry to reactions
    """

    # Dictionaries to map inchi <-> name
    name_ich_dct = mechanalyzer.parser.spc.name_inchi_dct(spc_dct)
    ich_name_dct = automol.util.dict_.invert(name_ich_dct)

    # Get lists of the reactions
    rxns = list(rxn_dct.keys())
    rxns_ich = _rxn_lst_by_ich(rxns, name_ich_dct)

    # Get a list of reactions with stereochemistry
    ste_rxns_ich = _ste_rxn_lsts(rxns_ich)

    # Use the reactions and convert them to names and build new spc dct
    _ = _build_new(ste_rxns_ich, spc_dct, ich_name_dct)

    return '\n\nEND'


# Build reaction lists
def _rxn_lst_by_ich(rxns, ich_dct):
    """ take a lst of reactions and convert them to inchis
    """

    rxns_ich = []
    for rxn in rxns:
        rxns_ich.append(
            [[ich_dct[rgt] for rgt in rgts]
             for rgts in rxn]
        )
    print('reg rxns')
    for i, r_rxn in enumerate(rxns_ich):
        print('rxn', i)
        for rgts in r_rxn:
            print(rgts)

    return rxns_ich


def _ste_rxn_lsts(rxns_ich):
    """ Build reaction onjects
    """

    # Obtain rxn objs with and without ste
    ste_rxns = []
    for rxn in rxns_ich:

        # Build reaction objects
        rxn_obj_sets = automol.reac.util.rxn_objs_from_inchi(
            rxn[0], rxn[1])
        rxn_obj = rxn_obj_sets[0][0]  # expand just with rxn object
        # rxn_obj = rxn_obj_sets[0]  # expand with rxn obj and geos

        # Use reaction object to get a list of stereo-inclusive reactions
        ste_rxns.append(_add_stereo_to_rxn(rxn_obj))

    for i, ste_rxn in enumerate(ste_rxns):
        print('rxn', i)
        for rgts in ste_rxn:
            print(rgts)

    return ste_rxns


# Add the stereochemistry to a reaction using automol
def _add_stereo_to_rxn(rxn_obj):
    """ Take a mechanism with reactions and species and
        expand it to include stereoselective reactions.
    """

    # Expand the reaction to get the stereochemistry
    srxns = automol.reac.expand_stereo(rxn_obj)
    # rxn, rct_geos, prd_geos = rxn_obj[0], rxn_obj[2], rxn_obj[3]
    # srxns = automol.reac.add_stereo_from_geometries(rxn, rct_geos, prd_geos)

    # Build a list of reactions
    ste_rxn_ichs = []
    for srxn in srxns:
        ste_rxn_ichs.append(_get_rxn_ichs(srxn))

    return ste_rxn_ichs


def _get_rxn_ichs(rxn):
    """ Get the reactant and product inchis using the rxn objecct
        Probably move to reac.util?
    """
    rct_gras = automol.reac.reactant_graphs(rxn)
    prd_gras = automol.reac.product_graphs(rxn)
    rct_ichs = list(map(automol.graph.stereo_inchi, rct_gras))
    prd_ichs = list(map(automol.graph.stereo_inchi, prd_gras))

    return (rct_ichs, prd_ichs)


# Functions to check and sort the reactions by stereochemistry
# def _stereo_changes(rxn):
#     """ Discern what is changing in the stereochemisry
#
#         How many changes in bond parity (E<->Z)?
#         How many changes in stereo parity (R<->S)?
#     """
# def _keep_reactions_from_set(rxn):
#     """ Take all reactions that occur from stereochemically
#         expanding the reaction
#         and determine which reactions should be kept and
#         whihc are unneccessary
#
#         Keep E<->E, Z<->Z
#         Remove RtoR and StoS if only one stereocenter.
#     """


# Functions to build new rxn and spc lists with mechanism names
def _build_new(ste_rxn_set_lst, spc_dct, ich_dct):
    """ Rxn list should probably be trimmed to the physical reactions

        maybe generate the names alongside the ste rxn generation?
    """

    # new_rxn_lst = []
    new_spc_dct = {}

    # print('ichdct\n', ich_dct)
    # print('spcdct\n', spc_dct)

    for ste_rxn_set in ste_rxn_set_lst:
        # Just looking at the reactants for now
        for rxn in ste_rxn_set:
            rct_ichs = rxn[0]
            print('rct_ichs', rct_ichs)
            for ich in rct_ichs:
                # Get the name associated with the ich,
                # either in dct or generated
                # figure out what the entry in the expanded spc dct should be
                if ich in ich_dct:
                    name = ich_dct[ich]
                    new_spc_dct[name] = spc_dct[name]
                else:
                    # Need a way to get ich dct with stereo removed to
                    # look up ich in non_ste ich
                    name = _generate_name

    return new_spc_dct


def _generate_name():
    """ generate a name
    """
    return 'random'
