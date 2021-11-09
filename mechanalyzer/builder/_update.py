""" Functions which handle updating the master and temporary
    objects describing the reactions and species of a mechanism
"""

import itertools
import automol
from autoreact.params import RxnParams
import thermfit
import mechanalyzer.parser


# Handles Species Object Updates
def update_spc_dct_from_reactions(rxns, spc_dct):
    """ Update a species with species from a set of reactions
    """

    spc_lst = _spc_from_reactions(rxns)
    spc_dct, new_spc_ichs = update_spc_dct(spc_lst, spc_dct)

    return spc_dct, new_spc_ichs


def update_spc_dct(spc_ichs, spc_dct):
    """ Update the species dictionary with a list of species
    """

    print('\nAdding new unique species to mechanism by',
          'adding to mechanism spc_dct...\n')

    # Generate the bookkeeping dictionaries to assign names
    fml_count_dct, ich_name_dct = _make_spc_bookkepping(spc_dct)

    # Add species dict to mech dct if it is not already in mechanism
    # Build a lst of species that have been added to the mechanism
    new_spc_ichs = ()
    for ich in spc_ichs:
        if ich not in ich_name_dct:
            # Generate unique name with formula, update formula dct
            fml = automol.inchi.formula_string(ich)
            name, fml_count_dct = mechanalyzer.parser.spc.assign_unique_name(
                fml, fml_count_dct, spc_dct)

            # Generate the data dct
            rgt_dct = thermfit.create_spec(ich)

            # Add to the overall mechanism spc_dct and new species lst
            smi = automol.inchi.smiles(ich)
            print(f'Adding species {name} = {smi} = {ich}')

            spc_dct.update({name: rgt_dct})
            new_spc_ichs += (ich,)

    return spc_dct, new_spc_ichs


# Handles Reaction Object Updates
def update_rxn_dct(rxn_lst, rxn_dct, spc_dct):
    """ Update the reaction dictionary with a list of reactions
    """

    print('\nAdding new unique reactions to mechanism...\n')

    for rxn in rxn_lst:
        if _unique_reaction(rxn, rxn_dct):

            # Convert to names and print message
            rxn_wname = _rxn_ich_to_name(rxn, spc_dct)
            print(f'Adding reaction {rxn_wname} to param dct')

            rxn_dct[rxn_wname] = RxnParams(
                arr_dct={'arr_tuples': ((1.0, 0.0, 0.0),)})

    return rxn_dct


# Removal functions
def remove_spc_not_in_reactions(rxn_param_dct, mech_spc_dct):
    """ Remove species from the spc dct not currently
        in the list of reactions in the rxn_dct
    """

    # spc_in_rxns = _spc_from_reactions(rxns)
    spc_in_rxns = ()
    for rxn in rxn_param_dct:
        spc_in_rxns += rxn[0]
        spc_in_rxns += rxn[1]
    spc_in_rxns = set(spc_in_rxns)

    new_mech_spc_dct = {}
    for name, dct in mech_spc_dct.items():
        if name in spc_in_rxns:
            new_mech_spc_dct[name] = dct
        elif any(x in name for x in ('cbh0_', 'cbh1_', 'cbh2_', 'cbh_3')):
            new_mech_spc_dct[name] = dct
        else:
            print(f'Remove species: {name}')

    return new_mech_spc_dct


def remove_improper_reactions(rxn_param_dct, mech_spc_dct, stereo=True):
    """ Remove reactions from the mechanism that do not correspond
        to proper, physical elementary step reactions.
    """

    ste_rxn_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        rcts_ich = tuple(mech_spc_dct[rct]['inchi'] for rct in rxn[0])
        prds_ich = tuple(mech_spc_dct[prd]['inchi'] for prd in rxn[1])

        rxn_obj_sets = automol.reac.rxn_objs_from_inchi(
            rcts_ich, prds_ich, stereo=stereo)
        if rxn_obj_sets is not None:
            ste_rxn_param_dct[rxn] = params
        else:
            print(f'Removing reaction {rcts_ich}->{prds_ich}')

    return ste_rxn_param_dct


def _unique_reaction(rxn, rxn_dct):
    """ Determine if a reaction is in the parameter_dictionary

        does not deal with the third body, so function does not work
    """
    rxns = _make_reaction_permutations(rxn)
    unique = not any(rxn in rxn_dct for rxn in rxns)
    return unique


# Make the bookkeeping objects
def _make_spc_bookkepping(spc_dct):
    """ make dictionaries that maintain info that is used for
        bookkeeping dictionaries
    """

    fml_count_dct = mechanalyzer.parser.spc.formula_count_dct(spc_dct)
    ich_name_dct = _ich_name_dct(spc_dct)

    return fml_count_dct, ich_name_dct


# Other helper functions
def _spc_from_reactions(rxns):
    """ Build a species dictionary from a list of reactions
        which define a reaction using the inchi strings
    """
    rgts = ()
    for rxn in rxns:
        rgts += tuple(itertools.chain(*rxn))
    uni_rgts = tuple(
        rgt for rgt in
        automol.util.remove_duplicates_with_order(rgts)
        if rgt not in (None, '+M', '(+M)'))
    return uni_rgts


def _make_reaction_permutations(rxn):
    """ reactions
    """
    # Get all reactions R=P including all permutations of R and P
    rct_perms = tuple(itertools.permutations(rxn[0]))
    prd_perms = tuple(itertools.permutations(rxn[1]))

    all_rxns = (
        tuple(itertools.product(rct_perms, prd_perms)) +
        tuple(itertools.product(prd_perms, rct_perms))
    )

    # Remove duplicates
    all_rxns = tuple(set(all_rxns))

    # Re-add the third body
    third_body = rxn[2]
    all_rxns = tuple((*rxn, third_body) for rxn in all_rxns)

    return all_rxns


def _rxn_ich_to_name(rxn, spc_dct):
    """ Set a reaction described by inchis to one where it is
        described by mechanism names
    """

    ich_name_dct = _ich_name_dct(spc_dct)
    rxn2 = (
        tuple(ich_name_dct[rgt] for rgt in rxn[0]),
        tuple(ich_name_dct[rgt] for rgt in rxn[1]),
        (None,)
    )

    return rxn2


def rxn_name_str(rxn, newline=False):
    """ get a reaction name string
    """
    if newline:
        rstr = ' =\n       '.join((' + '.join(rxn[0]), ' + '.join(rxn[1])))
    else:
        rstr = ' = '.join((' + '.join(rxn[0]), ' + '.join(rxn[1])))

    return rstr


def _ich_name_dct(spc_dct):
    """ get dct[ich] = name
    """

    ich_dct = {}
    for key in spc_dct.keys():
        ich_dct[spc_dct[key]['inchi']] = key

    return ich_dct
