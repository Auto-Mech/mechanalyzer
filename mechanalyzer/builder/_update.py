""" Functions which handle updating the master and temporary
    objects describing the reactions and species of a mechanism

    Need to add code to remove unstable species as reactants
    (going to bimol prods?)
"""

import sys
import itertools
import automol
from autoreact.params import RxnParams
import thermfit
from mechanalyzer.builder._names import rxn_ich_to_name
from mechanalyzer.builder._names import ich_name_dct
from mechanalyzer.builder._names import functional_group_name


# Handles Species Object Updates
def update_spc_dct_from_reactions(rxns, spc_dct):
    """ Update a species with species from a set of reactions
    """

    spc_lst = _spc_from_reactions(rxns)
    spc_dct = update_spc_dct(spc_lst, spc_dct)

    return spc_dct


def update_spc_dct(spc_ichs, spc_dct):
    """ Update the species dictionary with a list of species
    """

    print('\nAdding new unique species to mechanism by',
          'adding to mechanism spc_dct...\n')

    _ich_name_dct = ich_name_dct(spc_dct)

    # Add species dict to mech dct if it is not already in mechanism
    # Build a lst of species that have been added to the mechanism
    for ich in spc_ichs:
        if ich not in _ich_name_dct:
            # Generate a functional group name
            name = functional_group_name(ich, name='')

            if name in spc_dct:
                print(f'WARNING: GENERAED NAME {name} ALREADY IN DCT!!!')
                sys.exit()

            # Generate the data dct
            rgt_dct = thermfit.create_spec(ich)

            # Add to the overall mechanism spc_dct and new species lst
            smi = automol.inchi.smiles(ich)
            print(f'Adding species {name} = {smi} = {ich}')

            spc_dct.update({name: rgt_dct})

    return spc_dct


# Handles Reaction Object Updates
def update_rxn_dct(rxn_lst, rxn_dct, spc_dct):
    """ Update the reaction dictionary with a list of reactions
    """

    print('\nAdding new unique reactions to mechanism...\n')

    rxn_dct = rxn_dct if rxn_dct is not None else {}
    for rxn in rxn_lst:
        rxn_wname = rxn_ich_to_name(rxn, spc_dct)
        if _unique_reaction(rxn_wname, rxn_dct):

            # Convert to names and print message
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


def remove_improper_reactions(rxn_param_dct, mech_spc_dct,
                              stereo=True, reverse=True):
    """ Remove reactions from the mechanism that do not correspond
        to proper, physical elementary step reactions.
    """
    print('call improper')
    ste_rxn_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        print(f'Checking {rxn[0]}->{rxn[1]}')
        rcts_ich = tuple(mech_spc_dct[rct]['inchi'] for rct in rxn[0])
        prds_ich = tuple(mech_spc_dct[prd]['inchi'] for prd in rxn[1])

        rxn_obj_sets = automol.reac.rxn_objs_from_inchi(
            rcts_ich, prds_ich, stereo=stereo)
        if rxn_obj_sets is not None:
            rxn_class = rxn_obj_sets[0][0].class_
            print(f' - Keep: IDd {rxn_class} for reaction {rxn[0]}->{rxn[1]}')
            ste_rxn_param_dct[rxn] = params
        else:
            # Check if the reverse reaction cannot be ID'd
            if reverse:
                rxn_obj_sets = automol.reac.rxn_objs_from_inchi(
                    prds_ich, rcts_ich, stereo=stereo)
                if rxn_obj_sets is not None:
                    rev_rxn = (rxn[1], rxn[0], rxn[2])
                    ste_rxn_param_dct[rev_rxn] = params
                    print(
                        f' - Keep: (Reverse) IDd {rxn_class} for '
                        f'reaction {rxn[0]}->{rxn[1]}')
                else:
                    # print(f'Removing reaction {rcts_ich}->{prds_ich}')
                    print(
                        ' - Remove: No ID in either direction for '
                        f'reaction {rxn[0]}->{rxn[1]}')
            else:
                # print(f'Removing reaction {rcts_ich}->{prds_ich}')
                print(f' - Remove: No ID for reaction {rxn[0]}->{rxn[1]}')

    return ste_rxn_param_dct


def remove_unstable_reactions(rxn_param_dct, mech_spc_dct):
    """ Remove reaction A -> B + C where A is an unstable reactant.
        This is because these reactions are readded with other codes.
    """
    _rxn_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        reacs, prods, _ = rxn
        if len(reacs) == 1 and len(prods) == 2:
            rct_geo = automol.inchi.geometry(mech_spc_dct[reacs[0]]['inchi'])
            rct_zma = automol.geom.zmatrix(rct_geo)
            instab_zmas = automol.reac.instability_product_zmas(rct_zma)
            if instab_zmas:
                print('Removing reaction', reacs, prods)
                _rxn_param_dct[rxn] = params
            else:
                _rxn_param_dct[rxn] = params
        else:
            _rxn_param_dct[rxn] = params

    return _rxn_param_dct


def _unique_reaction(rxn, rxn_dct):
    """ Determine if a reaction is in the parameter_dictionary

        does not deal with the third body, so function does not work
    """
    return not any(rxn in rxn_dct for rxn in _make_reaction_permutations(rxn))


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
