"""
  Code to expand the mechanism using stereochemistry

    maybe remove all inchis that are not incomplete?
"""

import os
import automol
from autorun import execute_function_in_parallel
import mechanalyzer.parser
from mechanalyzer.builder._update import update_spc_dct_from_reactions
from mechanalyzer.builder._update import update_rxn_dct
from mechanalyzer.builder._update import rxn_name_str
from chemkin_io.writer._util import format_rxn_name


# MAIN CALLABLE
# def expand_mech_stereo(mech_rxn_dct, mech_spc_dct, nprocs='auto'):
def expand_mech_stereo(mech_rxn_dct, mech_spc_dct, nprocs=1):
    """ Build list of stereochemistry to reactions

        Currently, we assume that the species in them mech_spc_dct have
        stereochemistry already added to them.
    """

    def _expand(name_ich_dct, rxns, output_queue):
        """ Expand reactons
        """
        srxns = ()
        for rxn in rxns:

            # Split thrdbdy off, not needed for stereo code, add back later
            _rxn = (rxn[0], rxn[1])
            thrdbdy = rxn[2]

            log1 = f'\nExpanding Stereo for Reaction: {format_rxn_name(rxn)}\n'

            # Reformat reaction to use InChI instead of mechanism name
            rxn_ich = _rxn_ich(_rxn, name_ich_dct)

            # Build list of all stereochemically allowed versions of reaction
            ste_rxns_lst, log2 = _ste_rxn_lsts(rxn_ich)
            
            # Filter redundant reactions from each enantiomer pairs
            ste_rxns_lst, removed_rxns_lst = _remove_enantiomer_reactions(
                ste_rxns_lst)

            # Appropriately format the reactions with third body
            ste_rxns_lst = _add_third(ste_rxns_lst, thrdbdy)
            removed_rxns_lst = _add_third(removed_rxns_lst, thrdbdy)

            # Print final list of stereochemical reactions to potentially add
            log3 = _stereo_results(rxn, ste_rxns_lst, removed_rxns_lst)

            # Add to overall stereo reactions list
            srxns += ste_rxns_lst

            # Print log message for reaction
            log = log1 + log2 + log3
            print(log)

        output_queue.put(srxns)
        print(f'Processor {os.getpid()} finished')

    # Dictionaries to map name to inchi
    name_ich_dct = mechanalyzer.parser.spc.name_inchi_dct(mech_spc_dct)

    # Generate all stereo reactons from the initial set
    rxns = tuple(mech_rxn_dct.keys())
    args = (name_ich_dct,)
    ste_rxns = execute_function_in_parallel(_expand, rxns, args, nprocs=nprocs)

    # Update the mechanism objects with unique spc and rxns
    mech_spc_dct, _ = update_spc_dct_from_reactions(ste_rxns, mech_spc_dct)
    mech_rxn_dct = update_rxn_dct(ste_rxns, mech_rxn_dct, mech_spc_dct)

    return mech_rxn_dct, mech_spc_dct


# Build reaction lists
def _ste_rxn_lsts(rxn_ich):
    """ Build reaction onjects
    """

    # Build reaction objects
    rxn_obj_sets = automol.reac.util.rxn_objs_from_inchi(
        rxn_ich[0], rxn_ich[1])
    rxn_obj = rxn_obj_sets[0][0]  # expand just with rxn object

    # Build a list of stereo reactions
    ste_rxn_ichs = ()
    for ste_rxn in automol.reac.expand_stereo(rxn_obj):
        rct_gras = automol.reac.reactant_graphs(ste_rxn)
        prd_gras = automol.reac.product_graphs(ste_rxn)
        rct_ichs = tuple(map(automol.graph.stereo_inchi, rct_gras))
        prd_ichs = tuple(map(automol.graph.stereo_inchi, prd_gras))

        ste_rxn_ichs += ((rct_ichs, prd_ichs),)

    # Set log message
    log = f' - Reaction identified as {rxn_obj.class_}.\n'

    return ste_rxn_ichs, log


# Functions to check and sort the reactions by stereochemistry
def _remove_enantiomer_reactions(ste_rxn_lst):
    """ Take all reactions that occur from stereochemically
        and determine which reactions should be kept and
        whihc are unneccessary
    """

    # Remove redundant sets and rebuild proper list
    f_ste_rxn_lst = automol.inchi.filter_enantiomer_reactions(ste_rxn_lst)

    # Print the removed reactions
    removed_ste_rxn_lst = set(ste_rxn_lst) - set(f_ste_rxn_lst)

    return f_ste_rxn_lst, removed_ste_rxn_lst


# Formatters and printers
def _stereo_results(rxn, f_ste_rxns_lst, removed_ste_rxns_lst):
    """ Print the final filtered reactions and those removed
    """

    # Print final list of reactions
    log = f' - Stereochemical Versions of Reaction: {rxn_name_str(rxn)}\n'
    for ste_rxn in f_ste_rxns_lst:
        log += '    ' + rxn_name_str(ste_rxn, newline=True) + '\n'

    # Print removed reactions
    if removed_ste_rxns_lst:
        log += (' - Redundant, enantiomeric reactions '
                'precluded from final list\n')
        for ste_rxn in removed_ste_rxns_lst:
            log += '    ' + rxn_name_str(ste_rxn, newline=True) + '\n'

    return log


def _add_third(rxn_lst, thrdbdy):
    """ Format a rxn list to have the third-body added back
    """
    return tuple((rxn[0], rxn[1], thrdbdy) for rxn in rxn_lst)


def _rxn_ich(rxn, ich_dct):
    """ Convert a reacion written with spc names to spc inchis
        Third body list remains the same
    """
    return (
        tuple(ich_dct[rgt] for rgt in rxn[0]),
        tuple(ich_dct[rgt] for rgt in rxn[1]),
    )


def _rxn_smiles(rxn):
    """ write a reaction into smles
    """
    return (
        tuple(automol.inchi.smiles(rgt) for rgt in rxn[0]),
        tuple(automol.inchi.smiles(rgt) for rgt in rxn[1]),
    )
