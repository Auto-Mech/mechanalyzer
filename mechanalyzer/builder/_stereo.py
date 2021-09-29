"""
  Code to expand the mechanism using stereochemistry
"""

import os
import automol
from autorun import execute_function_in_parallel
import mechanalyzer
from mechanalyzer.builder._update import update_spc_dct_from_reactions
from mechanalyzer.builder._update import update_rxn_dct
from mechanalyzer.builder._update import rxn_name_str


# MAIN CALLABLE
def expand_mech_stereo(mech_rxn_dct, mech_spc_dct, nprocs='auto'):
    """ Build list of stereochemistry to reactions

        Currently, we assume that the species in them mech_spc_dct have
        stereochemistry already added to them.
    """

    def _expand(name_ich_dct, rxns, output_queue):
        """ Expand reactons
        """
        srxns = ()
        for rxn in rxns:
            print(f' - Expanding Stereo for Reaction: {rxn_name_str(rxn)}')
            rxn_ich = _rxn_ich(rxn, name_ich_dct)

            # Get a list of reactions with stereochemistry
            ste_rxns_ich_filt = _ste_rxn_lsts(rxn_ich)
            # ste_rxns_ich_filt = _ste_rxn_lsts2(rxn_ich)

            # Filter unneccseary reactions
            # ste_rxns_ich_filt = _remove_unneeded_reactions(ste_rxns_ich)
            print(' - Stereochemical Versions of Reaction: '
                  f'{rxn_name_str(rxn)}')
            for srxn in ste_rxns_ich_filt:
                # print(_rxn_smiles(x))
                print('    ', rxn_name_str(srxn, newline=True))

            # Add to overall stereo reactions list
            srxns += ste_rxns_ich_filt

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

    print(f' - Reaction identified as {rxn_obj.class_}.')

    # Build a list of
    ste_rxn_ichs = ()
    for ste_rxn in automol.reac.expand_stereo(rxn_obj):
        rct_gras = automol.reac.reactant_graphs(ste_rxn)
        prd_gras = automol.reac.product_graphs(ste_rxn)
        rct_ichs = tuple(map(automol.graph.stereo_inchi, rct_gras))
        prd_ichs = tuple(map(automol.graph.stereo_inchi, prd_gras))

        ste_rxn_ichs += ((rct_ichs, prd_ichs, (None,)),)

    return ste_rxn_ichs


def _ste_rxn_lsts2(rxn_ich):
    """ Build reaction onjects
    """

    # Build reaction objects
    rxn_obj_sets = automol.reac.util.rxn_objs_from_inchi(
        rxn_ich[0], rxn_ich[1])
    rxn, _, rct_geos, prd_geos = rxn_obj_sets[0]

    print(f' - Reaction identified as {rxn.class_}')

    srxn = automol.reac.add_stereo_from_geometries(rxn, rct_geos, prd_geos)

    # Build a list of
    ste_rxn_ichs = ()
    for ste_rxn in automol.reac.expand_product_stereo(srxn):
        rct_gras = automol.reac.reactant_graphs(ste_rxn)
        prd_gras = automol.reac.product_graphs(ste_rxn)
        rct_ichs = tuple(map(automol.graph.stereo_inchi, rct_gras))
        prd_ichs = tuple(map(automol.graph.stereo_inchi, prd_gras))

        ste_rxn_ichs += ((rct_ichs, prd_ichs, (None,)),)

    return ste_rxn_ichs


# Functions to check and sort the reactions by stereochemistry
def _remove_unneeded_reactions(srxns):
    """ Take all reactions that occur from stereochemically
        and determine which reactions should be kept and
        whihc are unneccessary

        Keep E<->E, Z<->Z
        Remove RtoR and StoS if only one stereocenter.
    """
    pnums = ()
    _srxns = ()
    for srxn in srxns:
        ste_lyr = automol.inchi.stereo_sublayers(srxn)

        if 'b' not in ste_lyr:
            if 't' in ste_lyr:
                tlyr = ste_lyr.get('t')
                # Only look for single stereo center
                if ',' not in tlyr:
                    num = tlyr.replace('+', '').replace('-', '')
                    pnums += (num,)
                    if num not in pnums:
                        _srxns += srxn

    return _srxns


def _rxn_ich(rxn, ich_dct):
    """ Convert a reacion written with spc names to spc inchis
    """
    return (
        tuple(ich_dct[rgt] for rgt in rxn[0]),
        tuple(ich_dct[rgt] for rgt in rxn[1]),
        (None,)
    )


def _rxn_smiles(rxn):
    """ write a reaction into smles
    """
    return (
        tuple(automol.inchi.smiles(rgt) for rgt in rxn[0]),
        tuple(automol.inchi.smiles(rgt) for rgt in rxn[1]),
        (None,)
    )
