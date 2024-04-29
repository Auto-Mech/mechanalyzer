""" Generate sets of reactions from reactants to construct mechanisms
"""

import itertools
import automol
from autorun import execute_function_in_parallel
from mechanalyzer.builder._update import update_spc_dct_from_reactions
from mechanalyzer.builder._update import update_rxn_dct
from mechanalyzer.builder._stereo import _add_third
from mechanalyzer.builder._stereo import _ste_rxn_lsts
from mechanalyzer.builder._stereo import _remove_enantiomer_reactions
from mechanalyzer.builder._stereo import _stereo_results
from chemkin_io.writer._util import format_rxn_name


# MAIN CALLABLE FUNCTIONS FOR GENERATING REACTION LISTS
def build_mechanism(mech_spc_dct, mech_rxn_dct, rxn_series, stereo=False, nprocs=1):
    """ Use the lst of reactions to build objects to describe mechanism
    """

    # Initialize new_spc (needed for sequential steps
    # new_spc_lst = None
    print('---------------------------------------------------------\n')

    # Loop over the reaction series consisting of reactants and reaction type
    for sidx, series in enumerate(rxn_series):
        rct1_set, rct2_set, allowed_prds, rxn_typs = series

        print(f'Generating Reactions for Series {sidx+1}')

        # Generate the reactions from the reactants
        ini_rxns = ()
        for rtyp in rxn_typs:

            # Determine reactants to generate reactions for
            rct_ichs, rct_names = _determine_reactants(
                mech_spc_dct, rct1_set, rct2_set, rtyp)
            allowed_prd_ichs = _determine_allowed_products(
                mech_spc_dct, allowed_prds)

            # Use SMILES to print info message for what reactions generating
            for ichs, names in zip(rct_ichs, rct_names):
                _rct_smis = tuple(map(automol.chi.smiles, ichs))
                print(f'{names}: {_rct_smis}')
            print('')
            if allowed_prd_ichs:
                allow_str = ' '.join(allowed_prd_ichs)
                print(f'\nEach reaction must produce {allow_str}')
                print('')

            # Generate Reactions
            # for ichs in rct_ichs:
            #    ini_rxns += generate_reactions(ichs, allowed_prd_ichs, rtyp)
            ini_rxns += execute_function_in_parallel(
                generate_reactions, rct_ichs, (allowed_prd_ichs, rtyp),
                nprocs=nprocs)

        # Add stereo
        if stereo:
            rxns = ()
            for _rxn in ini_rxns:
                log1 = ('\nExpanding Stereo for Reaction: '
                        f'{format_rxn_name(_rxn)}\n')
                _tmp_rxn = (_rxn[0], _rxn[1])
                thrdbdy = _rxn[2]
                ste_rxns_lst, log2 = _ste_rxn_lsts(_tmp_rxn)
                ste_rxns_lst, rem_rxns, log3 = _remove_enantiomer_reactions(
                    ste_rxns_lst, reacs_stereo_inchi=_rxn[0])
                rxns += _add_third(ste_rxns_lst, thrdbdy)
                log4 = _stereo_results(_rxn, ste_rxns_lst, rem_rxns)
                print(log1 + log2 + log3 + log4)
        else:
            rxns = ini_rxns

        # Update the mechanism objects with unique spc and rxns
        mech_spc_dct = update_spc_dct_from_reactions(
            rxns, mech_spc_dct)
        mech_rxn_dct = update_rxn_dct(
            rxns, mech_rxn_dct, mech_spc_dct)

        print('\n---------------------------------------------------------\n')

    return mech_spc_dct, mech_rxn_dct


def generate_reactions(allowed_prd_ichs, rtyp, rct_ich_lst, output_queue=None):
    """ For a given reactants
    """

    for rct_ichs in rct_ich_lst:
        _rct_smis = tuple(map(automol.chi.smiles, rct_ichs))
        log = f'Generating Reactions for {_rct_smis}...'

        # Generate the products with the desired reactant and reaction type
        rct_gras = _rct_gras(rct_ichs)
        prd_ichs = _prd_ichs(rct_gras, rtyp)

        # Build list of generated reactions including reactants and products
        rxn_ichs = ()
        for pidx, prds in enumerate(prd_ichs):
            # Don't add if requested products not found
            if allowed_prd_ichs:
                if not all(prd in prds for prd in allowed_prd_ichs):
                    continue
            # Don't add if self reaction was generated; weak check: ignores stereo
            if set(rct_ichs) == set(prds):
                continue
            # If continue not hit, save reaction to list and print to stdout
            _prd_smis = tuple(map(automol.chi.smiles, prds))
            log += f'\nFound Product(s) {pidx+1}: {_prd_smis}'

            rxn_ichs += ((rct_ichs, prds, (None,)),)

        if not prd_ichs:
            log += '\nNO Product(s) Found'
        print(log)
        output_queue.put(rxn_ichs)
        # return rxn_ichs


# Functions to set lists for mech building step
def _determine_reactants(spc_dct, rct1_set, rct2_set, rtyp):
    """ Determine a list of reactants for one or more reactions
        these could be unimolecular or bimolecular.
    """

    def _gen_set(spc_dct, spc_set):
        """ Make the list of species that will serve as
            reactant 1 or reactant 2
        """

        # Check if string id of class of species (e.g., 'all', 'radicals'), or
        # simply list of species name given
        if isinstance(spc_set, str):

            # Build ini list from input or spc dct
            spc_names, spc_ichs = (), ()
            for name, dct in spc_dct.items():
                spc_names += (name,)
                spc_ichs += (dct['inchi'],)

            # Trim the list if needed
            if spc_set == 'radicals':
                spc_ichs, spc_names = _radicals(spc_ichs, spc_names)

        else:
            spc_names = spc_set
            spc_ichs = tuple(spc_dct[name]['inchi'] for name in spc_set)

        return spc_ichs, spc_names

    rct1_ichs, rct1_names = _gen_set(spc_dct, rct1_set)
    if automol.ReactionClass.is_bimolecular(rtyp):
        rct2_ichs, rct2_names = _gen_set(spc_dct, rct2_set)
        rxn_ichs = tuple(itertools.product(rct1_ichs, rct2_ichs))
        rxn_names = tuple(itertools.product(rct1_names, rct2_names))
    else:
        rxn_ichs = tuple((ich,) for ich in rct1_ichs)
        rxn_names = tuple((name,) for name in rct1_names)

    return rxn_ichs, rxn_names


def _determine_allowed_products(spc_dct, allowed_prds):
    """ Takes a list of species mechanism names for the products
        that must be produced in a reaction for a generated reaction
        to be saved and generates the corresponding InChI strings.
    """
    if allowed_prds is not None:
        spc_ichs = tuple(spc_dct[name]['inchi'] for name in allowed_prds)
    else:
        spc_ichs = ()

    return spc_ichs


# Helper functions
def _radicals(ich_lst, name_lst):
    """ Determine the radicals
    """
    rad_ichs, rad_names = (), ()
    for ich, name in zip(ich_lst, name_lst):
        if automol.graph.is_radical_species(automol.chi.graph(ich)):
            rad_ichs += (ich,)
            rad_names += (name,)

    return rad_ichs, rad_names


def _rct_gras(rct_ichs):
    """ Get reactant graphs from smiles
    """

    rct_gras = list(map(automol.amchi.graph, rct_ichs))
    rct_gras = list(map(automol.graph.without_stereo, rct_gras))
    rct_gras = list(map(automol.graph.explicit, rct_gras))
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)

    return rct_gras


def _prd_ichs(rct_gras, rxn_class_typ, check=True):
    """ Check the products
    """

    # Enumerate all possible reactions checking they are correct
    rxns = automol.reac.enumerate_reactions(rct_gras, rxn_type=rxn_class_typ)

    # Get InChI strings from the enumerated reaction objects.
    # Check for validity if requested
    prd_ichs = ()
    for rxn in rxns:
        prd_gras_ = automol.reac.product_graphs(rxn, stereo=False)
        prd_ichs_ = tuple(map(automol.graph.chi, prd_gras_))
        prd_ichs += (prd_ichs_,)

        if check:
            rct_gras_ = automol.reac.reactant_graphs(rxn, stereo=False)
            rxns_ = automol.reac.find(rct_gras_, prd_gras_)
            try:
                assert rct_gras_ == rct_gras
                assert any(automol.reac.class_(r) == rxn_class_typ for r in rxns_)
            except AssertionError:
                print('WARNING: issue with reaction')

    # Remove duplicates (some reactions could have Reaction Objects, mult TS)
    prd_ichs = automol.util.remove_duplicates_with_order(prd_ichs)

    return prd_ichs
