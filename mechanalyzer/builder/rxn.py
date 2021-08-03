""" Generate sets of reactions from reactants to construct mechanisms
"""

import itertools
import automol
from mechanalyzer.builder._update import update_spc_dct_from_reactions
from mechanalyzer.builder._update import update_rxn_dct


# MAIN CALLABLE FUNCTIONS FOR GENERATING REACTION LISTS
def build_mechanism(mech_spc_dct, mech_rxn_dct, rxn_series):
    """ Use the lst of reactions to build objects to describe mechanism
    """

    # Initialize new_spc
    new_spc_lst = None

    print('---------------------------------------------------------\n')

    # Loop over the reaction series consisting of reactants and reaction type
    for sidx, series in enumerate(rxn_series):

        rct1_set, rct2_set, allowed_prds, rxn_typs = series
        print('allowed prds', allowed_prds)

        print('Generating Reactions for Series {}'.format(sidx+1))

        # Generate the reactions from the reactants
        rxns = ()
        for rtyp in rxn_typs:

            # Determine reactants to generate reactions for
            rct_ichs = _determine_reactants(
                mech_spc_dct, rct1_set, rct2_set, rtyp, rct1_lst=new_spc_lst)
            allowed_prd_ichs = _determine_allowed_products(
                mech_spc_dct, allowed_prds)

            # Print info message for generating reactions usings SMILES strings
            print('\nTrying to find {} products for reactants'.format(rtyp))
            for ichs in rct_ichs:
                _rct_smis = tuple(map(automol.inchi.smiles, ichs))
                print(_rct_smis)
            print('')
            if allowed_prd_ichs:
                allow_str = ' '.join(allowed_prd_ichs)
                print('\nEach reaction must produce {}'.format(allow_str))
                print('')

            # Generate Reactions
            for ichs in rct_ichs:
                rxns += generate_reactions(ichs, allowed_prd_ichs, rtyp)

        # Update the mechanism objects with unique spc and rxns
        mech_spc_dct, new_spc_lst = update_spc_dct_from_reactions(
            rxns, mech_spc_dct)
        mech_rxn_dct = update_rxn_dct(
            rxns, mech_rxn_dct, mech_spc_dct)

        print('\n---------------------------------------------------------\n')

    return mech_spc_dct, mech_rxn_dct


def generate_reactions(rct_ichs, allowed_prd_ichs, rtyp):
    """ For a given reactants
    """

    _rct_smis = tuple(map(automol.inchi.smiles, rct_ichs))
    print('Generating Reactions for {}...'.format(_rct_smis))

    # Generate the products with the desired reactant and reaction type
    rct_gras = _rct_gras(rct_ichs)
    prd_ichs = _prd_ichs(rct_gras, rtyp)

    # Build list of generated reactions including reactants and products
    rxn_ichs = ()
    for pidx, prds in enumerate(prd_ichs):
        # Move ahead in loop if requested products not found
        if allowed_prd_ichs:
            if not all(prd in prds for prd in allowed_prd_ichs):
                continue
        # If continue not hit, save reaction to list and print to stdout
        _prd_smis = tuple(map(automol.inchi.smiles, prds))
        print('Found Product(s) {}: {}'.format(pidx+1, _prd_smis))

        rxn_ichs += ((rct_ichs, prds, (None,)),)

    return rxn_ichs


# Functions to set lists for mech building step
def _determine_reactants(spc_dct, rct1_set, rct2_set, rtyp, rct1_lst=None):
    """ Determine a list of reactants for one or more reactions
        these could be unimolecular or bimolecular.
    """

    def _gen_set(spc_dct, spc_set):
        """ Make the list of species that will serve as
            reactant 1 or reactant 2
        """

        if isinstance(spc_set, str):

            # Build ini list from input or spc dct
            if rct1_lst is not None:
                spc_ichs = rct1_lst
            else:
                spc_ichs = tuple(dct['inchi'] for dct in spc_dct.values())

            # Trim the list if needed
            if spc_set == 'radicals':
                spc_ichs = _radicals(spc_ichs)

        else:
            spc_ichs = tuple(spc_dct[name]['inchi'] for name in spc_set)

        return spc_ichs

    rct1_ichs = _gen_set(spc_dct, rct1_set)
    if automol.par.isbimol(rtyp):
        rct2_ichs = _gen_set(spc_dct, rct2_set)
        rxn_ichs = tuple(itertools.product(rct1_ichs, rct2_ichs))
    else:
        rxn_ichs = tuple((ich,) for ich in rct1_ichs)

    return rxn_ichs


def _determine_allowed_products(spc_dct, allowed_prds):
    """ get the inchis for what products are allowed
    """
    if allowed_prds is not None:
        spc_ichs = tuple(spc_dct[name]['inchi'] for name in allowed_prds)
    else:
        spc_ichs = ()

    return spc_ichs


# Helper functions
def _radicals(ich_lst):
    """ Determine the radicals
    """
    return tuple(
        ich for ich in ich_lst
        if automol.graph.radical_species(automol.inchi.graph(ich)))


def _rct_gras(rct_ichs):
    """ Get reactant graphs from smiles
    """

    rct_geos = list(map(automol.inchi.geometry, rct_ichs))
    rct_gras = tuple(map(automol.geom.connectivity_graph, rct_geos))
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
        prd_gras_ = automol.reac.product_graphs(rxn)
        prd_ichs_ = tuple(map(automol.graph.inchi, prd_gras_))
        prd_ichs_ = tuple(map(automol.inchi.add_stereo, prd_ichs_))
        prd_ichs += (prd_ichs_,)

        if check:
            rct_gras_ = automol.reac.reactant_graphs(rxn)
            rxns_ = automol.reac.find(rct_gras_, prd_gras_)
            try:
                assert rct_gras_ == rct_gras
                assert any(r.class_ == rxn_class_typ for r in rxns_)
            except AssertionError:
                print('WARNING: issue with reaction')

    # Remove duplicates (some reactions could have Reaction Objects, mult TS)
    prd_ichs = automol.util.remove_duplicates_with_order(prd_ichs)

    return prd_ichs
