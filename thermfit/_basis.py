"""
  Therm Calculations
"""

import os
import automol.inchi
import automol.geom
from phydat import phycon
import autorun
from mechanalyzer.inf import rxn as rinfo
import thermfit.cbh


# Set up the references for building
def prepare_basis(ref_scheme, spc_dct, spc_names,
                  zrxn=None, print_log=True, nprocs=1):
    """ Main callable function to generate the reference
        basis molecules used in heat-of-formation calculations

        Executed in parallel using Python MultiProcessing
        module if nprocs set above 1. Parallelization occurs
        over the species list.
    """

    args = (ref_scheme, spc_dct, zrxn, print_log)
    basis_dcts = autorun.execute_function_in_parallel(
        _prepare_basis, spc_names, args, nprocs=nprocs)

    basis_dct = {}
    for dct in basis_dcts:
        basis_dct.update(dct)

    return basis_dct


def _prepare_basis(ref_scheme, spc_dct, zrxn, print_log,
                   spc_names, output_queue=None):
    """ Prepare references
    """

    # Get objects from spc_names
    spc_str = ', '.join(spc_names)
    spc_ichs = [spc_dct[spc]['inchi'] for spc in spc_names]

    # Begin prints
    if print_log:
        print('Process {} preppng species: {}'.format(os.getpid(), spc_str))

    # Print the message
    msg = '\nDetermining reference molecules for scheme: {}'.format(ref_scheme)
    msg += '\n'

    basis_dct = {}
    for spc_name, spc_ich in zip(spc_names, spc_ichs):

        msg += '\nDetermining basis for species: {}'.format(spc_name)

        # Build the basis set and coefficients for spc/TS
        if zrxn is not None:
            rcls = automol.reac.reaction_class(zrxn)
            radrad = automol.reac.is_radical_radical(zrxn)
            if (rcls, radrad) in thermfit.cbh.CBH_TS_CLASSES:
                scheme = ref_scheme
            else:
                scheme = 'basic'
            if '_' in scheme:
                scheme = 'cbh' + scheme.split('_')[1]
            spc_basis, coeff_basis = thermfit.cbh.ts_basis(
                zrxn, scheme)
        else:
            spc_basis, coeff_basis = thermfit.cbh.species_basis(
                spc_ich, ref_scheme)
        spc_basis = tuple(automol.inchi.add_stereo(bas) for bas in spc_basis
                          if isinstance(bas, str))

        msg += '\nInChIs for basis set: {}'.format('\n  '.join(spc_basis))

        # Add to the dct containing info on the species basis
        basis_dct[spc_name] = (spc_basis, coeff_basis)

    # Print log message if it is desired
    if print_log:
        print(msg)

    # Add results to multiprocessing.Queue obj passed in
    output_queue.put((basis_dct,))


def unique_references(basis_dct, spc_dct):
    """ Take one or more basis dictionaries and obtain the uniqueness of them
    """

    def _spc_ref_unique(ref, dct_ichs, bas_ichs):
        """ Determine if InChI is unique
        """
        return all(
            ref not in lst
            for lst in (spc_ichs, dct_ichs, bas_ichs))

    def _ts_ref_unique(ref, dct_ichs, bas_ichs):
        """ Determine if lst of inchis for ts is unique """
        return (
            _ich_unique(ref, dct_ichs, bas_ichs) or
            _ich_unique(ref[::-1], dct_ichs, bas_ichs[::-1])
        )

    # Generate list of all species currently in the spc dct
    dct_ichs = [spc_dct[spc]['inchi'] for spc in spc_dct.keys()
                if 'ts' not in spc]

    # Generate list of all prospective basis species
    bas_ichs = ()
    cnt = 1
    for ref_ich, bas_ichs in basis_dct.items():
        if isinstance(bas, str):
            # Generate info for the species 
            if _spc_ref_unique(ref, dct_ichs, bas_ichs):
                ref_name = 'REF_{}'.format(cnt)
                unique_refs_dct[ref_name] = create_spec(ref)
                bas_ichs += ichs
                cnt +=1
        else:
            # Generate info for TS if unique
            if _ts_ref_unique(ref, dct_ichs, bas_ichs):
                ref_name = 'TS_REF_{}'.format(cnt)
                unique_refs_dct[ref_name] = create_ts_spc(
                    ref, spc_dct, spc_dct[spc_name]['mult'],
                    rcls)
                bas_ichs += ichs
                cnt +=1

    return unique_refs_dct


# Create dictionries of information for species/TS
def create_ts_spc(ref, spc_dct, mult, rxnclass):
    """ add a ts species to the species dictionary
    """

    # Obtain the Reaction InChIs, Charges, Mults
    reacs, prods = ref[0], ref[1]
    rxn_ichs = (
        tuple(automol.inchi.add_stereo(ich) for ich in reacs if ich),
        tuple(automol.inchi.add_stereo(ich) for ich in prods if ich)
    )

    rxn_muls, rxn_chgs = (), ()
    for rgts in (reacs, prods):
        rgt_muls, rgt_chgs = (), ()
        for rgt in rgts:
            found = False
            for dct in spc_dct.values():
                if 'inchi' in dct:
                    if dct['inchi'] == rgt:
                        rgt_muls += (dct['mult'],)
                        rgt_chgs += (dct['charge'],)
                        found = True
                        break
            if not found:
                new_spc = create_spec(rgt)
                rgt_muls += (new_spc['mult'],)
                rgt_chgs += (new_spc['charge'],)
        rxn_muls += (rgt_muls,)
        rxn_chgs += (rgt_chgs,)

    return {
        'reacs': list(reacs),
        'prods': list(prods),
        'charge': 0,
        'inchi': '',
        'class': rxnclass,
        'mult': mult,
        'ts_locs': (0,),
        'rxn_info': rinfo.from_data(rxn_ichs, rxn_chgs, rxn_muls, mult)
    }


def create_spec(ich, charge=0,
                mc_nsamp=(True, 3, 1, 3, 100, 12),
                hind_inc=30.):
    """ add a species to the species dictionary
    """
    rad = automol.formula.electron_count(automol.inchi.formula(ich)) % 2
    mult = 1 if not rad else 2

    return {
        'smiles': automol.inchi.smiles(ich),
        'inchi': ich,
        'inchikey': automol.inchi.inchi_key(ich),
        'charge': charge,
        'mult': mult,
        'mc_nsamp': mc_nsamp,
        'hind_inc': hind_inc * phycon.DEG2RAD
    }
