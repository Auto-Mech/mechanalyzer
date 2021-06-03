"""
  Therm Calculations
"""

import os
import math
import multiprocessing
import random
import automol.inchi
import automol.geom
from phydat import phycon
from mechanalyzer.inf import rxn as rinfo
from mechroutines.pf.thermo import heatform


# FUNCTIONS TO PREPARE THE LIST OF REFERENCE SPECIES NEEDED FOR THERM CALCS #
REF_CALLS = {
    "basic": "get_basic",
    "cbh0": "get_cbhzed",
    "cbh1": "get_cbhone",
    "cbh1_0": "get_cbhone",
    "cbh1_1": "get_cbhone",
    "cbh2": "get_cbhtwo",
    "cbh2_0": "get_cbhtwo",
    "cbh2_1": "get_cbhtwo",
    "cbh2_2": "get_cbhtwo",
    "cbh3": "get_cbhthree"
}

TS_REF_CALLS = {
    "basic": "get_basic_ts",
    "cbh0": "get_cbhzed_ts",
    "cbh1": "get_cbhone_ts",
    "cbh1_0": "get_cbhzed_ts",
    "cbh1_1": "get_cbhone_ts",
    "cbh2": "get_cbhzed_ts",
    "cbh2_0": "get_cbhzed_ts",
    "cbh2_1": "get_cbhone_ts",
    "cbh3": "get_cbhone_ts"
}

IMPLEMENTED_CBH_TS_CLASSES = [
    'hydrogen abstraction high',
    # 'hydrogen migration',
    'beta scission',
    'elimination high',
    'radical radical hydrogen abstraction high',
    'addition high'
]


def prepare_refs(
        ref_scheme, spc_dct, spc_queue, run_prefix, save_prefix,
        repeats=False, parallel=False, zrxn=None):
    """ add refs to species list as necessary
    """
    spc_names = [spc[0] for spc in spc_queue]

    if parallel:
        nproc_avail = len(os.sched_getaffinity(0)) - 1

        num_spc = len(spc_names)
        spc_per_proc = math.floor(num_spc / nproc_avail)

        queue = multiprocessing.Queue()
        procs = []
        random.shuffle(spc_names)
        for proc_n in range(nproc_avail):
            spc_start = proc_n*spc_per_proc
            if proc_n == nproc_avail - 1:
                spc_end = num_spc
            else:
                spc_end = (proc_n+1)*spc_per_proc

            spc_lst = spc_names[spc_start:spc_end]

            proc = multiprocessing.Process(
                target=_prepare_refs,
                args=(queue, ref_scheme, spc_dct, spc_lst,
                      run_prefix, save_prefix,
                      repeats, parallel, zrxn))
            procs.append(proc)
            proc.start()

        basis_dct = {}
        unique_refs_dct = {}
        for _ in procs:
            bas_dct, unq_dct = queue.get()
            basis_dct.update(bas_dct)
            bas_ichs = [
                unique_refs_dct[spc]['inchi']
                if 'inchi' in unique_refs_dct[spc]
                else unique_refs_dct['reacs']
                for spc in unique_refs_dct]
            for spc in unq_dct:
                new_ich = (
                    unq_dct[spc]['inchi']
                    if 'inchi' in unq_dct[spc] else unq_dct[spc]['reacs'])
                if new_ich not in bas_ichs:
                    cnt = len(list(unique_refs_dct.keys())) + 1
                    if isinstance(new_ich, str):
                        ref_name = 'REF_{}'.format(cnt)
                        unique_refs_dct[ref_name] = unq_dct[spc]
                    else:
                        ref_name = 'TS_REF_{}'.format(cnt)
                        unique_refs_dct[ref_name] = unq_dct[spc]
        for proc in procs:
            proc.join()
    else:
        basis_dct, unique_refs_dct = _prepare_refs(
            None, ref_scheme, spc_dct, spc_names,
            run_prefix, save_prefix,
            repeats=repeats, parallel=parallel,
            zrxn=zrxn)

    return basis_dct, unique_refs_dct


def _prepare_refs(queue, ref_scheme, spc_dct, spc_names,
                  run_prefix, save_prefix,
                  repeats=False, parallel=False, zrxn=None):
    """ Prepare references
    """

    print(
        'Processor {} will prepare species: {}'.format(
            os.getpid(), ', '.join(spc_names)))
    spc_ichs = [spc_dct[spc]['inchi'] for spc in spc_names]
    dct_ichs = [spc_dct[spc]['inchi'] for spc in spc_dct.keys()
                if spc != 'global' and 'ts' not in spc]

    # Determine the function to be used to get the thermochemistry ref species
    if ref_scheme in REF_CALLS:
        get_ref_fxn = getattr(heatform, REF_CALLS[ref_scheme])
    if ref_scheme in TS_REF_CALLS:
        get_ts_ref_fxn = getattr(heatform, TS_REF_CALLS[ref_scheme])

    # Print the message
    msg = '\nDetermining reference molecules for scheme: {}'.format(ref_scheme)
    msg += '\n'

    basis_dct = {}
    unique_refs_dct = {}
    for spc_name, spc_ich in zip(spc_names, spc_ichs):
        msg += '\nDetermining basis for species: {}'.format(spc_name)
        if zrxn is not None:
            rxnclass = automol.reac.reaction_class(zrxn)
            if (rxnclass in IMPLEMENTED_CBH_TS_CLASSES and
               'basic' not in ref_scheme):
                spc_basis, coeff_basis = get_ts_ref_fxn(zrxn)
            else:
                # Use a basic scheme
                spc_basis = []
                coeff_basis = []
                ts_ref_scheme = ref_scheme
                if '_' in ts_ref_scheme:
                    ts_ref_scheme = 'cbh' + ref_scheme.split('_')[1]
                for spc_i in spc_dct[spc_name]['reacs']:
                    bas_dct_i, _ = prepare_refs(
                        ts_ref_scheme, spc_dct, [[spc_i, None]],
                        run_prefix, save_prefix)
                    spc_bas_i, coeff_bas_i = bas_dct_i[spc_i]
                    for bas_i, c_bas_i in zip(spc_bas_i, coeff_bas_i):
                        if bas_i not in spc_basis:
                            spc_basis.append(bas_i)
                            coeff_basis.append(c_bas_i)
                        else:
                            for j, bas_j in enumerate(spc_basis):
                                if bas_i == bas_j:
                                    coeff_basis[j] += c_bas_i
        else:
            spc_basis, coeff_basis = get_ref_fxn(spc_ich)
        for i,  _ in enumerate(spc_basis):
            if isinstance(spc_basis[i], str):
                spc_basis[i] = automol.inchi.add_stereo(spc_basis[i])

        msg += '\nInCHIs for basis set:'
        for base in spc_basis:
            msg += '\n  {}'.format(base)

        # Add to the dct containing info on the species basis
        basis_dct[spc_name] = (spc_basis, coeff_basis)

        # Add to the dct with reference dct if it is not in the spc dct
        for ref in spc_basis:
            print('ref test', ref)
            bas_ichs = [
                unique_refs_dct[spc]['inchi']
                if 'inchi' in unique_refs_dct[spc]
                else unique_refs_dct[spc]['reacs']
                for spc in unique_refs_dct]
            cnt = len(list(unique_refs_dct.keys())) + 1
            if isinstance(ref, str):
                if ((ref not in spc_ichs and ref not in dct_ichs)
                        or repeats) and ref not in bas_ichs:
                    ref_name = 'REF_{}'.format(cnt)
                    msg += (
                        '\nAdding reference species {}, InChI string:{}'
                    ).format(ref, ref_name)
                    unique_refs_dct[ref_name] = create_spec(ref)
            else:
                if _chk(ref, spc_ichs, dct_ichs, bas_ichs, repeats):
                    ref_name = 'TS_REF_{}'.format(cnt)
                    msg += (
                        '\nAdding reference species {}, InChI string:{}'
                    ).format(ref, ref_name)
                    unique_refs_dct[ref_name] = create_ts_spc(
                        ref, spc_dct, spc_dct[spc_name]['mult'],
                        run_prefix, save_prefix,
                        rxnclass)
    print(msg)

    ret = None
    if parallel:
        queue.put((basis_dct, unique_refs_dct))
    else:
        ret = (basis_dct, unique_refs_dct)

    return ret


def create_ts_spc(ref, spc_dct, mult, run_prefix, save_prefix, rxnclass):
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
            for name in spc_dct:
                if 'inchi' in spc_dct[name]:
                    if spc_dct[name]['inchi'] == rgt:
                        rgt_muls += (spc_dct[name]['mult'],)
                        rgt_chgs += (spc_dct[name]['charge'],)
                        found = True
                        break
            if not found:
                new_spc = create_spec(rgt)
                rgt_muls += (new_spc['mult'],)
                rgt_chgs += (new_spc['charge'],)
        rxn_muls += (rgt_muls,)
        rxn_chgs += (rgt_chgs,)

    rxn_info = rinfo.from_data(rxn_ichs, rxn_chgs, rxn_muls, mult)

    return {
        'reacs': list(reacs),
        'prods': list(prods),
        'charge': 0,
        'inchi': '',
        'class': rxnclass,
        'mult': mult,
        'rxn_info': rxn_info,
        'ts_locs': (0,),
        'rxn_fs': reaction_fs(run_prefix, save_prefix, rinfo.sort(rxn_info))
    }


def create_spec(ich, charge=0,
                mc_nsamp=(True, 3, 1, 3, 100, 12),
                hind_inc=30.):
    """ add a species to the species dictionary
    """
    rad = automol.formula.electron_count(automol.inchi.formula(ich)) % 2
    mult = 1 if not rad else 2

    return {
        'inchi': ich,
        'inchikey': automol.inchi.inchi_key(ich),
        'sens': 0.0,
        'charge': charge,
        'mult': mult,
        'mc_nsamp': mc_nsamp,
        'hind_inc': hind_inc * phycon.DEG2RAD
    }


def is_scheme(scheme):
    """ Return Boolean val if there is a scheme
    """
    return bool(scheme in REF_CALLS)


# Helpers
def _chk(ref, spc_ichs, dct_ichs, bas_ichs, repeats):
    """ a """

    ini = (
        ((ref not in spc_ichs and ref not in dct_ichs) or repeats) and
        (ref not in bas_ichs)
    )
    rref = ref[::-1]
    sec = (
        ((rref not in spc_ichs and rref not in dct_ichs) or repeats) and
        (ref not in bas_ichs[::-1])
    )

    return ini or sec
