"""
Read the mechanism file
"""

import os
import errno
from functools import wraps
import pandas
import math
import multiprocessing
import random
import signal
import pandas as pd
import automol
from mechanalyzer.parser.csv_ import csv_dct
from mechanalyzer.parser.csv_ import read_csv_headers
from mechroutines.pf.thermo import basis


# HACK FOR A MECHANISM
# put in bad names for the hack
BAD_NAMES = []


# Build a spc dct containing all info species from an input file
def build_spc_dct(spc_str, spc_type):
    """ Get a dictionary of all the input species
        indexed by InChi string
    """

    # new_headers = _set_headers(spc_str)

    if spc_type == 'csv':
        spc_dct = csv_dct(spc_str)
    else:
        raise NotImplementedError

    return spc_dct


# Modify new file
def order_species_by_atomcount(spc_dct):
    """ Returns a species dictionary ordered by increasing N of atoms
        in the species. Useful for nicer mechanism writing
    """

    natom_df = pd.Series(index=list(spc_dct.keys()))
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich = spc_dct[key]['inchi']
            fml_dct = automol.inchi.formula(ich)
            natoms = automol.formula.atom_count(fml_dct)
            natom_df[key] = natoms
    natom_df = natom_df.sort_values(ascending=True)

    # Rewrite spc_dct
    sorted_idx = list(natom_df.index)
    sorted_val = list(map(spc_dct.get, sorted_idx))
    spc_dct_ordered = dict(zip(sorted_idx, sorted_val))

    return spc_dct_ordered


def write_basis_csv(spc_str, outname='species_hof_basis.csv',
                    path='.', parallel=False):
    """ Read the species file in a .csv format and write a new one
        that has hof basis species added to it.
    """

    headers = [header for header in read_csv_headers(spc_str)
               if header != 'name']
    if 'inchi' not in headers:
        headers.append('inchi')
    headers_noich = [header for header in headers
                     if header not in ('inchi', 'inchikey')]
    new_headers = ['inchi', 'inchikey'] + headers_noich
    # Read in the initial CSV information (deal with mult stereo)
    init_dct = csv_dct(spc_str, values=headers, key='name')
    # Build a new dict
    names = list(init_dct.keys())
    spc_queue = []
    for name in names:
        spc_queue.append([name, ''])
    new_dct = {}
    # Get the hof basis molecules
    ref_schemes = ['cbh0', 'cbh1', 'cbh2']
    ref_dct = {}
    spc_str = ','.join(['name'] + new_headers) + '\n'
    for ref_scheme in ref_schemes:
        formula_dct = {}
        _, uniref_dct = basis.prepare_refs(
            ref_scheme, init_dct, spc_queue, repeats=True, parallel=parallel)
        for name in uniref_dct:
            spc_str += ref_scheme + '_'
            smiles = automol.inchi.smiles(uniref_dct[name]['inchi'])
            uniref_dct[name]['smiles'] = smiles
            formula = automol.inchi.formula_string(uniref_dct[name]['inchi'])
            if formula in formula_dct:
                formula_dct[formula] += 1
                formula = formula + '({})'.format(formula_dct[formula])
            else:
                formula_dct[formula] = 0
            spc_str += formula + ','
            for idx, header in enumerate(new_headers):
                val = uniref_dct[name][header]
                if isinstance(val, str):
                    val = "'{}'".format(val)
                spc_str += str(val)
                if idx+1 < len(new_headers):
                    spc_str += ','
            spc_str += '\n'

    basis_file = os.path.join(path, outname + '_basis')
    with open(basis_file, 'w') as file_obj:
        file_obj.write(spc_str)

    new_names = []
    for ref_scheme in ref_schemes:
        _, uniref_dct = basis.prepare_refs(
            ref_scheme, init_dct, spc_queue, repeats=False, parallel=parallel)
        for newn in list(uniref_dct.keys()):
            tempn_smiles = automol.inchi.smiles(uniref_dct[newn]['inchi'])
            tempn = ref_scheme + '_' + tempn_smiles
            if tempn not in new_names:
                new_names.append(tempn)
                new_dct[tempn] = uniref_dct[newn]
                new_dct[tempn]['smiles'] = tempn_smiles
                ref_dct.update(new_dct)
                new_dct.update(init_dct)
                init_dct = new_dct
    # Writer string
    spc_str = ','.join(['name'] + new_headers) + '\n'
    for name in init_dct:
        formula_dct = {}
        if 'cbh' in name:
            namelabel = name.split('_')[0]
            if namelabel not in formula_dct:
                formula_dct[namelabel] = {}
            frmdct = formula_dct[namelabel]
            spc_str += namelabel + '_'
            formula = automol.inchi.formula_string(init_dct[name]['inchi'])
            if formula in frmdct:
                frmdct[formula] += 1
                formula = formula + '({})'.format(frmdct[formula])
            else:
                frmdct[formula] = 0
            spc_str += formula + ','
        else:
            spc_str += '{},'.format(name)
        for idx, header in enumerate(new_headers):
            val = init_dct[name][header]
            if isinstance(val, str):
                val = "'{}'".format(val)
            spc_str += str(val)
            if idx+1 < len(new_headers):
                spc_str += ','
        spc_str += '\n'

    # Write the file
    basis_file = os.path.join(path, outname)
    with open(basis_file, 'w') as file_obj:
        file_obj.write(spc_str)


def write_stereo_csv(spc_str, outname='species_stereo.csv', path='.',
                     allstereo=False):

    # Build a stereochemical dictionary
    new_dct, new_headers, names_in_order = write_stereo_dct(
        spc_str, allstereo=allstereo)

    # Writer string
    spc_str = ','.join(['name'] + new_headers) + '\n'
    for name in names_in_order:
        if name in new_dct:
            spc_str += '{},'.format(name)
            for idx, header in enumerate(new_headers):
                val = new_dct[name][header]
                if isinstance(val, str):
                    val = "'{}'".format(val)
                spc_str += str(val)
                if idx+1 < len(new_headers):
                    spc_str += ','
            spc_str += '\n'

    # Write the file
    stereo_file = os.path.join(path, outname)
    with open(stereo_file, 'w') as file_obj:
        file_obj.write(spc_str)


def write_stereo_dct(spc_str, allstereo=False):
    """ read the species file in a .csv format and write a new one
        that has stero information
    """

    # Read the headers
    headers = [header for header in read_csv_headers(spc_str)
               if header != 'name']
    if 'inchi' not in headers:
        headers.append('inchi')
    headers_noich = [header for header in headers
                     if header not in ('inchi', 'inchikey')]
    new_headers = ['inchi', 'inchikey'] + headers_noich

    # Read in the initial CSV information (deal with mult stereo)
    init_dct = csv_dct(spc_str, values=headers, key='name')

    # Build a new dict
    new_dct = {}
    names_in_order = list(init_dct.keys())
    randomized_names = list(init_dct.keys())
    nproc_avail = max(len(os.sched_getaffinity(0)) - 1, 1)
    num_spc = len(randomized_names)
    spc_per_proc = math.floor(num_spc / nproc_avail)

    queue = multiprocessing.Queue()
    procs = []
    random.shuffle(randomized_names)
    for proc_n in range(nproc_avail):
        spc_start = proc_n*spc_per_proc
        if proc_n == nproc_avail - 1:
            spc_end = num_spc
        else:
            spc_end = (proc_n+1)*spc_per_proc
        names = randomized_names[spc_start:spc_end]

        proc = multiprocessing.Process(
            target=_add_stereo_to_dct,
            args=(queue, names, init_dct, headers_noich, allstereo,))
        procs.append(proc)
        proc.start()

    for _ in procs:
        new_dct.update(queue.get())
    for proc in procs:
        proc.join()

    return new_dct, new_headers, names_in_order


def _add_stereo_to_dct(queue, names, init_dct, headers_noich, allstereo):

    print('Processor {} will check species: {}'.format(os.getpid(),
                                                       ', '.join(names)))
    good_names = []
    for name in names:
        if name not in BAD_NAMES:
            print('getting nrings for {}'.format(name))
            nrings = get_nrings(name, init_dct[name])
            if nrings < 2:
                good_names.append(name)
            else:
                print('{} has {:g} rings, '.format(name, nrings),
                      'stereo not implemented')
    new_dct = {}

    print('Processor {} will add stereo to species: {}'.format(
        os.getpid(), ', '.join(good_names)))
    for name in good_names:
        # Get the inchi key
        ich = init_dct[name]['inchi']

        # Generate ichs with stereo and hashes
        ichs_wstereo, succeeded = _generate_stereo(ich, allstereo=allstereo)
        if not succeeded:
            print('{} timed out in stereo generation'.format(name))
            continue

        # Add name and inchi info to string
        for idx, ich_wstereo in enumerate(ichs_wstereo):

            # Augment name if needed
            if idx == 0:
                sname = name
            else:
                sname = name + '-{}'.format(str(idx+1))

            # Initialize
            new_dct[sname] = {}

            # Generate hash key from InChI
            hashkey = automol.inchi.inchi_key(ich_wstereo)

            # Add vals to dct
            new_dct[sname].update({'inchi': ich_wstereo, 'inchikey': hashkey})

            for header in headers_noich:
                new_dct[sname][header] = init_dct[name][header]
    queue.put(new_dct)
    print('Processor {} finished'.format(os.getpid()))


# HELPER FUNCTIONS
def _set_headers(spc_str):
    """
    """
    headers = [header for header in read_csv_headers(spc_str)
               if header != 'name']
    if 'inchi' not in headers:
        headers.append('inchi')
    headers_noich = [header for header in headers
                     if header not in ('inchi', 'inchikey')]
    new_headers = ['inchi', 'inchikey'] + headers_noich

    return new_headers


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            print(error_message)
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator


@timeout(200)
def _generate_stereo(ich, allstereo=False):
    """ stereo
    """
    print('start ich {}'.format(ich))
    ret_ichs = []
    worked = True
    try:
        if not automol.inchi.is_complete(ich):
            if allstereo:
                ret_ichs = automol.inchi.expand_stereo(ich)
            else:
                ret_ichs = [automol.inchi.add_stereo(ich)]
        else:
            ret_ichs = [ich]
    except:
        print('Time out on ich {}'.format(ich))
        ret_ichs = ich
        worked = False
    return ret_ichs, worked


@timeout(30)
def get_nrings(name, dct):
    """ Count the number of rings
    """
    nrings = 1000
    try:
        nrings = len(automol.graph.rings(
            automol.inchi.graph(dct['inchi'])))
    except:
        print('Cannot produce graph for {} '.format(name))
        nrings = 2000
    return nrings
