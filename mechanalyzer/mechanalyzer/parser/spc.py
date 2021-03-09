"""
Read the mechanism file
"""

import os
import pandas as pd
import errno
from functools import wraps
import pandas
import math
import multiprocessing
import random
import signal
import automol
from mechanalyzer.parser.csv_ import csv_dct
from mechanalyzer.parser.csv_ import read_csv_headers
from mechroutines.pf.thermo import basis


# HACK FOR A MECHANISM
BAD_NAMES = ['C7H7(1132)', 'C5H6(1767)', 'S(2226)', 'C10H8(661)', 'C10H9(638)', 'C7H6(1119)', 'S(1060)', 'S(2355)', 'S(2477)', 'S(2468)', 'S(530)', 'C6H6(888)', 'C8H11(511)', 'C7H6(1542)', 'S(276)', 'N4-1(2604)', 'S(280)', 'S(1332)', 'S(2441)', 'C7H6(1512)', 'S(2454)', 'N5-2(2616)', 'S(271)', 'S(275)', 'W7_2(2646)', 'S(239)', 'S(1262)', 'S(2739)', 'S(2527)', 'S(1817)', 'S(2494)', 'W6_1(2639)', 'W4(2767)', 'ANT(2617)', 'S(2744)', 'S(2619)', 'S(1647)', 'S(1576)', 'W3_2(2647)', 'C9H9(2264)', 'S(1570)', 'W4_2(2651)', 'S(1892)', 'C8H8(1447)', 'S(1715)', 'C8H9(2357)', 'S(2540)', 'C7H6(1355)', 'S(1565)', 'S(2244)', 'S(2621)', 'C8H12(364)', 'P1(2755)', 'C9H9(2296)', 'W4(2759)', 'S(1273)', 'W1_2(2649)', 'C7H6(1867)', 'S(1824)', 'C8H9(517)', 'C9H13(931)', 'W1(2754)', 'S(2413)', 'N3-2(2614)', 'C5H5(1441)', 'C7H6(1835)', 'W5(2768)', 'S(2452)', 'S(1719)', 'C7H7(987)', 'S(2620)', 'S(2623)', 'S(270)', 'C7H7(1088)', 'C7H6(1560)', 'N3-1(2603)', 'P4_2(2648)', 'PTR(2605)', 'C8H12(352)', 'C8H7(2197)', 'C7H9(1172)', 'C9H9(2271)', 'C8H8(1156)', 'C7H6(1540)', 'S(2613)', 'S(2743)', 'S(277)', 'S(2514)', 'C10H9(617)', 'S(2399)', 'S(1626)', 'S(284)', 'S(2425)', 'S(1702)', 'S(2393)', 'S(2002)', 'S(2740)', 'S(2071)', 'S(1629)', 'S(286)', 'W1(2764)', 'W2(2756)', 'S(2742)', 'S(1630)', 'S(1428)', 'S(1727)', 'C9H9(2260)', 'S(2411)', 'C8H8(1423)', 'CHN(2607)', 'C8H11(513)', 'S(1891)', 'S(2392)', 'S(2237)', 'S(1666)', 'C10H9(615)', 'S(1648)', 'S(2236)', 'C7H7(1362)', 'W3(2758)', 'N9-1(2608)', 'S(1431)', 'P3_2(2650)', 'S(194)', 'C9H9(2295)', 'S(2630)', 'C7H6(1529)', 'S(2330)', 'S(1653)', 'S(1819)', 'S(2084)', 'C9H9(1846)', 'C7H6(990)', 'S(1264)', 'S(2105)', 'W7_1(2642)', 'C8H7(2195)', 'C4H5(833)', 'S(279)', 'S(2610)', 'S(2745)', 'S(2609)', 'C7H6(1545)', 'S(1659)', 'C8H12(355)', 'S(267)',
             'S(196)', 'S(2474)', 'S(2385)', 'S(2470)', 'S(1627)', 'P4_1(2641)', 'W4_1(2634)', 'W4_n(2635)', 'S(1276)', 'S(539)', 'S(2442)', 'C10H9(651)', 'N2-1(2558)', 'S(2631)', 'C9H9(2251)', 'C7H6(1527)', 'S(1070)', 'S(1652)', 'S(2629)', 'C7H8(1051)', 'N5-1(2606)', 'S(2423)', 'S(2485)', 'C16H18(54)', 'S(2644)', 'C18H22(57)', 'C10H9(639)', 'S(2384)', 'S(2696)', 'S(1759)', 'W3_1(2640)', 'C7H6(1324)', 'P3_1(2633)', 'S(2746)', 'S(281)', 'N4-2(2615)', 'S(1643)', 'S(1915)', 'C7H7(1099)', 'S(2343)', 'S(283)', 'C6H8(862)', 'C8H9(1155)', 'S(2391)', 'S(1664)', 'C9H13(936)', 'W7_n(2645)', 'S(1885)', 'S(282)', 'S(2612)', 'C8H9(725)', 'S(1828)', 'S(2611)', 'N9-2(2618)', 'C8H9(724)', 'W5(2760)', 'S(2738)', 'S(533)', 'S(2223)', 'C7H6(1494)', 'S(1748)', 'C18H22(58)', 'S(2536)', 'S(2490)', 'W6_n(2638)', 'P2(2757)', 'P1(2765)', 'S(2165)', 'C7H7(1089)', 'W3(2766)', 'S(2302)', 'S(1736)', 'S(1638)', 'S(1771)', 'C9H11(994)', 'S(2643)', 'W1_1(2632)', 'N7-1(2561)', 'S(2339)', 'C9H9(2276)', 'S(2622)', 'S(2695)', 'S(218)', 'S(2741)', 'C8H7(2415)', 'N2-2(2570)', 'S(278)', 'C7H6(1538)', 'S(2146)', 'C7H7(1111)', 'S(195)', 'W6_2(2653)', 'C8H8(1157)', 'S(1661)', 'pdt7(2700)', 'C7H9(1178)', 'S(1731)', 'S(1687)', 'S(1724)', 'C9H11(960)', 'C7H6(1353)', 'C8H9(1837)', 'S(2145)', 'S(1651)', 'S(1310)', 'S(1685)', 'S(1689)', 'S(2516)', 'S(2453)', 'S(1788)', 'C7H9(1209)', 'S(1931)', 'S(1996)', 'S(266)', 'S(2382)', 'S(1954)', 'C7H7(1354)', 'S(1361)', 'S(1728)', 'S(1657)', 'S(1893)', 'S(268)', 'S(2172)', 'S(1682)', 'S(2173)', 'S(2729)', 'S(1713)', 'S(240)', 'S(1716)', 'S(2567)', 'S(1402)', 'S(2697)', 'S(2255)', 'S(1683)', 'C10H9(755)', 'C7H8(937)', 'C7H7(1161)', 'S(238)', 'C8H9(2358)', 'C7H6(1492)', 'S(1403)', 'S(2447)', 'S(1722)', 'S(138)', 'C18H22(59)', 'S(1425)', 'S(1742)', 'S(1342)', 'S(532)']


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
def modify_csv(spc_dct, mod_lst=()):
    """ Take a CSV and add species and info as desired
    """

    # Modify the entries of the dct

    # Add species


def order_species_by_atomcount(spc_dct):
    """
    Returns a species dictionary ordered by increasing N of atoms in the species
    Useful for nicer mechanism writing
    """
    spc_atom_N = pd.Series(index=list(spc_dct.keys()))
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich = spc_dct[key]['inchi']
            fml_dct = automol.inchi.formula(ich)
            N_atoms = automol.formula.atom_count(fml_dct)
            spc_atom_N[key] = N_atoms

    spc_atom_N = spc_atom_N.sort_values(ascending=True)
    # rewrite spc_dct
    sorted_idx = list(spc_atom_N.index)
    sorted_val = list(map(spc_dct.get, sorted_idx))
    spc_dct_ordered = dict(zip(sorted_idx, sorted_val))

    return spc_dct_ordered


def write_basis_csv(spc_str, outname='species_hof_basis.csv', path='.', parallel=False):
    """read the species file in a .csv format and write a new one
    that has hof basis species added to it
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
            if not namelabel in formula_dct:
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
                print('{} has {:g} rings, stereo not implemented'.format(name, nrings))
    new_dct = {}

    print('Processor {} will add stereo to species: {}'.format(os.getpid(),
                                                               ', '.join(good_names)))
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
                ret_ichs = automol.inchi.add_stereo(ich)
            else:
                ret_ichs = [automol.inchi.add_stereo(ich)[0]]
        else:
            ret_ichs = [ich]
    except:
        print('Time out on ich {}'.format(ich))
        ret_ichs = ich
        worked = False
    return ret_ichs, worked


@timeout(30)
def get_nrings(name, dct):
    nrings = 1000
    try:
        nrings = len(automol.graph.rings(
            automol.inchi.graph(dct['inchi'])))
    except:
        print('Cannot produce graph for {} '.format(name))
        nrings = 2000
    return nrings
