""" Handles the construction, parsing, and manipulation of
    mechanism species dictionaries structured as

    {name: spc_data_dct}

    Also handles ancillary dictionaries used contain and
    organize information related to the species dictionary
"""

import os
import copy
import csv
import pandas as pd
import automol
import thermfit
from mechanalyzer.parser.csv_ import csv_dct
from autorun import timeout, execute_function_in_parallel
import ioformat.pathtools as text_parser


# Write a spc_dct to a CSV string
def csv_string(spc_dct, headers):
    """ Convert a mechanism species dictionary to a CSV file formatted
        string. The headers provided correspond to the physical data
        or parameters that will be parsed from the dictionary and written
        into the file.

        Note that the first column of the CSV string will always
        be 'name' and so that should be precluded from the list. The
        remaining headers will be written in the order provided.

        :param spc_dct: mechanism species dictionary
        :type spc_dct: dict[str: dict]
        :param headers: headers of the CSV file
        :type headers: tuple(str)
        :rtype: str
    """

    # Build spc sub dct of just info in the headers
    _csv_dct = {}
    for name, dct in spc_dct.items():
        _csv_dct[name] = {}
        for header in headers:
            _csv_dct[name][header] = dct.get(header, None)

    # Build the datagrame and resultant CSV string
    dframe = pd.DataFrame.from_dict(_csv_dct, orient='index')
    _csv_str = dframe.to_csv(
        index_label='name',
        quotechar="'",
        quoting=csv.QUOTE_NONNUMERIC)

    return _csv_str


# Build spc_dct from file/string i/o
def load_spc_dcts(spc_csv_filenames, direc):
    """ Reads one or more spc.csv files (in the AutoMech format). Outputs a
        list of spc_dcts.

        :param spc_csv_filenames: filenames containing spc info
        :type spc_csv_filenames: list [filename1, filename2, ...]
        :param direc: directory with file(s) (all must be in same directory)
        :type direc: str
        :return spc_dcts: list of spc_dcts
        :rtype: list of dcts [spc_dct1, spc_dct2, ...]
    """
    spc_dcts = []
    for spc_csv_filename in spc_csv_filenames:
        print('spc csv test:', spc_csv_filename)
        spc_csv_str = text_parser.read_file(direc, spc_csv_filename)
        spc_dct = build_spc_dct(spc_csv_str, 'csv')
        spc_dcts.append(spc_dct)

    return spc_dcts


def build_spc_dct(spc_str, spc_type):
    """ Get a dictionary of all the input species
        indexed by InChi string
    """

    if spc_type == 'csv':
        spc_dct = csv_dct(spc_str)
    else:
        raise NotImplementedError

    return spc_dct


# Build a spc dct from constituent information
def spc_dct_from_smiles(smiles_lst, stereo=False):
    """ Build a spc dct from a set of smiles
    """

    # Initialize empty formula dct
    fml_count_dct = {}

    spc_dct = {}
    for smi in smiles_lst:
        # Generate InChI string and formula
        ich = automol.smiles.inchi(smi)
        if stereo:
            ich = automol.inchi.add_stereo(ich)

        # Generate Name
        fml = automol.inchi.formula_string(ich)
        name, fml_count_dct = assign_unique_name(fml, fml_count_dct)

        # Add species dictionary
        spc_dct.update({name: thermfit.create_spec(ich)})

    return spc_dct


# Modify an existing spc_dct
def reorder_by_atomcount(spc_dct):
    """ Returns a species dictionary ordered by increasing N of atoms
        in the species. Useful for nicer mechanism writing
    """

    natom_df = pd.Series(index=list(spc_dct.keys()))
    for key in spc_dct.keys():
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


def add_heat_of_formation_basis(spc_dct,
                                ref_schemes=('cbh0', 'cbh1', 'cbh2'),
                                parallel=False):
    """ Adds all species required to form the basis set for requested CBHn
        heat-of-formation calculations for all species in the mechanism
        species dictionary.

        Function will loop over all species in the input dictionary,
        determine the basis and check if each basis currently in dictionary.

        If species is missing it is assigned a name:
        cbh{N}_{smiles} where smiles is the smiles string, and N is the
        lowest order cbh scheme for which a basis species was found.
    """

    def _filter_redundant_basis(cbh_ref_dct, ref_scheme, cbh_smiles):
        """ Get only spc of cbh_ref_dct that have smiles not already added
            to spc_dctfrom earlier CBHn schemes
        """
        new_dct = {}
        for dct in cbh_ref_dct.values():
            tempn_smiles = automol.inchi.smiles(dct['inchi'])
            if tempn_smiles not in cbh_smiles:
                # Add to new dictionry
                tempn = ref_scheme + '_' + tempn_smiles
                new_dct[tempn] = dct
                new_dct[tempn]['smiles'] = tempn_smiles
                # Add smiles to new smiles bookkeeping list
                cbh_smiles.append(tempn_smiles)

        return new_dct, cbh_smiles

    # Find all references
    cbh_smiles = []
    for ref_scheme in ref_schemes:
        _, cbh_ref_dct = thermfit.prepare_refs(
            ref_scheme, spc_dct, list(spc_dct.keys()),
            repeats=True, parallel=parallel)
        nonred_dct, cbh_smiles = _filter_redundant_basis(
            cbh_ref_dct, ref_scheme, cbh_smiles)
        spc_dct.update(nonred_dct)

    return spc_dct


def stereochemical_spc_dct(spc_dct, nprocs='auto', all_stereo=False):
    """ read the species file in a .csv format and write a new one
        that has stero information
    """

    # Generate names and procs for parallel stereochem add'n
    init_names = list(spc_dct.keys())

    # Add stereo using multiple processes
    args = (spc_dct, all_stereo)
    ste_dcts = execute_function_in_parallel(
        _add_stereo_to_dct, init_names, args, nprocs=nprocs)

    # Build and reoroder the stereo dct (only works if all_stereo=False)
    ste_spc_dct = {}
    for dct in ste_dcts:
        ste_spc_dct.update(dct)
    if not all_stereo:
        ste_spc_dct_ord = {name: ste_spc_dct[name] for name in init_names}
    else:
        ste_spc_dct_ord = copy.deepcopy(ste_spc_dct)

    return ste_spc_dct_ord


def _add_stereo_to_dct(init_dct, all_stereo, names, output_queue):
    """ Builds a modified species dictionary for a set of names where
        each sub species dictionary contains an InChI string with
        stereochemical layers being added.

        The function can either add all of the stereochemical species or
        just pick one.
    """

    # Functions for calling automol for ring counts and stereo add'n
    @timeout(30)
    def _nrings(name, dct):
        """ Count the number of rings
        """
        try:
            nrings = len(automol.graph.rings(
                automol.inchi.graph(dct['inchi'])))
        except:
            print('Cannot produce graph for {} '.format(name))
            nrings = 2000
        return nrings

    @timeout(200)
    def _generate_stereo(name, ich, all_stereo=False):
        """ stereo
        """
        ret_ichs, worked = [ich], True
        try:
            if not automol.inchi.is_complete(ich):
                ret_ichs = (
                    [automol.inchi.add_stereo(ich)] if not all_stereo else
                    automol.inchi.expand_stereo(ich))
        except:
            print('{} timed out in stereo generation'.format(name))
            worked = False
        return ret_ichs, worked

    # Assess the species the code is able to add stereochemistry to
    print('Processor {} will check species: {}'.format(os.getpid(),
                                                       ', '.join(names)))
    good_names = []
    for name in names:
        if _nrings(name, init_dct[name]) < 2:
            good_names.append(name)
        else:
            print('{} has >2 rings; stereo not implemented'.format(name))

    # Construct a new dictionary with stereochemical inchi strings
    new_dct = {}
    print('Processor {} will add stereo to species: {}'.format(
        os.getpid(), ', '.join(good_names)))
    for name in good_names:

        # Generate ichs with stereo and hashes
        ich = init_dct[name]['inchi']
        ste_ichs, success = _generate_stereo(name, ich, all_stereo=all_stereo)
        if not success:
            continue

        # Add name and inchi info to string
        spc_dct_keys = tuple(key for key in init_dct[name].keys()
                             if key != 'inchi')
        for idx, ich_wstereo in enumerate(ste_ichs):
            # name gen wrong, should use formula builder
            sname = name+'({})'.format(str(idx+1)) if idx != 0 else name
            new_dct[sname] = {'inchi': ich_wstereo}
            for key in spc_dct_keys:
                new_dct[sname][key] = init_dct[name][key]

    output_queue.put(new_dct)
    print('Processor {} finished'.format(os.getpid()))


# HELPER FUNCTIONS
# Handle Formula Count Objects
def formula_count_dct(spc_dct):
    """ get the spc_fml dct

        Loop over the spc dct, obtain the formula
        and use it to build a count of each formula
        in the dictionary
    """

    fml_count_dct = {}
    for dct in spc_dct.values():
        fml_str = automol.inchi.formula_string(dct['inchi'])
        _, fml_count_dct = assign_unique_name(
            fml_str, fml_count_dct)

    return fml_count_dct


def assign_unique_name(fml, fml_count_dct):
    """ Generate a unique name for a species with the
        given formula that corresponds to

        formula(N) where formula is the string N
        is the Nth iteration of said formula in the
        dictionary.

        Also, updates the overall formula dictionary
        which contains a count of how many times the
        formula appears in some mechanism.
    """

    if fml in fml_count_dct:
        fml_count_dct[fml] += 1
        fml = fml + '({})'.format(fml_count_dct[fml])
    else:
        fml_count_dct[fml] = 1

    return fml, fml_count_dct


# add functions like adding the hashkey
def add_hashkey(spc_dct):
    """ Add the inchi hashkey
    """

    for dct in spc_dct.values():
        ich = dct.get('inchi')
        ick = automol.inchi.inchi_key(ich) if ich is not None else None
        dct.update({'inchikey': ick})

    return spc_dct
