""" Modifies the species.csv file in ways:
    Transform existing rows:
    (0) adds inchis if only SMILES are present
    (1) (optional)adds canonical enantiomer inchis to the species in the file csv
        (the canonical enantiomer is which enantiomer will be used for calculations to
        avoid reduandant calculations)
    (2) (optional) adds stereochemistry to species in the file csv (chooses one enantiomer/diastereomer)
        (see ste_mech script for expansion of all stereochemistry)
    (3) (optional) sorts the species in the csv file by stoichiometries
    (4) (optional) turns inchis into amchis that cannot be described by inchis (e.g., resonance)
    Add new rows:
    (5) (optional) adds required heat-of-formation basis species not present in csv file
    (6) (optional) adds instability product species to the file csv
"""

import os
import time
from ioformat import pathtools
from mechanalyzer import parser


def main(
        input_fname: str='species.csv',
        output_fname: str='mod_species.csv',
        sort: bool=False,
        include_stereo: bool=False,
        include_canonical: bool=False,
        use_amchi: bool=False,
        expand_hof_basis: bool=False,
        expand_instability: bool=False,
        ncpus: int=1,
        ):
    # Initialize the start time for script execution
    t0 = time.time()
    cwd = os.getcwd()

    # Read input species file into a species dictionary and add
    # necessary information like inchis if only smiles are present
    # and useful information like inchikey
    print(f'Reading species from {input_fname}...')
    spc_str = pathtools.read_file(cwd, input_fname)
    mech_spc_dct = parser.new_spc.parse_mech_spc_dct(spc_str)
    new_mech_headers = ('smiles', 'inchi', 'inchikey', 'mult', 'charge')

    # Add species that unstable will directly decompose into
    if expand_instability:
        print('Adding instability products to species dictionary...')
        mech_spc_dct = parser.spc.add_instability_products(
            mech_spc_dct, nprocs=ncpus, stereo=True)

    # Convert inchis to amchis for species that cannot be represented by inchis
    if use_amchi:
        print('Converting InChIs to AMChIs where necessary...')
        mech_spc_dct = parser.new_spc.mech_inchi_to_amchi(mech_spc_dct)

    # Add a stereochemical label to any stereochemical species
    if include_stereo:
        print('Adding stereochemical information to species dictionary...')
        mech_spc_dct = parser.spc.stereochemical_spc_dct(
            mech_spc_dct, nprocs=ncpus, all_stereo=False)

    # Sort the species dictionary, if requested
    if sort:
        print('Sorting species dictionary by atom count...')
        mech_spc_dct = parser.spc.reorder_by_atomcount(
            mech_spc_dct)

    if include_canonical:
        print('Adding canonical enantiomer InChIs to species dictionary...')
        mech_spc_dct = parser.new_spc.add_canonical_enantiomer(
            mech_spc_dct)
        new_mech_headers += ('canon_enant_ich',)

    # Add the thermochemical species to the species dictionary
    if expand_hof_basis:
        print('Adding heat-of-formation basis species to species dictionary...')
        if not include_canonical:
            mech_spc_dct = parser.new_spc.add_canonical_enantiomer(
                mech_spc_dct, dummy=True)
        # mech_spc_dct = parser.spc.add_heat_of_formation_basis(
        #    mech_spc_dct, ref_schemes=('cbh0', 'cbh1'),
        mech_spc_dct = parser.spc.add_heat_of_formation_basis(
            mech_spc_dct, ref_schemes=('cbh0', 'cbh1', 'cbh2'),
            nprocs=ncpus)

    # Write the new species dictionary to a string
    csv_str = parser.spc.csv_string(mech_spc_dct, new_mech_headers)

    # Write the string to a file
    pathtools.write_file(csv_str, cwd, output_fname)

    # Compute script run time and print to screen
    tf = time.time()
    print(f'\nSuccess: {input_fname} has been processed an updated in {output_fname}.')
    print(f'Time to complete: {tf-t0:.2f}')
