""" Script to generate files with reactions
"""

import os
import time

import ioformat
import chemkin_io
import mechanalyzer
from mechanalyzer.builder import sorter


def main(
        out_loc, out_spc, out_mech,
        inp_spc_str, inp_mech_str, sort_str,
        check_mechanism=False):
    """ carry out all the mechanism things you wanna do
    """
    # Build the initial dictionaries
    mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(inp_spc_str, 'csv')
    rxn_param_dct = mechanalyzer.parser.mech.parse_mechanism(
        inp_mech_str, 'chemkin')
    isolate_spc, sort_lst = mechanalyzer.parser.mech.parse_sort(sort_str)

    # Remove reactions that should not be there
    if check_mechanism:
        print('\n Removing improper reactions')
        rxn_param_dct = mechanalyzer.builder.remove_improper_reactions(
            rxn_param_dct, mech_spc_dct)

    # ADD STEREO TO SPC AND REACTIONS
    print('\n---- Adding stereochemistry to InChIs of mechanism'
          ' species where needed ---\n')
    mech_spc_dct = mechanalyzer.parser.spc.stereochemical_spc_dct(
        mech_spc_dct, nprocs='auto', all_stereo=True)
    print('Mechanism species with stereo added')
    for name, dct in mech_spc_dct.items():
        print(f'Name: {name:<25s} InChI: {dct["inchi"]}')

    # Expand the reactions in the mechanism to include stereochemical variants
    print('\n---- Expanding the list of mechanism reactions to include all'
          ' valid, stereoselective permutations ---\n')
    sccs_rxn_dct_lst = mechanalyzer.builder.expand_mech_stereo(
        rxn_param_dct, mech_spc_dct, nprocs='auto')

    # Loop over all sccs to write their mechfiles and keep track of
    # the good ones to reduce from
    ccs_sccs_spc_dct = {}
    for ccs_idx, sccs_rxn_dct in enumerate(sccs_rxn_dct_lst):
        for sccs_idx, sccs_rxn_lst in sccs_rxn_dct.items():

            ret_dcts = dictionaries_from_rxn_lst(
                sccs_rxn_lst)
            ste_mech_spc_dct, ste_rxn_dct = ret_dcts
            is_valid = mechanalyzer.builder.valid_enantiomerically(
                ste_mech_spc_dct)

            if is_valid:
                if ccs_idx not in ccs_sccs_spc_dct:
                    ccs_sccs_spc_dct[ccs_idx] = {0: [
                        ste_mech_spc_dct[spc_name][
                            'inchi'] for spc_name in ste_mech_spc_dct.keys()]}
                else:
                    ccs_sccs_spc_dct[ccs_idx][sccs_idx] = [
                        ste_mech_spc_dct[spc_name][
                            'inchi'] for spc_name in ste_mech_spc_dct.keys()]
            # OBTAIN SORTED SPECIES AND MECHANISMS
            write_mechanism(
                ste_mech_spc_dct, ste_rxn_dct, out_loc, out_spc, out_mech,
                sort_lst.copy(), isolate_spc, ccs_idx, sccs_idx)

    # Choose the sccs for each ccs that works with the rest of the mech
    chosen_idx_lst = ()
    all_chosen_ichs = ()
    for ccs_idx, sccs_dct in ccs_sccs_spc_dct.items():
        if not chosen_idx_lst:
            chosen_idx_lst += ((ccs_idx, 0),)
            all_chosen_ichs += tuple(sccs_dct[0])
            continue
        max_overlap = 0
        best_idx = 0
        for sccs_idx, spc_dct in sccs_dct.items():
            num_overlap = len(set(all_chosen_ichs) & set(spc_dct))
            if num_overlap > max_overlap:
                max_overlap = num_overlap
                best_idx = sccs_idx
        chosen_idx_lst += ((ccs_idx, best_idx),)
        all_chosen_ichs += tuple(sccs_dct[best_idx])

    print('Identifying diastereomer abstractions')
    dias_idxs = mechanalyzer.builder.diastereomer_abstractions(
        sccs_rxn_dct_lst, ccs_sccs_spc_dct,
        chosen_idx_lst, all_chosen_ichs)

    if dias_idxs:
        dia_str = ', '.join(
            (f'({ccs}, {sccs})' for (ccs, sccs) in dias_idxs)
        )
        print(f'Found new diastereomers to add: {dia_str}')

        chosen_idx_lst += dias_idxs

    print('Reduced reactions are from (CCS,SCCS):', chosen_idx_lst)
    reduced_rxn_lst = ()
    for ccs_idx, sccs_idx in chosen_idx_lst:
        reduced_rxn_lst += sccs_rxn_dct_lst[ccs_idx][sccs_idx]

    ste_mech_spc_dct, ste_rxn_dct = dictionaries_from_rxn_lst(
        reduced_rxn_lst)
    write_mechanism(
        ste_mech_spc_dct, ste_rxn_dct, out_loc, out_spc, out_mech,
        sort_lst.copy(), isolate_spc)


def write_mechanism(
        ste_mech_spc_dct, ste_rxn_dct, out_loc, out_spc, out_mech,
        sort_lst, isolate_spc, ccs_idx=None, sccs_idx=None):
    """ write mechanism
    """
    # Write the new (sorted) species dictionary to a string
    headers = ('smiles', 'inchi', 'inchikey', 'mult', 'charge')
    ste_mech_spc_dct_sort = mechanalyzer.parser.spc.reorder_by_atomcount(
        ste_mech_spc_dct)
    csv_str = mechanalyzer.parser.spc.csv_string(
        ste_mech_spc_dct_sort, headers)

    # Write initial string to call the sorter
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=None,
        mech_spc_dct=ste_mech_spc_dct_sort,
        spc_nasa7_dct=None,
        rxn_param_dct=ste_rxn_dct,
        rxn_cmts_dct=None)

    # Use strings to generate ordered objects
    param_dct_sort, ste_mech_spc_dct_sort, cmts_dct, _, _ = sorter.sorted_mech(
        csv_str, mech_str, isolate_spc, sort_lst)
    rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
        rxn_sort_dct=cmts_dct)

    # Write the dictionaries to ordered strings
    header_lst = list(file_dct['headers'])
    csv_str = mechanalyzer.parser.spc.csv_string(
        ste_mech_spc_dct_sort, header_lst)
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=(),
        mech_spc_dct=ste_mech_spc_dct_sort,
        spc_nasa7_dct=None,
        rxn_param_dct=param_dct_sort,
        rxn_cmts_dct=rxn_cmts_dct)

    # Write the species and mechanism files
    if ccs_idx is not None and sccs_idx is not None:
        ioformat.pathtools.write_file(
            csv_str, out_loc, out_spc + '_{:g}_{:g}'.format(ccs_idx, sccs_idx))
        ioformat.pathtools.write_file(
            mech_str, out_loc,
            out_mech + '_{:g}_{:g}'.format(ccs_idx, sccs_idx))
    else:
        ioformat.pathtools.write_file(csv_str, out_loc, out_spc)
        ioformat.pathtools.write_file(mech_str, out_loc, out_mech)


def input_from_location_dictionary(cwd, loc_dct):
    """ read mechanism info using file dictionary
        and current working directory
    """
    return input_info_from_file(
        cwd, loc_dct['inp_spc'], loc_dct['inp_mech'], loc_dct['sort'])


def input_info_from_file(cwd, species_file, mech_file, sort_file):
    """ read mechanism info using filenames
    """
    inp_spc_str = ioformat.pathtools.read_file(
        cwd, species_file,
        remove_comments='!', remove_whitespace=True)
    inp_mech_str = ioformat.pathtools.read_file(
        cwd, mech_file,
        remove_comments='#', remove_whitespace=True)
    sort_str = ioformat.pathtools.read_file(
        cwd, sort_file,
        remove_comments='#', remove_whitespace=True)
    return inp_spc_str, inp_mech_str, sort_str


def dictionaries_from_rxn_lst(sccs_rxn_lst):
    """ transform the reaction list to the dictionaries the writer likes
    """
    ste_mech_spc_dct, ste_rxn_dct = {}, {}
    ste_mech_spc_dct = mechanalyzer.builder.update_spc_dct_from_reactions(
        sccs_rxn_lst, ste_mech_spc_dct)
    ste_rxn_dct = mechanalyzer.builder.update_rxn_dct(
        sccs_rxn_lst, ste_rxn_dct, ste_mech_spc_dct)
    ste_mech_spc_dct = mechanalyzer.builder.remove_spc_not_in_reactions(
        ste_rxn_dct, ste_mech_spc_dct)
    return ste_mech_spc_dct, ste_rxn_dct


if __name__ == '__main__':

    # Initialize the start time for script execution
    t0 = time.time()

    # Set useful global variables
    oscwd = os.getcwd()

    # READ AND PARSE INPUT
    print('\n---- Parsing the species and mechanism files ---\n')

    # Read the build input file and set up info
    bld_str = ioformat.pathtools.read_file(
        oscwd, 'build.dat',
        remove_comments='#', remove_whitespace=True)
    file_dct, _ = mechanalyzer.parser.build_input_file(bld_str)

    # Read input species and mechanism files into dictionary
    mech_info = input_from_location_dictionary(oscwd, file_dct)
    main(oscwd, file_dct['out_spc'], file_dct['out_mech'], *mech_info)

    # Compute script run time and print to screen
    tf = time.time()
    print('\n\nScript executed successfully.')
    print(f'Time to complete: {tf-t0:.2f}')
    print('Exiting...')
