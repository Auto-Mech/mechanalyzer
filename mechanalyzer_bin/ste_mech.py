""" Script to generate files with reactions
"""

import os
import time
import itertools as it

import ioformat
import chemkin_io
import mechanalyzer
from autofile import io_ as io
import automol
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


def read_all_sccs_spc(cwd, spc_file_name):
    """ create a nested dictionary with ccs, sccs, spc_lst
        where spc lst is the inchis of that ccs, sccs combo
    """
    files = os.listdir(cwd)
    spc_files = [fil for fil in files if spc_file_name + '_' in fil]

    spc_ccs_dct = {}
    for fil in spc_files:
        ccs, sccs = [int(idx) for idx in fil.split('_')[-2:]]
        spc_lines = io.read_file(fil).splitlines()
        ich_col = [
            idx for idx, val in enumerate(spc_lines[0].split("'"))
            if val == 'inchi'][0]
        spc_ichs = tuple([line.split("'")[ich_col] for line in spc_lines[1:]])
        if ccs not in spc_ccs_dct:
            spc_ccs_dct[ccs] = {sccs: spc_ichs}
        else:
            spc_ccs_dct[ccs][sccs] = spc_ichs
    return spc_ccs_dct


def group_sccs(spc_ccs_dct):
    """ using the keys of the spc_ccs_dct group all
        ccs combos into a tuple of ccs_sccs tuples
    """
    all_idxs = ()
    ccs_lst = sorted(list(spc_ccs_dct.keys()))
    for ccs in ccs_lst:
        sccs_lst = sorted(list(spc_ccs_dct[ccs].keys()))
        sccs_tup = tuple([(ccs, sccs,) for sccs in sccs_lst])
        all_idxs += (sccs_tup,)
    return all_idxs


def seperate_single_ccs(all_idxs):
    sccs_lsts = ()
    ccs_lsts = ()
    for idx_lst in all_idxs:
        if len(idx_lst) == 1:
            ccs_lsts += idx_lst
        else:
            sccs_lsts += (idx_lst,)
    return sccs_lsts, ccs_lsts


def all_sccs_combos(spc_ccs_dct):
    """ build all possible combinations of sccs
    """
    combo_lst = ((),)
    all_idxs = group_sccs(spc_ccs_dct)
    sccs_lsts, ccs_lsts = seperate_single_ccs(all_idxs)
    for ccs_tups in sccs_lsts:
        new_combo_lst = ()
        print(ccs_tups)
        for idxs in ccs_tups:
            for combo in combo_lst:
                new_combo_lst += (combo + (idxs,),)
        combo_lst = new_combo_lst

    new_combo_lst = ()
    for combo in combo_lst:
        new_combo_lst += (combo + ccs_lsts,)
    return new_combo_lst


def find_best_combination(spc_ccs_dct, combo_lst):
    """ find the sccs combo that has the shortest unique spc lst
        which means it has the least number of enantiomers
    """
    uniq_spc_lst_lst = ()
    num_uniq_spc_for_combo = ()
    print('Checking {:g} combinations'.format(len(combo_lst)))
    for combo in combo_lst:
        spc_ich_lst = ()
        for idxs in combo:
            spc_ich_lst += spc_ccs_dct[idxs[0]][idxs[1]]
        uniq_spc_lst = ()
        for spc in spc_ich_lst:
            if spc not in uniq_spc_lst:
                uniq_spc_lst += (spc,)
        uniq_spc_lst_lst += (uniq_spc_lst,)
        num_uniq_spc_for_combo += (len(uniq_spc_lst),)
        # enant_count = 0
        # for spc_a, spc_b in it.combinations(spc_ich_lst, 2):
        #     if automol.inchi.are_enantiomers(spc_a, spc_b):
        #         enant_count += 1
        # combo_ent_count += (enant_count,)
        # print('found {:g} enantiomers for this combo'.format(enant_count))
    # print(combo_ent_count)
    # min_enants = min(combo_ent_count)
    # print('minimum overlap', min_enants)
    # best_combo_idx = combo_ent_count.index(min_enants)
    # best_combo = combo_ent_count[best_combo_idx]
    # print('best combo', best_combo)
    shortest_spc_lst = min(num_uniq_spc_for_combo)
    print('greatest overlap', shortest_spc_lst)
    best_combo_idx = num_uniq_spc_for_combo.index(shortest_spc_lst)
    best_combo = combo_lst[best_combo_idx]
    print('best combo', best_combo)
    enant_count = 0
    for spc_a, spc_b in it.combinations(uniq_spc_lst_lst[best_combo_idx], 2):
        if automol.inchi.are_enantiomers(spc_a, spc_b):
            enant_count += 1
    print('found {:g} enantiomers for this combo'.format(enant_count))
    return best_combo


def write_best_combination(best_combo, loc_dct, cwd):
    """ make combined species and mechanism dictionaries
        from a set of sccs
    """
    mech_spc_dct = {}
    rxn_param_dct = {}
    for (ccs, sccs,) in best_combo:
        mech_file = '{}_{:g}_{:g}'.format(loc_dct['out_mech'], ccs, sccs)
        spc_file = '{}_{:g}_{:g}'.format(loc_dct['out_spc'], ccs, sccs)
        inp_info = input_info_from_file(
            cwd, spc_file, mech_file, loc_dct['sort'])
        next_spc_dct = mechanalyzer.parser.spc.build_spc_dct(
            inp_info[0], 'csv')
        next_rxn_param_dct = mechanalyzer.parser.mech.parse_mechanism(
            inp_info[1], 'chemkin')
        mech_spc_dct.update(next_spc_dct)
        rxn_param_dct.update(next_rxn_param_dct)
    isolate_spc, sort_lst = mechanalyzer.parser.mech.parse_sort(inp_info[2])
    write_mechanism(
        mech_spc_dct, rxn_param_dct,
        cwd, 'red_species.csv', 'red_mechanism.dat',
        sort_lst, isolate_spc)


def reduction(cwd, loc_dct):
    """ reduce by doing all possible combinations
        of sccss and figuring out
        which combination results in the shortest unique species list
    """
    spc_sccs_dct = read_all_sccs_spc(cwd, loc_dct['out_spc'])
    sccs_combo_lst = all_sccs_combos(spc_sccs_dct)
    best_combo = find_best_combination(spc_sccs_dct, sccs_combo_lst)
    write_best_combination(best_combo, loc_dct, cwd)


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
    # main(oscwd, file_dct['out_spc'], file_dct['out_mech'], *mech_info)
    reduction(oscwd, file_dct)
    # Compute script run time and print to screen
    tf = time.time()
    print('\n\nScript executed successfully.')
    print(f'Time to complete: {tf-t0:.2f}')
    print('Exiting...')
