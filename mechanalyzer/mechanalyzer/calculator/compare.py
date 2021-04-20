"""
  Take data dictionaries from mechanisms and combine them under a common index
"""
#import mechlib.amech_io.parser.mechanism as parser
import ioformat.pathtools as parser
import itertools
import numpy as np
import ratefit
from phydat import phycon
from mechanalyzer.parser import spc as parser_spc
from chemkin_io.parser import mechanism as parser_mech
from chemkin_io.parser import thermo as parser_thermo
from chemkin_io.parser import reaction as parser_rxn
from chemkin_io.writer import _util as writer_util
from mechanalyzer.calculator import thermo as calc_thermo
from mechanalyzer.calculator import rates as calc_rates
import copy

RC_CAL = phycon.RC_cal  # universal gas constant in cal/mol-K


def get_aligned_rxn_ktp_dct(rxn_ktp_dcts, spc_thermo_dcts, spc_ident_dcts, temps,
                            rev_rates=True, remove_loners=True, write_file=False):
    """ Create an aligned_rxn_ktp_dct, which contains a single set of rxn keys, where each reaction
        has a ktp_dct for each mechanism. Options exist for reversing rates or not, removing rates
        that don't have entries from all mechanisms, and writing outputs to a file.
    
        :param rxn_ktp_dcts: list of rxn_ktp_dcts
        :type: list of dcts [rxn_ktp_dct1, rxn_ktp_dct2, ...]
        :param spc_thermo_dcts: list of spc_thermo_dcts
        :type: list of dcts [spc_thermo_dct1, spc_thermo_dct2, ...]
        :param spc_ident_dcts: list of spc_ident_dcts
        :type: list of dcts [spc_ident_dct1, spc_ident_dct2, ...]
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :param rev_rates: whether or not rates should be reversed
        :type rev_rates: Bool
        :param remove_loners: whether or not reactions with any None entries should be removed
        :type remove_loners: Bool
        :param write_file: whether or not to write output to a text file
        :type write_file: Bool
        :return aligned_rxn_ktp_dct: single dct with a list of ktp_dcts for each rxn
        :rtype: dct {rxn1: [ktp_dct1, ktp_dct2, ...], rxn2: ...}
    """
    assert len(rxn_ktp_dcts) == len(spc_thermo_dcts) == len(spc_ident_dcts), (
        f'Lengths of rxn_ktp_dcts ({len(rxn_ktp_dcts)}), spc_thermo_dcts ({len(spc_thermo_dcts)})' +
        f', and spc_ident_dcts ({len(spc_ident_dcts)}) should all be the same.'
    )

    # Get the renamed dictionaries
    renamed_rxn_ktp_dcts, rename_instructions_lst = rename_dcts(rxn_ktp_dcts, spc_ident_dcts,
                                                                target_type='rxn')
    if rev_rates:
        renamed_spc_thermo_dcts, _ = rename_dcts(spc_thermo_dcts, spc_ident_dcts,
                                                 target_type='spc')
    else:
        renamed_spc_thermo_dcts = []

    # Get the reversed_rxn_ktp_dcts
    reversed_rxn_ktp_dcts = reverse_rxn_ktp_dcts(renamed_rxn_ktp_dcts, renamed_spc_thermo_dcts,
                                                 temps, rev_rates=rev_rates)

    #print('inside get_aligned_rxn')
    #print('original:\n', rxn_ktp_dcts)
    #print('renamed:\n', renamed_rxn_ktp_dcts)
    #print('reversed:\n', reversed_rxn_ktp_dcts)

    # Get the aligned_rxn_ktp_dct
    aligned_rxn_ktp_dct = align_dcts(reversed_rxn_ktp_dcts)

    #print('aligned_rxn_ktp_dct:\n', aligned_rxn_ktp_dct)

    # If indicated, remove rxns that don't have rates from all the mechanisms
    if remove_loners:
        aligned_rxn_ktp_dct = remove_incomplete_items(aligned_rxn_ktp_dct)

    return aligned_rxn_ktp_dct


def get_aligned_spc_thermo_dct(spc_thermo_dcts, spc_ident_dcts, remove_loners=True,
                               write_file=False, print_output=False):
    """ Create an aligned_spc_thermo_dct, which contains a single set of spc keys, where each 
        species has a thermo_array for each mechanism. Options exist for removing species
        that don't have entries from all mechanisms and writing outputs to a file.

        :param spc_thermo_dcts: list of spc_thermo_dcts
        :type: list of dcts [spc_thermo_dct1, spc_thermo_dct2, ...]
        :param spc_ident_dcts: list of spc_ident_dcts
        :type: list of dcts [spc_ident_dct1, spc_ident_dct2, ...]
        :param remove_loners: whether or not species with any None entries should be removed
        :type remove_loners: Bool
        :param write_file: whether or not to write output to a text file
        :type write_file: Bool
        :return aligned_spc_thermo_dct: single dct with a list of thermo_arrays for each spc
        :rtype: dct {spc1: [thermo_array1, thermo_array2, ...], spc2: ...}
    """
    assert len(spc_thermo_dcts) == len(spc_ident_dcts), (
        f'Lengths of spc_thermo_dcts ({len(spc_thermo_dcts)}) ' +
        f'and spc_ident_dcts ({len(spc_ident_dcts)}) should be the same.'
        )

    # Get the renamed spc_thermo_dct
    renamed_spc_thermo_dcts, rename_instructions_lst = rename_dcts(
        spc_thermo_dcts, spc_ident_dcts, target_type='spc'
    )

    # Get the aligned_rxn_ktp_dct
    aligned_spc_thermo_dct = align_dcts(renamed_spc_thermo_dcts)

    # If indicated, remove spcs that don't have thermo from all the mechanisms
    if remove_loners:
        aligned_spc_thermo_dct = remove_incomplete_items(aligned_spc_thermo_dct)

    return aligned_spc_thermo_dct


def align_dcts(list_of_dcts):
    """ Takes a list of dcts and aligns it to be a single dct. Can handle any list of dcts in
        principle, but intended for either a reversed_rxn_ktp_dct or a renamed_spc_thermo_dct.

        :param list_of_dcts: a list of dcts
        :type list_of_dcts: [dct1, dct2, ...]
        :return aligned_dct: single dct with a list of items in each value
        :rtype: dct {key1: [item1, item2, ...], key2: ...}
    """
    num_mechs = len(list_of_dcts)
    aligned_dct = {}
    for mech_idx, dct in enumerate(list_of_dcts):
        for key, value in dct.items():

            # If the key does not yet exist, add it
            if key not in aligned_dct.keys():
                dct_list = [None] * mech_idx  # add empty entries to account for previous mechs that are missing
                dct_list.append(value)
                aligned_dct[key] = dct_list

            # If the key already exists, append the new dct
            else:
                dct_list = aligned_dct[key]  # get the current list of dcts
                # If any of the previous entries were blank, add None entries to fill
                if len(dct_list) < mech_idx:
                    dct_list.extend([None] * (mech_idx - len(dct_list)))
                dct_list.append(value)

    # Fill up the aligned dct by adding Nones for missing mechs
    for key, dct_list in aligned_dct.items():
        if len(dct_list) < num_mechs:
            dct_list.extend([None] * (num_mechs - len(dct_list)))
            aligned_dct[key] = dct_list

    return aligned_dct


def write_output_file(rename_instructions_lst, aligned_rxn_ktp_dct=None):
    """ Write the output of a comparison to a text file
    """
    with open("comparison_output.txt", "w") as f:
        for mech_idx, rename_instructions in enumerate(rename_instructions_lst):
            f.write(f"Rename instructions for converting mech {mech_idx+2} to mech {mech_idx+1}:\n")
            f.write(f"First column: mech {mech_idx+2} name, Second column: mech {mech_idx+1} name\n")
            for spc_name2, spc_name1 in rename_instructions.items():
                f.write(('{0:<1s}, {1:<5s}\n').format(spc_name2, spc_name1))
            f.write("\n\n")

    # NOTICE TO SELF: Commenting out for now...need to redo when I finally fix the handling of third bodies
#    if aligned_rxn_ktp_dct:
#        f.write(f"Combined_rxn_ktp_dct\n\n")
#        f.write(('{0:<64s}{1:<15s}{2:<15s}').format('Rxn name', 'In mech 1?', 'In mech 2?'))
#        for rxn, ktp_dct_lst in aligned_rxn_ktp_dct.items():
#            rxn_name = format_rxn_name(rxn, aligned_rxn_em_dct[rxn])  # note: call this from chemkinio
#            present = []
#            for entry in ktp_dct_lst:
#                if entry:
#                    present.append('Yes')
#                else: 
#                    present.append('No')
#
#            f.write(('\n{0:<64s}{1:<15s}{2:<15s}').format(rxn_name, present[0], present[1]))
       
    f.close()


def remove_incomplete_items(aligned_dct):
    """ Takes an aligned_dct and removes any entries that don't have values from
        all the mechs (i.e., any entries that have any None entries).
    
        :param aligned_dct: aligned dct with all keys (rxns or spcs) written in the same way 
        :type aligned_dct: dct {key1: [val1, val2, ...], key2: ...}
        :return filtered_aligned_dct: aligned dct with only entries that have no None values
        :rtype: dct {key1: [val1, val2, ...], key2: ...}
    """
    filtered_aligned_dct = {}
    for key, values in aligned_dct.items():
        num_values = len(values)
        num_actual_values = len([value for value in values if value is not None])

        # Only add the spc if all mechanisms have an entry for the spc
        if num_actual_values == num_values:
            filtered_aligned_dct[key] = values

    return filtered_aligned_dct


def rename_dcts(target_dcts, spc_ident_dcts, target_type):
    """ Takes a list of dictionaries and renames all the species. The species are renamed in order
        of the preference specified by the order of the list (first dct is unchanged, second is
        only changed by first, third is changed by first and then second, etc.).

        :param target_dcts: list of dictionaries to be renamed
        :type target_dcts: list of dcts; [{dct1}, {dct2}, ...]
        :param spc_ident_dcts: list of species dictionaries corresponding to dcts
        :type spc_ident_dcts: list of dcts; [{dct1}, {dct2}, ...]
        :param target_type: either 'rxn' or 'spc'; refers to what the key of the dct is
        :type target_type: str
        :return renamed_dcts: dcts with species renamed
        :rtype: list of dcts [dct1, dct2, ...]
    """
    renamed_target_dcts = copy.deepcopy(target_dcts)  # deepcopy to prevent external changes
    renamed_spc_ident_dcts = copy.deepcopy(spc_ident_dcts)
    num_mechs = len(target_dcts)
    assert num_mechs == len(spc_ident_dcts), (
        f'Length of dct_list is {num_mechs} while length of spc_dct_lst is {len(spc_ident_dcts)}.'
        )

    # Loop through each item in the list of dictionaries
    rename_instructions_lst = []
    for mech_idx in range(num_mechs-1):
        spc_ident_dct1 = renamed_spc_ident_dcts[mech_idx]
        for idx2 in range(mech_idx+1, num_mechs):
            spc_ident_dct2 = renamed_spc_ident_dcts[idx2]

            # Get the rename instructions from the species dcts
            rename_instructions, _ = get_rename_instructions(spc_ident_dct1, spc_ident_dct2)

            # Rename and store the current spc_dct
            renamed_spc_ident_dct2 = rename_species(spc_ident_dct2, rename_instructions, target_type='spc')
            renamed_spc_ident_dcts[idx2] = renamed_spc_ident_dct2

            # Rename and store the current dct
            renamed_dct = rename_species(renamed_target_dcts[idx2], rename_instructions, target_type)
            renamed_target_dcts[idx2] = renamed_dct
            rename_instructions_lst.append(rename_instructions)
    
    return renamed_target_dcts, rename_instructions_lst


def get_rename_instructions(spc_ident_dct1, spc_ident_dct2):
    """ Get instructions for renaming spc_ident_dct2 to be consistent with spc_ident_dct1. Also,
        output a combined spc_ident_dct that contains all species from spc_ident_dct1 along with
        any species unique to spc_ident_dct2
        
        :param spc_ident_dct1: the reference spc_ident_dct
        :type spc_ident_dct1: dct {spc1: ident_array1, spc2: ...}
        :param spc_ident_dct2: the spc_ident_dct to be renamed (and also added)
        :type spc_ident_dct2: dct {spc1: ident_array1, spc2: ...}
        :return rename_instructions: instructions for renaming the species in spc_ident_dct2
        :rtype: dct {spc_to_be_renamed1: new_spc_name1, spc_to_be_renamed2: ...}
        :return combined_spc_dct: spc_ident_dct1 plus any species unique to spc_ident_dct2
        :rtype: dct {spc1: ident_array1, spc2: ...}
    """
    combined_spc_dct = copy.deepcopy(spc_ident_dct1)  # deepcopy to prevent external changes
    rename_instructions = {}
    rename_str = '-zz'
    
    # Loop through each species in mech1
    for spc_name1, spc_vals1 in spc_ident_dct1.items():
        i1 = spc_vals1['inchi']
        m1 = spc_vals1['mult']
        c1 = spc_vals1['charge']

        for spc_name2, spc_vals2 in spc_ident_dct2.items():
            i2 = spc_vals2['inchi']
            m2 = spc_vals2['mult']
            c2 = spc_vals2['charge']
             
            # If species are identical      
            if i1  == i2 and m1 == m2 and c1 == c2:
                if spc_name1 == spc_name2: # do nothing
                    break
                else:
                    rename_instructions[spc_name2] = spc_name1

            # If species are different but have same name
            elif spc_name1 == spc_name2:  
                rename_instructions[spc_name2] = spc_name2 + rename_str                 

    # Add species that are only in mech2
    for spc_name2, spc_vals2 in spc_ident_dct2.items():
        unique = True
        i2 = spc_vals2['inchi']
        m2 = spc_vals2['mult']
        c2 = spc_vals2['charge']

        for spc_name1, spc_vals1 in spc_ident_dct1.items():
            i1 = spc_vals1['inchi']
            m1 = spc_vals1['mult']
            c1 = spc_vals1['charge']

            if i1 == i2 and m1 == m2 and c1 == c2:
                if spc_name2 in rename_instructions.keys():
                    unique = False
                break
                
        if unique:
            if spc_name2 in rename_instructions.keys():
                combined_spc_dct[spc_name2 + rename_str] = spc_vals2
            else:
                combined_spc_dct[spc_name2] = spc_vals2

    return rename_instructions, combined_spc_dct


def rename_species(target_dct, rename_instructions, target_type='rxn'):
    """ Rename the species inside a rxn_ktp_dct, rxn_param_dct, or thermo_dct according to the
        instructions inside the rename_instructions.

        :param target_dct: the dct whose species are to be renamed
        :type target_dct: dct; either a rxn_ktp, rxn_param, or thermo dct
        :param rename_instructions: instructions for renaming the species in the target_dct
        :type rename_instructions: dct {spc_to_be_renamed1: new_spc_name1, spc_to_be_renamed2: ...}
        :param target_type: the type of dictionary; either "rxn", "thermo", or "spc"
        :type target_type: str
        :return renamed_dct: dct with all species renamed according to the rename_instructions
        :rtype: dct; either a rxn_ktp, rxn_param, or thermo dct
    """
    def strip_third_bod(third_bod):
        """ Strip away the '(', '+', and ')' from a third body
        """
        if third_bod == '(+M)' or third_bod == '+M' or third_bod is None:
            stripped_third_bod = third_bod  # return unchanged third_bod
            addition = None
        elif third_bod[0] == '(':
            stripped_third_bod = third_bod[2:(len(third_bod) - 1)]  # get rid of '(', '+', and ')'
            addition = 'paren'
        else:
            stripped_third_bod = third_bod[1:]  # just get rid of the '+'
            addition = 'plus only'

        return stripped_third_bod, addition

    assert target_type in ('rxn', 'thermo', 'spc'), (
        f'The target_type is {target_type}, but should be either "rxn", "thermo", or "spc"'
        )
    renamed_dct = {}

    # If a rxn_ktp_dct
    if target_type == 'rxn':
        for rcts, prds, third_bods in target_dct.keys():
            new_rcts = []
            new_prds = []
            new_third_bods = []
            for spc in rcts:
                if spc in rename_instructions.keys():
                    new_rcts.append(rename_instructions[spc])
                else:
                    new_rcts.append(spc)
            for spc in prds:
                if spc in rename_instructions.keys():
                    new_prds.append(rename_instructions[spc])
                else:
                    new_prds.append(spc)
            for spc in third_bods:
                spc, addition = strip_third_bod(spc)  # if not '(+M)' or '+M', strip '(+)'
                if spc in rename_instructions.keys():
                    new_third_bod = rename_instructions[spc]
                else:  # this condition will occur if third_bod is '(+M)' or '+M'
                    new_third_bod = spc
                
                # Add back on the '(+...)' or '+' if indicated
                if addition == 'paren':
                    new_third_bod = f'(+{new_third_bod})'
                elif addition == 'plus only':
                    new_third_bod = f'+{new_third_bod}'
                new_third_bods.append(new_third_bod)
    
            new_rcts = tuple(new_rcts)
            new_prds = tuple(new_prds)
            new_third_bods = tuple(new_third_bods)
            renamed_dct[new_rcts, new_prds, new_third_bods] = target_dct[rcts, prds, third_bods]

    # If a spc_thermo_dct or spc_ident_dct
    else:
        for spc, data in target_dct.items():
            if spc in rename_instructions.keys():
                new_spc_name = rename_instructions[spc]
                renamed_dct[new_spc_name] = data
            else:
                renamed_dct[spc] = data

    return renamed_dct


def reverse_rxn_ktp_dcts(renamed_rxn_ktp_dcts, renamed_spc_thermo_dcts, temps, rev_rates=True):
    """ Takes a list of rxn_ktp_dcts that have all been renamed to have consistent species
        names and reverses any reactions if requested (if rev_rates=True) and also reorders all
        reactant/product ordering to be consistent.

        :param renamed_rxn_ktp_dcts: list of rxn_ktp_dcts with all species renamed to be matching
        :type renamed_rxn_ktp_dcts: list of dcts [rxn_ktp_dct1, rxn_ktp_dct2, ...]
        :param renamed_spc_thermo_dcts: list of spc_thermo_dcts with all  species renamed to be matching
        :type renamed_spc_thermo_dcts: list of dcts [spc_thermo_dct1, spc_thermo_dct2, ...]
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :param rev_rates: whether or not rates should be reversed
        :type rev_rates: Bool
        :return reversed_rxn_ktp_dcts: list of rxn_ktp_dcts that have the rxns written with reacts and
            prods in the same order (and reversed if that was indicated)
        :rtype: list of dcts [rxn_ktp_dct1, rxn_ktp_dct2, ...]
    """
    #renamed_rxn_ktp_dcts = copy.deepcopy(renamed_rxn_ktp_dcts)  # deepcopy to prevent external changes
    num_mechs = len(renamed_rxn_ktp_dcts)
    reversed_rxn_ktp_dcts = copy.deepcopy(renamed_rxn_ktp_dcts)  # deepcopy so no external changes
    #reversed_rxn_ktp_dcts = renamed_rxn_ktp_dcts
    for mech_idx in range(num_mechs-1):
        rxn_ktp_dct1 = renamed_rxn_ktp_dcts[mech_idx]
        for idx2 in range(mech_idx+1, num_mechs):
            rxn_ktp_dct2 = renamed_rxn_ktp_dcts[idx2]

            # If reversing rates, get thermo; otherwise, just get a blank list
            if rev_rates:
                spc_thermo_dct1 = renamed_spc_thermo_dcts[mech_idx]
            else:
                spc_thermo_dct1 = []

            # Either reverse rxns (if rev_rates=True) or just flip the products and reactants to
            # be in the same order
            reversed_rxn_ktp_dct = reverse_rxn_ktp_dct(
                rxn_ktp_dct1, rxn_ktp_dct2, spc_thermo_dct1, temps, rev_rates=rev_rates
            )
            reversed_rxn_ktp_dcts[idx2] = reversed_rxn_ktp_dct

    return reversed_rxn_ktp_dcts


def reverse_rxn_ktp_dct(rxn_ktp_dct1, rxn_ktp_dct2, spc_thermo_dct1, temps, rev_rates=True):
    """ Takes two rxn_ktp_dcts whose species have already been renamed to be identical and reverses
        any reactions *in the second dct* that need to be reversed

        :param rxn_ktp_dct1: rxn_ktp_dct for mech1
        :type rxn_ktp_dct1: dict {rxn1: ktp_dct1, rxn2: ...}
        :param rxn_ktp_dct2: rxn_ktp_dct for mech2
        :type rxn_ktp_dct2: dict {rxn1: ktp_dct1, rxn2: ...}
        :param spc_thermo_dct1: spc_thermo_dct for mech1
        :type spc_thermo_dct1: dict {spc1: thermo_array1, spc2: ...}
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :param rev_rates: whether or not rates should be reversed
        :type rev_rates: Bool
    """
#    rxn_ktp_dct1 = copy.deepcopy(rxn_ktp_dct1)  # deepcopy to prevent external changes
#    rxn_ktp_dct2 = copy.deepcopy(rxn_ktp_dct2)  # deepcopy to prevent external changes
    rev_rxn_ktp_dct2 = copy.deepcopy(rxn_ktp_dct2)  # deepcopy to prevent external changes
    for rxn1, ktp_dct1 in rxn_ktp_dct1.items():
        rxn2, rev_rate = assess_rxn_match(rxn1, rxn_ktp_dct2)
        # Only do something if a match was found
        if rxn2 is not None:
            # If the user indicated to reverse rates, check if they need to be
            if rev_rates:
                if rev_rate:
                    ktp_dct2 = rxn_ktp_dct2[rxn2]
                    rev_ktp_dct2 = reverse_ktp_dct(ktp_dct2, spc_thermo_dct1, rxn2, temps)
                    rev_rxn_ktp_dct2.pop(rxn2)
                    rev_rxn_ktp_dct2[rxn1] = rev_ktp_dct2
                # If a match was found that does not need to be reversed but has rcts and prds written
                # differently, align rcts and prds
                elif rxn1 != rxn2:
                    rev_rxn_ktp_dct2[rxn1] = rev_rxn_ktp_dct2[rxn2]
                    rev_rxn_ktp_dct2.pop(rxn2)
    
            # Otherwise, rename rxn2 so that rcts and prds are in the same order if they're not already
            # Only do so if the rxns should not be flipped!
            elif rxn1 != rxn2 and not rev_rate: 
                    rev_rxn_ktp_dct2[rxn1] = rev_rxn_ktp_dct2[rxn2]
                    rev_rxn_ktp_dct2.pop(rxn2)

    return rev_rxn_ktp_dct2


def reverse_ktp_dct(ktp_dct, spc_thermo_dct, rxn, temps):
    """ For a given reaction, use the thermochemistry values of its constituent species to
        calculate the equilibrium constant and reverse the rate constants.

        :param ktp_dct: k(T,P) dct for a single reaction
        :type ktp_dct: dict {pressure1: (temp_array1, rates_array1), pressure2: ...}
        :param spc_thermo_dct: thermochemical values for all species in mechanism
        :type spc_thermo_dct: dict {spc1: thermo_array1, spc2: ...}
        :param rxn: rxn key
        :type rxn: tuple (rcts, prds, third_bods)
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :return rev_ktp_dct: k(T,P) dct of the reversed reaction
        :type ktp_dct: dict {pressure1: (temp_array1, rates_array1), pressure2: ...}
    """
    [rcts, prds, third_bods] = rxn
    k_equils = _calculate_equilibrium_constant(spc_thermo_dct, rcts, prds, temps)
    rev_ktp_dct = {}
    for pressure, (_, kts) in ktp_dct.items():

        # Calculate density to handle units, if needed
        if len(rcts) > 1 and len(prds) == 1:
            densities = ratefit.calc.p_to_m(1.0, temps)
            kts *= densities
        elif len(rcts) == 1 and len(prds) > 1:
            densities = ratefit.calc.p_to_m(1.0, temps)
            kts /= densities

        # Calculate the reverse rates with K_equil
        rev_rates = []
        for forw_k, k_equil in zip(kts, k_equils):
            rev_rates.append(forw_k / k_equil)

        # Add reversed rates to dict
        rev_ktp_dct[pressure] = (temps, rev_rates)

    return rev_ktp_dct


def assess_rxn_match(rxn1, rxn_ktp_dct2):
    """ Assess whether the reaction should be flipped. Takes a rxn_name from mech1 and searches
        through all of mech2 in search of a matching rxn. If a matching rxn is found, returns the
        matching rxn name and whether the rxn should be flipped

        Note: it is possible that a poorly constructed mechanism will have more than one instance
        of the same reaction. This function will only return the first instance of any matching
        reaction. However, it will print out a warning if duplicate matching reactions are
        found.

        :param rxn1: rxn key for which a match is being sought
        :type rxn1: tuple (rcts, prds, third_bods)
        :param rxn_ktp_dct2: rxn_ktp_dct for mech2
        :type rxn_ktp_dct2: dict {rxn1: ktp_dct1, rxn2: ...}
        :return matching_rxn: rxn key for the matching reaction
        :rtype: tuple (rcts, prds, third_bods)
        :return rev_rate: whether or not the rate should be reversed
        :rtype: Bool
    """
    # Get all possible orderings of the reactants and products for mech1
    [rcts1, prds1, third_bods1] = rxn1
    third_bod1 = third_bods1[0]
    rcts1_perm = list(itertools.permutations(rcts1, len(rcts1)))
    prds1_perm = list(itertools.permutations(prds1, len(prds1)))

    matching_rxn_name = None
    rev_rate = None
    already_found = False
    for rxn2 in rxn_ktp_dct2.keys():
        [rcts2, prds2, third_bods2] = rxn2
        third_bod2 = third_bods2[0]
        if rcts2 in rcts1_perm and prds2 in prds1_perm and third_bod1 == third_bod2:
            matching_rxn_name = rxn2
            rev_rate = False
            if already_found:
                rxn_name1 = writer_util.format_rxn_name(rxn1)
                rxn_name2 = writer_util.format_rxn_name(rxn2)
                print(f'For the reaction {rxn_name1}, more than one match was found: {rxn_name2}')
                print('This will cause errors!')
            else:
                already_found = True

        if rcts2 in prds1_perm and prds2 in rcts1_perm and third_bod1 == third_bod2:
            matching_rxn_name = rxn2
            rev_rate = True
            if already_found:
                rxn_name1 = writer_util.format_rxn_name(rxn1)
                rxn_name2 = writer_util.format_rxn_name(rxn2)
                print(f'For the reaction {rxn_name1}, more than one match was found: {rxn_name2}')
                print('This will cause errors!')
            else:
                already_found = True

    return matching_rxn_name, rev_rate


def _calculate_equilibrium_constant(spc_thermo_dct, rcts, prds, temps):
    """ Calculate the equilibrium constant for a given reaction at
        a set of temperatures using constituent species' thermochemistry.

        :param spc_thermo_dct: thermochemical values for all species in mechanism
        :type spc_thermo_dct: dict {spc1: thermo_array1, spc2: ...}
        :param rcts: name(s) of reactant(s)
        :type rcts: tuple (rct1, rct2, ...)
        :param prds: name(s) of product(s)
        :type prds: tuple (prd1, prd2, ...)
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :return k_equils: equilibrium constant at each temperature
        :rtype: list [float]
    """
    k_equils = []
    for temp_idx, temp in enumerate(temps):
        rct_gibbs = 0.0
        for rct in rcts:
            rct_gibbs += spc_thermo_dct[rct][4][temp_idx]  # [4] accesses Gibbs

        prd_gibbs = 0.0
        for prd in prds:
            prd_gibbs += spc_thermo_dct[prd][4][temp_idx]  # [4] accesses Gibbs

        rxn_gibbs = prd_gibbs - rct_gibbs
        k_equils.append(np.exp(-rxn_gibbs / (RC_CAL * temp)))

    return k_equils


def load_rxn_ktp_dcts_chemkin(mech_filenames, direc, temps, pressures):
    """ Read one or more Chemkin-formatted mechanisms files and calculate rates at the indicated
        pressures and temperatures. Return a list of rxn_ktp_dcts.

        :param mech_filenames: filenames containing Chemkin-formatted kinetics information
        :type mech_filenames: list [filename1, filename2, ...]
        :param direc: directory with file(s) (all files must be in the same directory)
        :type direc: str
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :param pressures: pressures at which to do calculations (atm)
        :type pressures: list [float]
        :return rxn_ktp_dcts: list of rxn_ktp_dcts
        :rtype: list of dcts [rxn_ktp_dct1, rxn_ktp_dct2, ...]
    """
    rxn_ktp_dcts = []
    for idx in range(len(mech_filenames)):
        mech_str = parser.read_file(direc, mech_filenames[idx])
        ea_units, a_units = parser_mech.reaction_units(mech_str)
        rxn_block_str = parser_mech.reaction_block(mech_str)
        rxn_param_dct = parser_rxn.param_dct(rxn_block_str, ea_units, a_units)
        rxn_ktp_dct = calc_rates.eval_rxn_param_dct(rxn_param_dct, pressures, temps)
        rxn_ktp_dcts.append(rxn_ktp_dct)

    return rxn_ktp_dcts


def load_spc_thermo_dcts_chemkin(thermo_filenames, direc, temps):
    """ Reads one or more Chemkin-formatted thermo files and calculates thermo at the indicated
        temperatures. Outputs a list of spc_thermo_dcts.

        :param thermo_filenames: filenames containing Chemkin-formatted thermo information
        :type thermo_filenames: list [filename1, filename2, ...]
        :param direc: directory with file(s) (all files must be in the same directory)
        :type direc: str
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :return spc_thermo_dcts: list of spc_thermo_dcts
        :rtype: list of dcts [spc_thermo_dct1, spc_thermo_dct2, ...]
    """
    spc_thermo_dcts = []
    for idx in range(len(thermo_filenames)):
        thermo_str = parser.read_file(direc, thermo_filenames[idx])
        thermo_block_str = parser_mech.thermo_block(thermo_str)
        spc_nasa7_dct = parser_thermo.create_spc_nasa7_dct(thermo_block_str)
        spc_thermo_dct = calc_thermo.create_spc_thermo_dct(spc_nasa7_dct, temps)
        spc_thermo_dcts.append(spc_thermo_dct)

    return spc_thermo_dcts


def load_spc_ident_dcts(spc_csv_filenames, direc):
    """ Reads one or more spc.csv files (in the AutoMech format). Outputs a list of
        spc_ident_dcts.

        :param spc_csv_filenames: filenames containing species identity information (spc.csv files)
        :type spc_csv_filenames: list [filename1, filename2, ...]
        :param direc: directory with file(s) (all files must be in the same directory)
        :type direc: str
        :return spc_ident_dcts: list of spc_ident_dcts
        :rtype: list of dcts [spc_ident_dct1, spc_ident_dct2, ...]
    """
    spc_ident_dcts = []
    for idx in range(len(spc_csv_filenames)):
        spc_csv_str = parser.read_file(direc, spc_csv_filenames[idx])
        spc_ident_dct = parser_spc.build_spc_dct(spc_csv_str, 'csv')
        spc_ident_dcts.append(spc_ident_dct)

    return spc_ident_dcts


def combine_spc_nasa7(spc_nasa7_dct1, spc_nasa7_dct2, rename_instructions):
    """

    """
    combined_spc_nasa7_dct = copy.copy(spc_nasa7_dct1)
    renamed_spc_nasa7_dct2 = rename_species(spc_nasa7_dct2, rename_instructions, target_type='spc')

    for spc_name2, nasa7_params2 in renamed_spc_nasa7_dct2.items():
        if spc_name2 not in spc_nasa7_dct1.keys():
            combined_spc_nasa7_dct[spc_name2] = nasa7_params2

    return combined_spc_nasa7_dct


def combine_mech_params(rxn_param_dct1, rxn_param_dct2, rename_instructions):
    """ Combine together two rxn_param_dcts to create a single dct,
        renaming all species as indicated in the rename_instructions.
    """
    combined_rxn_param_dct = copy.copy(rxn_param_dct1)
    renamed_rxn_param_dct2 = rename_species(rxn_param_dct2, rename_instructions)

    for rxn_name2, param2 in renamed_rxn_param_dct2.items():
        rxn_name, flip_rxn, p_dep_same = assess_rxn_match(rxn_name2, param2, rxn_param_dct1)

        # If rxn is unique to param2, add to param 1. Note: this does not
        # consider whether the pressure dependence is the same. It just keeps
        # the rate from mechanism 1.
        if not rxn_name:
            combined_rxn_param_dct[rxn_name2] = param2

    return combined_rxn_param_dct
