"""
  Take data dictionaries from mechanisms and combine them under a common index
"""
import lib.amech_io.parser as parser
import itertools
import numpy as np
import ratefit
from ioformat import remove_whitespace
from ioformat import phycon
from mechanalyzer.parser import spc as parser_spc
from chemkin_io.parser import mechanism as parser_mech
from chemkin_io.parser import reaction as parser_rxn
from mechanalyzer.calculator import thermo as calc_thermo
from mechanalyzer.calculator import rates as calc_rates
from chemkin_io.parser.mechanism import reaction_units
import copy
import os
import sys

# Read the job path, which should have been input as the current directory
JOB_PATH = sys.argv[1]
    

def align_mechs(mech_filenames, thermo_filenames, spc_csv_filenames, temps, pressures, rev_rates=True):
    """ Takes any number of mechanisms and renames the species to all be the same,
        reverses any reactions to be in the same direction, and outputs a dictionary
        of all the reaction names and the corresponding k(T,P) dictionaries.

        :param mech_filenames: names of the Chemkin-formatted mechanism text files
        :type mech_filenames: list
        :param thermo_filenames: names of the Chemkin-formatted thermo text files; if the thermo
            info is contained in the mechanism.txt, just put the same filename here
        :type thermo_filenames: list
        :param spc_csv_filenames: names of the spc_csv filenames, which should have the following info:
            names, SMILES, inchi, multiplicity, charge, sensitivity (sensitivity not required)
        :type spc_csv_filenames: list
        :param temps: array of temperatures in Kelvin
        :type temps: Numpy 1-D array
        :param pressures: array of pressures in atm
        :type pressures: Numpy 1-D arrayi
        :param rev_rates: whether or not to reverse reactions; thermo filenames not required if False
        :type rev_rates: Bool
        :return combined_rxn_ktp_dcts: dictionary of k(T,P) values for each mechanism, all under the same
            reaction name index
        :rtype: dict {((rcts), (prds)): [{ktp_dct1, ktp_dct2, ...], ...}
        :return combined_rxn_em_dcts: dictionary of Boolean parameters indicating whether or not the rxn
            has a '+ M' term (i.e., third body) associated with it
        :rtype: dict {((rcts), (prds)): Boolean, ...}
    """
    if rev_rates:
        assert len(mech_filenames) == len(spc_csv_filenames) == len(thermo_filenames), (
            f"""The lengths of the mechanism, spc_csv, and thermo filename inputs should 
                all be the same, but are instead {len(mech_filenames)}, {len(spc_csv_filenames)}, 
                and {len(thermo_filenames)}, respectively."""
            )
    else: 
        assert len(mech_filenames) == len(spc_csv_filenames), (
            f"""The lengths of the mechanism and spc_csv inputs should 
                all be the same, but are instead {len(mech_filenames)} and {len(spc_csv_filenames)}, 
                respectively."""
            )
        
    num_mechs = len(mech_filenames)
    
    # Load rxn_ktp_dcts, spc_dcts, and thermo_dcts
    rxn_ktp_dcts = []
    rxn_param_dcts = []
    spc_dcts = []
    thermo_dcts = []
    for idx in range(num_mechs):
        # Get the rxn_ktp_dct and the rxn_param_dct
        mech_str = parser.ptt.read_inp_str(JOB_PATH,mech_filenames[idx],remove_comments=False)
        rxn_block_str = parser_mech.reaction_block(mech_str)
        if idx == 1:  # if on the second mechanism
            print('inside compare.py, changing units manually to kcal/mol')
            rxn_param_dct = parser_rxn.param_dct(rxn_block_str, ea_units='kcal/mol')
        else:
            rxn_param_dct = parser_rxn.param_dct(rxn_block_str)

        rxn_ktp_dct = calc_rates.eval_rxn_param_dct(rxn_param_dct, pressures, temps)
        rxn_ktp_dcts.append(rxn_ktp_dct)
        rxn_param_dcts.append(rxn_param_dct)

        # Get the spc_dct
        spc_csv_str = parser.ptt.read_inp_str(JOB_PATH, spc_csv_filenames[idx], remove_comments=False)
        spc_dct = parser_spc.build_spc_dct(spc_csv_str, 'csv')
        spc_dcts.append(spc_dct)

        # Get the thermo_dct
        if rev_rates:
            thermo_str = parser.ptt.read_inp_str(JOB_PATH,thermo_filenames[idx],remove_comments=False)
            thermo_block_str = parser_mech.thermo_block(thermo_str)
            thermo_dct = calc_thermo.mechanism(thermo_block_str, temps)
            thermo_dcts.append(thermo_dct)

#    print('rxn_ktp_dcts\n', rxn_ktp_dcts)
#    print('spc_dts\n', spc_dcts)
#    print('thermo_dcts\n', thermo_dcts)

    # Get the renamed dictionaries
    renamed_rxn_ktp_dcts = rename_dcts(rxn_ktp_dcts, 'rxn', spc_dcts)
    renamed_rxn_param_dcts = rename_dcts(rxn_param_dcts, 'rxn', spc_dcts)
    if rev_rates:
        renamed_thermo_dcts = rename_dcts(thermo_dcts, 'spc', spc_dcts)
#    print('renamed rxn_ktp_dcts\n', renamed_rxn_ktp_dcts) 
   
    # Get the em_param_dcts                    
    renamed_rxn_em_dcts = []
    for dct in renamed_rxn_param_dcts:
        renamed_rxn_em_dct = get_rxn_em_dct(dct)
        renamed_rxn_em_dcts.append(renamed_rxn_em_dct)
    
#    print('param_dcts\n', renamed_rxn_param_dcts)
#    print('em_dcts\n', renamed_rxn_em_dcts)

    # Loop over each reaction and reverse if necessary
    # This creates the aligned dct, which has the species and rxn names all the same and in the same order 
    aligned_rxn_ktp_dcts = copy.copy(renamed_rxn_ktp_dcts)
    aligned_rxn_em_dcts = copy.copy(renamed_rxn_em_dcts)
    for mech_idx in range(num_mechs-1):
        for idx2 in range(mech_idx+1, num_mechs):    
            
            # If indicated, reverse the rxn_ktp_dcts and rxm_em_dcts
            if rev_rates:
                aligned_rxn_ktp_dct = reverse_rxn_ktp_dct(renamed_rxn_ktp_dcts[mech_idx], renamed_rxn_ktp_dcts[idx2],
                    renamed_rxn_param_dcts[mech_idx], renamed_rxn_param_dcts[idx2], renamed_thermo_dcts[mech_idx], temps, rev_rates)
                aligned_rxn_em_dct = reverse_rxn_em_dct(renamed_rxn_em_dcts[idx2], renamed_rxn_param_dcts[mech_idx], renamed_rxn_param_dcts[idx2], rev_rates)

            # Otherwise, just rename so the reactants and products are in the same order
            else:
                aligned_rxn_ktp_dct = reverse_rxn_ktp_dct(renamed_rxn_ktp_dcts[mech_idx], renamed_rxn_ktp_dcts[idx2],
                    renamed_rxn_param_dcts[mech_idx], renamed_rxn_param_dcts[idx2], [], temps, rev_rates)
                aligned_rxn_em_dct = reverse_rxn_em_dct(renamed_rxn_em_dcts[idx2], renamed_rxn_param_dcts[mech_idx], renamed_rxn_param_dcts[idx2], rev_rates)
 
            # Store   
            aligned_rxn_ktp_dcts[idx2] = aligned_rxn_ktp_dct
            aligned_rxn_em_dcts[idx2] = aligned_rxn_em_dct
            
#    print('aligned_rxn_ktp_dcts\n', aligned_rxn_ktp_dcts)

    # Sort the aligned dictionaries into a single output
    # Loop over each mechanism
    combined_rxn_ktp_dct = {}
    combined_rxn_em_dct = {}
    for mech_idx, dct in enumerate(aligned_rxn_ktp_dcts):
        #print('mech_idx\n', mech_idx)
        # Loop over each rxn in the mechanism
        for rxn, ktp_dct in dct.items():
            #print('ktp_dct in loop\n', ktp_dct)

            # If the reaction does not yet exist, add it
            if rxn not in combined_rxn_ktp_dct.keys():
                #print('if statement!')
                ktp_dct_list = [None] * mech_idx  # add empty entries to account for previous mechs that are missing
                #print(ktp_dct_list)
                ktp_dct_list.append(ktp_dct)
                combined_rxn_ktp_dct[rxn] = ktp_dct_list
                #print(ktp_dct_list.append(ktp_dct))
                combined_rxn_em_dct[rxn] = aligned_rxn_em_dcts[mech_idx][rxn]  # add the em_param

            # If the reaction already exists, append the new ktp_dct
            else:
                #print('else statement!')
                ktp_dct_list = combined_rxn_ktp_dct[rxn]  # get the current list of ktp_dcts
                # If any of the previous entries were blank, add None entries to fill
                if len(ktp_dct_list) < mech_idx:
                     ktp_dct_list.extend([None] * (mech_idx - len(ktp_dct_list)))
                ktp_dct_list.append(ktp_dct) 
            
    # Clean up the combined dct
    for rxn, ktp_dct_list in combined_rxn_ktp_dct.items():
        # Add None entries so that all are the same length
        if len(ktp_dct_list) < num_mechs:
            ktp_dct_list.extend([None] * (num_mechs - len(ktp_dct_list))) 
            combined_rxn_ktp_dct[rxn] = ktp_dct_list
#    print(combined_rxn_ktp_dct)
            
    return combined_rxn_ktp_dct, combined_rxn_em_dct


def remove_loners(combined_rxn_ktp_dct, combined_rxn_em_dct):
    """ Removes any reactions that don't have rates from all the involved mechs


    """
    filtered_rxn_ktp_dct = {}
    filtered_rxn_em_dct = {}
    for rxn, ktp_dcts in combined_rxn_ktp_dct.items():
        num_mechs = len(ktp_dcts)
        num_rates = len([ktp_dct for ktp_dct in ktp_dcts if ktp_dct is not None])

        # Only add the rxn if all mechanisms have an entry for the rate
        if num_rates == num_mechs: 
            filtered_rxn_ktp_dct[rxn] = ktp_dcts
            filtered_rxn_em_dct[rxn] = combined_rxn_em_dct[rxn]

    #print('filtered_rxn_ktp_dct\n', filtered_rxn_ktp_dct)

    return filtered_rxn_ktp_dct, filtered_rxn_em_dct 
        

def reverse_rxn_ktp_dct(rxn_ktp_dct1, rxn_ktp_dct2, rxn_param_dct1, rxn_param_dct2, thermo_dct1, temps, rev_rates):
    """ This takes two rxn_ktp_dcts whose species have already been renamed
        to be identical and reverses any reactions *in the second dct* that 
        need to be reversed

        :param rxn_ktp_dct1: rxn_ktp_dct for mech1 
        :type rxn_ktp_dct1: dict of dicts {((rcts,prds)):({pressure:(temp_array,rate_array)}, ...)}
        :param rxn_ktp_dct2: rxn_ktp_dct for mech2 
        :type rxn_ktp_dct2: dict of dicts {((rcts,prds)):({pressure:(temp_array,rate_array)}, ...)}
        :param rxn_param_dct1: rxn_param_dct for mech1
        :type rxn_param_dct1: dict{((rcts,prds)):(
    """

    rev_rxn_ktp_dct2 = copy.copy(rxn_ktp_dct2)
    for rxn1, params1 in rxn_param_dct1.items():
        rxn2, flip_rxn, p_dep_same = assess_rxn_match(rxn1, params1, rxn_param_dct2)

        # If the user indicated to reverse rates, check if they need to be 
        if rev_rates:
            try:
                if flip_rxn and p_dep_same:
                    ktp_dct2 = rxn_ktp_dct2[rxn2]
                    rev_ktp_dct2 = reverse_ktp_dct(ktp_dct2, thermo_dct1, rxn2, temps)
                    rev_rxn_ktp_dct2.pop(rxn2)
                    rev_rxn_ktp_dct2[rxn1] = rev_ktp_dct2            
                elif p_dep_same:  # if the reactions should not be reversed, just rename to get reactants and products ordered properly
                    rev_rxn_ktp_dct2[rxn1] = rev_rxn_ktp_dct2[rxn2]            
                    rev_rxn_ktp_dct2.pop(rxn2)
            except KeyError: 
                print('KeyError in compare.reverse_rxn_ktp_dct with the following rxns:\n', rxn1, '\n', rxn2)
        # Otherwise, just rename the rxn so that the reactants and products are in the same order
        else:
#            print('inside compare.align_mechs, else statement\n', rxn2, '\n', rxn1)    
            if p_dep_same and not flip_rxn:  # only do this if the rxns should not be flipped!
                try:
#                    print('inside compare.align_mechs, rename\n', rxn2, '\n', rxn1) 
                    rev_rxn_ktp_dct2[rxn1] = rev_rxn_ktp_dct2[rxn2]            
                    rev_rxn_ktp_dct2.pop(rxn2)
                except KeyError:
                    print('KeyError in compare.reverse_rxn_ktp_dct with the following rxns:\n', rxn1, '\n', rxn2)
                    pass
            
    return rev_rxn_ktp_dct2


def reverse_rxn_em_dct(rxn_em_dct2, rxn_param_dct1, rxn_param_dct2, rev_rates):
    """ Takes a rxn_em_dct and reverses the names of any reactions that need to be reversed

    :param rxn_em_dct2: dict of True/False for whether a reaction does/doesn't have '+M'
    :type rxn_em_dct2: dict {((rcts), (prds)): em}
    :param rxn_param_dct1: rxn_param_dct for mech1
    :type rxn_param_dct1: dict {((rcts), (prds)): [params]}
    :param rxn_param_dct2: rxn_param_dct for mech2
    :type rxn_param_dct2: dict {((rcts), (prds)): [params]}
    :return rev_rxn_em_dct2: rxn_em_dct with reversed rxn names as needed to match rxn_param_dct1
    :rtype: dict
    """
    rev_rxn_em_dct2 = copy.copy(rxn_em_dct2)
    for rxn1, params1 in rxn_param_dct1.items():
        rxn2, flip_rxn, p_dep_same = assess_rxn_match(rxn1, params1, rxn_param_dct2)
        # If the user indicated to reverse rates and it is necessary, reverse the rxn name
        if rev_rates:
            if flip_rxn and p_dep_same:
                em = rxn_em_dct2[rxn2]
                rev_rxn_em_dct2.pop(rxn2)
                rev_rxn_em_dct2[rxn1] = em
        # Otherwise, just rename the rxn so that the reactants and products are in the same order
        else:        
            if rxn2 is not None and p_dep_same and not flip_rxn:
                try:     
                    rev_rxn_em_dct2[rxn1] = rev_rxn_em_dct2[rxn2]            
                    rev_rxn_em_dct2.pop(rxn2)
                except KeyError:
                    print('KeyError in compare.reverse_rxn_em_dct with the following rxns:\n', rxn1, '\n', rxn2)
                    pass    
        

    return rev_rxn_em_dct2


def rename_dcts(dcts, dcts_type, spc_dcts):
    """ Takes a list of dictionaries and renames all the species. The species are renamed in order of the preference
        specified by the order of the list (first dct is unchanged, second is only changed by first, etc.).

        :param dcts: list of dictionaries to be renamed
        :type dcts: list of dcts; [{dct1}, {dct2}, ...]
        :param dcts_type: either 'rxn' or 'spc'; refers to what the key of the dct is
        :type dcts_type: str
        :param spc_dcts: list of species dictionaries corresponding to dcts
        :type spc_dcts: list of dcts; [{dct1}, {dct2}, ...]
        :return renamed_dcts: dcts with species renamed
        :rtype: list of dcts; [{dct1}, {dct2}, ...]
    """
    assert len(dcts) == len(spc_dcts), (
        f'Length of dct_list is {len(dcts)}, while length of spc_dct_lst is {len(spc_dcts)}.'
        )
    num_mechs = len(dcts)

    # Copy both sets of dcts
    renamed_dcts = copy.copy(dcts)
    renamed_spc_dcts = copy.copy(spc_dcts)

    # Loop through each item in the list of dictionaries
    for mech_idx in range(num_mechs-1):
        for idx2 in range(mech_idx+1, num_mechs):

            # Get the rename instructions from the species dcts
            _, rename_instructions = combine_species(renamed_spc_dcts[mech_idx], renamed_spc_dcts[idx2])

            # Rename and store the current spc_dct
            renamed_spc_dct = rename_species(renamed_spc_dcts[idx2],rename_instructions, target_type='spc')
            renamed_spc_dcts[idx2] = renamed_spc_dct

            # Rename and store the current dct
            renamed_dct = rename_species(renamed_dcts[idx2], rename_instructions, dcts_type)
            renamed_dcts[idx2] = renamed_dct

    return renamed_dcts


def combine_species(mech1_spc_dct, mech2_spc_dct):
    """ Combine two species dictionaries
    """

    combined_spc_dct = copy.copy(mech1_spc_dct)
    rename_instructions = {}
    rename_str = '-zz'
    
    # Loop through each species in mech1
    for spc_name1, spc_vals1 in mech1_spc_dct.items():
        i1 = spc_vals1['inchi']
        m1 = spc_vals1['mult']
        c1 = spc_vals1['charge']

        for spc_name2, spc_vals2 in mech2_spc_dct.items():
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
    for spc_name2, spc_vals2 in mech2_spc_dct.items():
        unique = True
        i2 = spc_vals2['inchi']
        m2 = spc_vals2['mult']
        c2 = spc_vals2['charge']

        for spc_name1, spc_vals1 in mech1_spc_dct.items():
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
    #print('rename_instructions\n', rename_instructions)        

    return combined_spc_dct, rename_instructions


def combine_mech_ktps(rxn_ktp_dct1, rxn_ktp_dct2, rename_instructions):
    """ Combine together two rxn_ktp_dcts to create a single dct,
        renaming all species as indicated in the rename_instructions.

        This is largely a useless function.
    """
    combined_rxn_ktp_dct = copy.copy(rxn_ktp_dct1)
    renamed_rxn_ktp_dct2 = rename_species(rxn_ktp_dct2, rename_instructions)

    for rxn_name2, ktp2 in renamed_rxn_ktp_dct2.items():
        rxn_name, flip_rxn, p_dep_same = assess_rxn_match(rxn_name2, rxn_ktp_dct1)    
        
        # If rxn is unique to ktp2, add to ktp1        
        if not rxn_name:
            combined_rxn_ktp_dct[rxn_name2] = ktp2

    return combined_rxn_ktp_dct


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


def rename_species(target_dct, rename_instructions, target_type='rxn'):
    """ Rename the species inside a rxn_ktp_dct, rxn_param_dct, or
        thermo_dct according to the instructions inside the 
        rename_instructions.

        :param target_dct: the dct whose species are to be renamed
        :type target_dct: dct; either a rxn_ktp, rxn_param, or thermo dct
        :param rename_instructions 
    """
    assert target_type in ('rxn', 'thermo', 'spc'), (
        f'The target_type is {target_type}, but should be either "rxn", "thermo", or "spc"'
        )
    renamed_dct = {}

    # If a rxn_ktp or rxn_param dct
    if target_type == 'rxn':
        for rcts, prds in target_dct.keys():
            new_rcts = []
            new_prds = []
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
    
            new_rcts = tuple(new_rcts)
            new_prds = tuple(new_prds)
            renamed_dct[new_rcts, new_prds] = target_dct[rcts, prds]

    # If a thermo or species dct
    else:
        for spc, data in target_dct.items():
            if spc in rename_instructions.keys():
                new_spc_name = rename_instructions[spc]
                renamed_dct[new_spc_name] = data
            else:
                renamed_dct[spc] = data

    return renamed_dct


def assess_rxn_match(rxn_name1, params1, rxn_param_dct2, print_output=False):
    """ assess whether the reaction should be flipped
    """
    # Get all possible orderings of the reactants and products for mech1
    [rcts1, prds1] = rxn_name1
    rcts1_perm = list(itertools.permutations(rcts1, len(rcts1)))
    prds1_perm = list(itertools.permutations(prds1, len(prds1)))
    em_param1 = params1[6]  #'+M', '(+M)', or None

    rxn_name2 = None
    flip_rxn = None
    p_dep_same = None
    for rxn_name, params in rxn_param_dct2.items():
        [rcts2, prds2] = rxn_name
        em_param2 = params[6]
        if rcts2 in rcts1_perm and prds2 in prds1_perm:
            rxn_name2 = rxn_name
            flip_rxn = False

            if em_param2 == em_param1:
                p_dep_same = True
            else:
                p_dep_same = False
            break

        if rcts2 in prds1_perm and prds2 in rcts1_perm:
            rxn_name2 = rxn_name
            flip_rxn = True

            if em_param2 == em_param1:
                p_dep_same = True
            else:
                p_dep_same = False
            break

    if print_output:
        if rxn_name2:
            print('\n- M2 RCTS', '+'.join(mech2_key[0]))
            print('- M2 PRDS', '+'.join(mech2_key[1]))
        else:
            print('\n- NO M2 MATCH')

    return rxn_name2, flip_rxn, p_dep_same


def reverse_ktp_dct(ktp_dct, thermo_dct, rxn, temps):
    """ For a given reaction, use the thermochemistry values of its
        constituent species to calculate the equilibrium constant
        and reverse the rate constants.

        :param ktp_dct: k(T,P) dictionary for a single reaction
        :type ktp_dct: dict[pressure: (temp_array, rates_array)]
        :param thermo_dct: thermochemical values of all species in mechanism
        :type thermo_dct: dict[spc name: [thermo vals]]
        :param rxn: reactant-product pair for the reaction
        :type rxn: tuple(tuple(str), tuple(str))
        :param temps: list of temperatures 
        :type temps: list(float)
        :return: rev_ktp_dct: reversed rates of the reaction
        :rtype: dict[pressure: (temp_array, reversed_rates_array)]
    """

    [rcts, prds] = rxn
    k_equils = _calculate_equilibrium_constant(
        thermo_dct, rcts, prds, temps)

    rev_ktp_dct = {}
    for pressure, (_, rate_ks) in ktp_dct.items():

        # Calculate density to handle units, if needed
        if len(rcts) > 1 and len(prds) == 1:
            densities = ratefit.calc.p_to_m(1.0, temps)
            rate_ks *= densities
            # print('flip1')
        elif len(rcts) == 1 and len(prds) > 1:
            densities = ratefit.calc.p_to_m(1.0, temps)
            rate_ks /= densities
            # print('flip2')

        # Calculate the reverse rates with K_equil
        rev_rates = []
        for forw_k, k_equil in zip(rate_ks, k_equils):
            rev_rates.append(forw_k / k_equil)

        # Add reversed rates to dict
        rev_ktp_dct[pressure] = (temps, rev_rates)

    return rev_ktp_dct


def _calculate_equilibrium_constant(thermo_dct, rct_idxs, prd_idxs, temps):
    """ Calculate the equilibrium constant for a given reaction at
        a set of temperatures using constituent species' thermochemistry.

        :param thermo_dct: thermochemical values of all species in mechanism
        :type thermo_dct: dict[spc name: [thermo vals]]
        :param rct_idxs: name(s) of the reactants
        :type rct_idxs: tuple(str)
        :param prd_idxs: name(s) of the products
        :type prd_idxs: tuple(str)
        :param temps: Temperatures to calculate thermochemistry (K)
        :type temps: list(float)
        :return: equilibrium constants
        :rtype: numpy.ndarray
    """

    k_equils = []
    for temp_idx, temp in enumerate(temps):
        rct_gibbs = 0.0
        for rct in rct_idxs:
            rct_gibbs += _grab_gibbs(thermo_dct[rct], temp_idx)
            # print('idv rct g', _grab_gibbs(thermo_dct[rct], temp_idx))
        prd_gibbs = 0.0
        for prd in prd_idxs:
            prd_gibbs += _grab_gibbs(thermo_dct[prd], temp_idx)
            # print('idv prd g', _grab_gibbs(thermo_dct[prd], temp_idx))

        rxn_gibbs = prd_gibbs - rct_gibbs

        k_equils.append(
            np.exp(-rxn_gibbs / (phycon.RC * temp)))
        ktest = np.exp(-rxn_gibbs / (phycon.RC * temp))

#        print('rct gibbs', rct_gibbs)
#        print('prd gibbs', prd_gibbs)
#        print('rxn gibbs', rxn_gibbs)
#        print('kequil', ktest)

    return k_equils


def _grab_gibbs(thermo_vals, temp_idx):
    """ calculate the Gibbs Free energy value
    """
    gibbs = thermo_vals[3][temp_idx]
    return gibbs


def get_rxn_em_dct(rxn_param_dct):
    
    rxn_em_dct = {}
    for rxn, params in rxn_param_dct.items():
        em_param = False
        if params[6] == '+M':
            em_param = True
        
        rxn_em_dct[rxn] = em_param
 
    return rxn_em_dct
   
def format_rxn_name(rxn_key, em):
    """ Receives a rxn_key and the corresponding param_vals
        from a rxn_param_dct and writes it to a string that
        the above functions can handle. Adds +M or if
        applicable.
    """
    rcts = rxn_key[0]
    prds = rxn_key[1]
    for idx, rct in enumerate(rcts):
        if idx == 0:
            rct_str = rct
        else:
            rct_str += ' + ' + rct
    for idx, prd in enumerate(prds):
        if idx == 0:
            prd_str = prd
        else:
            prd_str += ' + ' + prd

    # Add '+ M' if it is applicable
    if em:
        rct_str += ' + M'
        prd_str += ' + M'

    rxn_name = rct_str + ' = ' + prd_str

    return rxn_name
 
