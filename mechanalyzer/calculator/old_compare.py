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
from mechanalyzer.calculator import combine as calc_combine
from chemkin_io.parser.mechanism import reaction_units
import copy
import os

# Read the job path, which should have been input as the current directory
JOB_PATH = sys.argv[1]


def plot_aligned_mechs(aligned_rxn_ktp_dcts):
    """ Produces plots of aligned rxn_ktp_dcts
        
        :param aligned_rxn_ktp_dcts: a list of aligned rxn_ktp_dcts,
            all with aligned species names and directions
        :type aligned_rxn_ktp_dcts: list({dict{dict}}) [{((rcts),(prds)):{pres:(temps_array,rates_array)}, ...}, ...]
    """
    

def align_mechs(mech_filenames, spc_csv_filenames, thermo_filenames, temps, pressures):
    """ Calculates kTP dictionaries for any number of mechanisms and then
        flips the reactions within the mechanisms to all be in the same direction


    """
    assert len(mech_filenames) == len(spc_csv_filenames) == len(thermo_filenames), (
        f"""The lengths of the mechanism, spc_csv, and thermo filename inputs should 
            all be the same, but are instead {len(mech_filenames)}, {len(spc_csv_filenames)}, 
            and {len(thermo_filenames)}, respectively."""
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
        rxn_param_dct = parser_rxn.param_dct(rxn_block_str)
        rxn_ktp_dct = calc_rates.eval_rxn_param_dct(rxn_param_dct, pressures, temps)
        rxn_ktp_dcts.append(rxn_ktp_dct)
        rxn_param_dcts.append(rxn_param_dct)

        # Get the spc_dct
        spc_csv_str = parser.ptt.read_inp_str(JOB_PATH, spc_csv_filenames[idx], remove_comments=False)
        spc_dct = parser_spc.build_spc_dct(spc_csv_str, 'csv')
        spc_dcts.append(spc_dct)

        # Get the thermo_dct
        thermo_str = parser.ptt.read_inp_str(JOB_PATH,mech_filenames[idx],remove_comments=False)
        thermo_block_str = parser_mech.thermo_block(thermo_str)
        thermo_dct = calc_thermo.mechanism(thermo_block_str, temps)
        thermo_dcts.append(thermo_dct)

    # Get the renamed dictionaries
    renamed_rxn_ktp_dcts = []
    renamed_rxn_param_dcts = []
    renamed_spc_dcts = []
    renamed_thermo_dcts = []
    for mech_idx in range(num_mechs-1):
        # If on first mech, don't rename, just copy
        if mech_idx == 0:
            renamed_rxn_ktp_dcts.append(rxn_ktp_dcts[mech_idx])
            renamed_rxn_param_dcts.append(rxn_param_dcts[mech_idx])
            renamed_spc_dcts.append(spc_dcts[mech_idx])            
            renamed_thermo_dcts.append(thermo_dcts[mech_idx])
        # Loop through the mechs remaining after the current one
        for idx2 in range(mech_idx+1, num_mechs):
            # If on the first time through, do things a bit differently
            if mech_idx == 0:
                # Get the instructions on renaming species, aka the rename_spc_dct
                _, rename_spc_dct = calc_combine.combine_species(spc_dcts[mech_idx], spc_dcts[idx2])

                # Rename the species in the various dictionaries
                renamed_rxn_ktp_dct = calc_combine.rename_species(rxn_ktp_dcts[idx2], rename_spc_dct, target_type='rxn')  
                renamed_rxn_param_dct = calc_combine.rename_species(rxn_param_dcts[idx2], rename_spc_dct, target_type='rxn')  
                renamed_spc_dct = calc_combine.rename_species(spc_dcts[idx2], rename_spc_dct, target_type='spc')  
                renamed_thermo_dct = calc_combine.rename_species(thermo_dcts[idx2], rename_spc_dct, target_type='thermo')
                
                # Store the results
                renamed_rxn_ktp_dcts.append(renamed_rxn_ktp_dct)
                renamed_rxn_param_dcts.append(renamed_rxn_param_dct)
                renamed_spc_dcts.append(renamed_spc_dct)
                renamed_thermo_dcts.append(renamed_thermo_dct)

            else:
                # Get the instructions on renaming species, aka the rename_spc_dct
                _, rename_spc_dct = calc_combine.combine_species(renamed_spc_dcts[mech_idx], renamed_spc_dcts[idx2])

                # Rename the species in the various dictionaries
                renamed_rxn_ktp_dct = calc_combine.rename_species(renamed_rxn_ktp_dcts[idx2], rename_spc_dct, target_type='rxn')  
                renamed_rxn_param_dct = calc_combine.rename_species(renamed_rxn_param_dcts[idx2], rename_spc_dct, target_type='rxn')  
                renamed_spc_dct = calc_combine.rename_species(renamed_spc_dcts[idx2], rename_spc_dct, target_type='spc')  
                renamed_thermo_dct = calc_combine.rename_species(renamed_thermo_dcts[idx2], rename_spc_dct, target_type='thermo')

                # Store the results                
                renamed_rxn_ktp_dcts[idx2] = renamed_rxn_ktp_dct
                renamed_rxn_param_dcts[idx2] = renamed_rxn_param_dct
                renamed_spc_dcts[idx2] = renamed_spc_dct
                renamed_thermo_dcts[idx2] = renamed_thermo_dct
    
    # Get the em_param_dcts                    
    renamed_em_param_dcts = []
    for dct in renamed_rxn_param_dcts:
        renamed_em_param_dct = get_em_param_dct(dct)
        renamed_em_param_dcts.append(renamed_em_param_dct)

    # Loop over each reaction and reverse if necessary
    aligned_rxn_ktp_dcts = []
    aligned_em_param_dcts = []
    for mech_idx in range(num_mechs-1):
        if mech_idx == 0:
            aligned_rxn_ktp_dcts.append(renamed_rxn_ktp_dcts[mech_idx])
            aligned_em_param_dcts.append(renamed_em_param_dcts[mech_idx])
        for idx2 in range(mech_idx+1, num_mechs):    
            if mech_idx == 0:
                aligned_rxn_ktp_dct = reverse_ktp_dct(renamed_rxn_ktp_dcts[mech_idx], renamed_rxn_ktp_dcts[idx2],
                    renamed_rxn_param_dcts[mech_idx], renamed_rxn_param_dcts[idx2], renamed_thermo_dcts[mech_idx], temps)
                    
                aligned_rxn_ktp_dcts.append(aligned_rxn_ktp_dct)

            else:
                aligned_rxn_ktp_dct = reverse_ktp_dct(renamed_rxn_ktp_dcts[mech_idx], renamed_rxn_ktp_dcts[idx2],
                    renamed_rxn_param_dcts[mech_idx], renamed_rxn_param_dcts[idx2], renamed_thermo_dcts[mech_idx], temps)
                # Note: the preceding line does not necessarily use the thermo of the first mechanism; it only is sure
                # to use the thermo of the mechanism currently being used as the reference

                aligned_rxn_ktp_dcts[idx2] = aligned_rxn_ktp_dct

    # Sort the aligned dictionaries into a single output
    # Loop over each mechanism
    combined_rxn_ktp_dct = {}
    for mech_idx, dct in enumerate(aligned_rxn_ktp_dcts):

        # Loop over each rxn in the mechanism
        for rxn, ktp_dct in dct.items():

            # If the reaction does not yet exist, add it
            if rxn not in combined_rxn_ktp_dct.keys():
                ktp_dct_list = [None] * mech_idx
                combined_rxn_ktp_dct[rxn] = ktp_dct_list.append(ktp_dct)

                # Also, add the em_param
                

            # If the reaction already exists, append the new ktp_dct
            else:
                ktp_dct_list = combined_rxn_ktp_dct[rxn]  # get the current list of ktp_dcts

                # If any of the previous entries were blank, add None entries to fill
                if len(ktp_dct_list) < mech_idx:
                     ktp_dct_list.extend([None] * (mech_idx - len(ktp_dct_list)))
                                
                ktp_dct_list.append(ktp_dct) 
            
    # Clean up the combined dct in two ways: 
    for rxn, ktp_dct_list in combined_rxn_ktp_dct.items():
        # 1: Add None entries so that all are the same length
        if len(ktp_dct_list) < num_mechs:
            ktp_dct_list.extend([None] * (num_mechs - len(ktp_dct_list))) 
            combined_rxn_ktp_dct[rxn] = ktp_dct_list
        # 2: Create the em_param_dct, which denotes whether a reaction has 
        
            
    return combined_rxn_ktp_dct


def reverse_rxn_ktp_dct(rxn_ktp_dct1, rxn_ktp_dct2, rxn_param_dct1, rxn_param_dct2, thermo_dct1, temps)
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
    rev_rxn_em_dct2 = {}
    for rxn1, params1 in rxn_param_dct1.items():
        rxn2, flip_rxn, p_dep_same = _assess_rxn_match(rxn1, params1, rxn_param_dct2)
        if flip_rxn and p_dep_same:
            ktp_dct2 = rxn_ktp_dct2[rxn2]
            rev_ktp_dct2 = _reverse_ktp_dct(ktp_dct2, thermo_dct1, rxn_name2, temps)
            rev_rxn_ktp_dct2.pop(rxn2)
            rev_rxn_ktp_dct2[rxn1] = rev_ktp_dct2            
            # Update the 
            if params1[6] == '+M':
                rev_rxn_em_dct2[rxn1] = True
            else: 
                rev_rxn_em_dct2[rxn1] = False
                
        else:
            if params1[6] == '+M':
                rev_rxn_em_dct2[rxn1] = True
            else: 
                rev_rxn_em_dct2[rxn1] = False
            
    return rev_rxn_ktp_dct2


def reverse_ktp_values(rxn_name2, rxn_ktp_dct2, thermo_dct1):
    """ Takes the k(T,P) values for a single reaction in 
        rxn_ktp_dct2 and reverses them according to the thermo values 
        in thermo
    """
    ktps = rxn_ktp_dct2[rxn_name2]
    
     


def combine_species(mech1_spc_dct, mech2_spc_dct):
    """ Combine two species dictionaries
    """

    combined_spc_dct = copy.copy(mech1_spc_dct)
    rename_spc_dct = {}
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
                    rename_spc_dct[spc_name2] = spc_name1

            # If species are different but have same name
            elif spc_name1 == spc_name2:  
                print('loop 1: ' + spc_name2)
                rename_spc_dct[spc_name2] = spc_name2 + rename_str                 

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
                if spc_name2 in rename_spc_dct.keys():
                    unique = False
                break
                
        if unique:
            if spc_name2 in rename_spc_dct.keys():
                combined_spc_dct[spc_name2 + rename_str] = spc_vals2
            else:
                combined_spc_dct[spc_name2] = spc_vals2
            
    return combined_spc_dct, rename_spc_dct


def combine_mech_ktps(rxn_ktp_dct1, rxn_ktp_dct2, rename_spc_dct):
    """ Combine together two rxn_ktp_dcts to create a single dct,
        renaming all species as indicated in the rename_spc_dct.

        This is largely a useless function.
    """
    combined_rxn_ktp_dct = copy.copy(rxn_ktp_dct1)
    renamed_rxn_ktp_dct2 = rename_species(rxn_ktp_dct2, rename_spc_dct)

    for rxn_name2, ktp2 in renamed_rxn_ktp_dct2.items():
        rxn_name, flip_rxn, p_dep_same = _assess_reaction_match(rxn_name2, rxn_ktp_dct1)    
        
        # If rxn is unique to ktp2, add to ktp1        
        if not rxn_name:
            combined_rxn_ktp_dct[rxn_name2] = ktp2

    return combined_rxn_ktp_dct


def combine_mech_params(rxn_param_dct1, rxn_param_dct2, rename_spc_dct):
    """ Combine together two rxn_param_dcts to create a single dct,
        renaming all species as indicated in the rename_spc_dct.
    """
    combined_rxn_param_dct = copy.copy(rxn_param_dct1)
    renamed_rxn_param_dct2 = rename_species(rxn_param_dct2, rename_spc_dct)

    for rxn_name2, param2 in renamed_rxn_param_dct2.items():
        rxn_name, flip_rxn, p_dep_same = _assess_reaction_match(rxn_name2, param2, rxn_param_dct1)    
        
        # If rxn is unique to param2, add to param 1. Note: this does not
        # consider whether the pressure dependence is the same. It just keeps
        # the rate from mechanism 1. 
        if not rxn_name:
            combined_rxn_param_dct[rxn_name2] = param2

    return combined_rxn_param_dct


def rename_species(rxn_dct, rename_spc_dct):
    """ Rename the species inside a rxn_ktp_dct OR
        a rxn_param_dct according to the instructions
        inside the rename_spc_dct
    """
    renamed_rxn_dct = {}
    for rcts, prds in rxn_dct.keys():
        new_rcts = []
        new_prds = []
        for spc in rcts:
            if spc in rename_spc_dct.keys():
                new_rcts.append(rename_spc_dct[spc])
            else:
                new_rcts.append(spc)
        for spc in prds:
            if spc in rename_spc_dct.keys():
                new_prds.append(rename_spc_dct[spc])
            else:
                new_prds.append(spc)
    
        new_rcts = tuple(new_rcts)
        new_prds = tuple(new_prds)
        renamed_rxn_dct[new_rcts, new_prds] = rxn_dct[rcts, prds]

    return renamed_rxn_dct


def mechanism_thermo(mech1_thermo_dct, mech2_thermo_dct):
    """ Combine together the thermo dictionaries for two mechanisms.

        :param mech1_thermo_dct: thermo dct for mechanism 1
        :type mech1_thermo_dct: dict[spc: [[H(t)], [Cp(T)], [S(T)], [G(T)]]
        :param mech2_thermo_dct: thermo dct for mechanism 2
        :type mech2_thermo_dct: dict[spc: [[H(t)], [Cp(T)], [S(T)], [G(T)]]
        :return total_thermo_dct: dict with thermo from both mechanisms
        :rtype: dict[mech: mech_thermo_dct]
    """

    total_thermo_dct = {}

    # Build full thermo dictionary with common index
    # First loop through mech1: add common and mech1-unique species
    for mech1_name, mech1_vals in mech1_thermo_dct.items():
        mech2_vals = mech2_thermo_dct.get(mech1_name, [None, None, None, None])
        total_thermo_dct[mech1_name] = {
            'mech1': mech1_vals,
            'mech2': mech2_vals
        }

    # Now add the entries where mech2 exists, but mech1 does not
    uni_mech2_names = [name
                       for name in mech2_thermo_dct.keys()
                       if name not in mech1_thermo_dct]

    for mech2_name in uni_mech2_names:
        total_thermo_dct[mech2_name] = {
            'mech1': [None, None, None, None],
            'mech2': mech2_thermo_dct[mech2_name]
        }

    return total_thermo_dct


def mechanism_rates(rxn_ktp_dct1, rxn_ktp_dct2,
                    mech2_thermo_dct=None, ignore_reverse=False):
    """ Combine together the rate dictionaries for two mechanisms.
        Currently, there are two assumptions to combining the dicts.
          (1) k(T,P) values for both mechanisms have same temperature range.
          (2) Reaction directions in combined dict match those from Mech1.

        :param mech1_ktp_dct: rate constant dict for mechanism 1
        :type mech1_ktp_dct: dict[reaction: dict[pressure: k(T,P)s]]
        :param mech2_ktp_dct: rate constant dict for mechanism 2
        :type mech2_ktp_dct: dict[reaction: dict[pressure: k(T,P)s]]
        :param mech2_thermo_dct: thermo dict for mechanism 2; to reverse rxns
        :type mech2_thermo_dct: dict[spc: [[H(t)], [Cp(T)], [S(T)], [G(T)]]
        :param temps: Temperatures the k(T,P) values were calculated for (K)
        :type temps: list(float)
        :return total_ktp_dct: dict with thermo from both mechanisms
        :rtype: dict[mech: mech_ktp_dct]
    """

    total_rxn_ktp_dct = copy.copy(rxn_ktp_dct1)

    for rxn_name2, ktp2 in rxn_ktp_dct2.items():
        rxn_name, flip_rxn = _assess_reaction_match(rxn_name2, rxn_ktp_dct1)
        
        # If the rxn doesn't exist in ktp1, it's unique to ktp2, so add
        # it to the total ktp_dct
        if not rxn_name:
            total_rxn_ktp_dct[rxn_name2] = ktp2


        print('---------')
        print('- M1 RCTS', '+'.join(mech1_name[0]))
        print('- M1 PRDS', '+'.join(mech1_name[1]))
        # Check what (if/any) combination of mech2 matches with mech1

        rxn_name2, reverse_rates = _assess_reaction_match(
            rxn_name2, rxn_ktp_dct1)

        # Calculate reaction rates, reverse if needed
        if rxn_name2:
            if not reverse_rates:
                print('\n - Match in forward direction')
                mech2_ktp = mech2_ktp_dct[mech2_name_match]
                # mech2_ktp = mech2_ktp_dct[mech1_name]
            else:
                if not ignore_reverse:
                    print('\n - Match in reverse direction, rev k(T) w/ therm')
                    assert mech2_thermo_dct is not None
                    mech2_ktp = _reverse_reaction_rates(
                        mech2_ktp_dct, mech2_thermo_dct,
                        mech2_name_match, temps)
                else:
                    continue
        else:
            mech2_ktp = {}

        # print('\npost flip mech2 ktp')
        # print(mech2_ktp)

        # Add data_entry to overal thermo dictionary
        total_ktp_dct[mech1_name] = {
            'mech1': mech1_ktp,
            'mech2': mech2_ktp
        }

        # Now add the entries where mech2 exists, but mech1 does not
        # add the code to do this

    return total_ktp_dct


# Thermo functions
def build_thermo_name_dcts(mech1_str, mech2_str, temps):
    """ Builds the thermo dictionaries indexed by names.

        :param mech1_str: string of mechanism 1 input file
        :type mech1_str: str
        :param mech2_str: string of mechanism 2 input file
        :type mech2_str: str
        :param temps: Temperatures to calculate thermochemistry (K)
        :type temps: list(float)
        :return: mech1_thermo_dct
        :rtype: dict[name: [thermo]]
        :return: mech2_thermo_dct
        :rtype: dict[name: [thermo]]
    """

    mech1_thermo_block = mech_parser.thermo_block(mech1_str)
    if mech1_thermo_block is not None:
        mech1_thermo_block = remove_whitespace(
            mech1_thermo_block)
        mech1_thermo_dct = thermo.mechanism(
            mech1_thermo_block, temps)
    else:
        mech1_thermo_dct = None

    mech2_thermo_block = mech_parser.thermo_block(mech2_str)
    if mech2_thermo_block is not None:
        mech2_thermo_block = remove_whitespace(
            mech2_thermo_block)
        mech2_thermo_dct = thermo.mechanism(
            mech2_thermo_block, temps)
    else:
        mech2_thermo_dct = None

    return mech1_thermo_dct, mech2_thermo_dct


def build_thermo_inchi_dcts(mech1_str, mech2_str,
                            mech1_csv_str, mech2_csv_str,
                            temps):
    """ Builds the thermo dictionaries indexed by InChI strings.

        :param mech1_str: string of mechanism 1 input file
        :type mech1_str: str
        :param mech2_str: string of mechanism 2 input file
        :type mech2_str: str
        :param mech1_csv_str: species.csv file string for mechanism 1
        :type mech1_csv_str: str
        :param mech2_csv_str: species.csv file string for mechanism 2
        :type mech2_csv_str: str
        :param temps: Temperatures to calculate thermochemistry (K)
        :type temps: list(float)
        :return: mech1_thermo_dct
        :rtype: dict[name: [thermo]]
        :return: mech2_thermo_dct
        :rtype: dict[name: [thermo]]
    """

    # Get dicts: dict[name] = thm_dstr
    mech1_thermo_dct, mech2_thermo_dct = build_thermo_name_dcts(
        mech1_str, mech2_str, temps)

    # Build the inchi dicts where:
    # Get dicts: dict[name] = inchi
    # Convert name dict to get: dict[inchi] = name
    if mech1_thermo_dct is not None:
        mech1_name_inchi_dct = mech_parser.spc_name_dct(
            mech1_csv_str, 'inchi')
        mech1_thermo_ich_dct = {}
        for name, data in mech1_thermo_dct.items():
            ich = mech1_name_inchi_dct[name]
            mech1_thermo_ich_dct[ich] = data
    else:
        mech1_thermo_ich_dct = None

    if mech2_thermo_dct is not None:
        mech2_name_inchi_dct = mech_parser.spc_name_dct(
            mech2_csv_str, 'inchi')
        mech2_thermo_ich_dct = {}
        for name, data in mech2_thermo_dct.items():
            ich = mech2_name_inchi_dct[name]
            mech2_thermo_ich_dct[ich] = data
    else:
        mech2_thermo_ich_dct = None

    return mech1_thermo_ich_dct, mech2_thermo_ich_dct


def spc_name_from_inchi(mech1_csv_str, mech2_csv_str, ich):
    """ uses dict[inchi]=name dicts to get
        the mechanism name for a given InChI string
    """

    mech1_inchi_dct = mech_parser.spc_inchi_dct(mech1_csv_str)
    mech2_inchi_dct = mech_parser.spc_inchi_dct(mech2_csv_str)
    mech_name = mech1_inchi_dct.get(ich)
    if mech_name is not None:
        mech_name = mech2_inchi_dct.get(ich)

    return mech_name


# Rate functions
def _assess_reaction_match(rxn_name1, params1, rxn_param_dct2, print_output=False):
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

        if mech2_rcts in prds1_perm and mech2_prds in rcts1_perm:
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


def _reverse_ktp_dct(ktp_dct, thermo_dct, rxn, temps):
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

        print('rct gibbs', rct_gibbs)
        print('prd gibbs', prd_gibbs)
        print('rxn gibbs', rxn_gibbs)
        print('kequil', ktest)

    return k_equils


def _grab_gibbs(thermo_vals, temp_idx):
    """ calculate the Gibbs Free energy value
    """
    gibbs = thermo_vals[3][temp_idx]
    return gibbs


# Functions to build dictionaries
def build_reaction_name_dcts(mech1_str, mech2_str, t_ref, temps, pressures,
                             ignore_reverse=True, remove_bad_fits=False):
    """ Parses the strings of two mechanism files and calculates
        rate constants [k(T,P)]s at an input set of temperatures and pressures.

        :param mech1_str: string of mechanism 1 input file
        :type mech1_str: str
        :param mech2_str: string of mechanism 2 input file
        :type mech2_str: str
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: List of Temperatures (K)
        :type temps: numpy.ndarray
        :return mech1_ktp_dct: rate constants for mechanism 1
        :rtype: dict[pressure: rates]
        :return mech2_ktp_dct: rate constants for mechanism 2
        :rtype: dict[pressure: rates]
    """

    mech1_reaction_block = remove_whitespace(
        mech_parser.reaction_block(mech1_str))
    mech1_units = reaction_units(mech1_str)
    mech1_ktp_dct = rates.mechanism(
        mech1_reaction_block, mech1_units, t_ref, temps, pressures,
        ignore_reverse=ignore_reverse, remove_bad_fits=remove_bad_fits)

    if mech2_str:
        mech2_reaction_block = remove_whitespace(
            mech_parser.reaction_block(mech2_str))
        mech2_units = reaction_units(mech2_str)
        mech2_ktp_dct = rates.mechanism(
            mech2_reaction_block, mech2_units, t_ref, temps, pressures,
            ignore_reverse=ignore_reverse, remove_bad_fits=remove_bad_fits)
    else:
        mech2_ktp_dct = {}

    return mech1_ktp_dct, mech2_ktp_dct


def build_reaction_inchi_dcts(mech1_str, mech2_str,
                              mech1_csv_str, mech2_csv_str,
                              t_ref, temps, pressures,
                              ignore_reverse=True,
                              remove_bad_fits=False):
    """ builds new reaction dictionaries indexed by inchis
    """
    # Get dicts: dict[name] = rxn_dstr
    mech1_reaction_dct, mech2_reaction_dct = build_reaction_name_dcts(
        mech1_str, mech2_str, t_ref, temps, pressures,
        ignore_reverse=ignore_reverse)

    # Get dicts: dict[name] = inchi
    mech1_name_inchi_dct = mech_parser.spc_name_dct(
        mech1_csv_str, 'inchi')
    mech2_name_inchi_dct = mech_parser.spc_name_dct(
        mech2_csv_str, 'inchi')

    # Convert name dict to get: dict[inchi] = rxn_data
    mech1_reaction_ich_dct = {}
    for names, data in mech1_reaction_dct.items():
        [rct_names, prd_names] = names
        rct_ichs, prd_ichs = (), ()
        for rcts in rct_names:
            rct_ichs += ((mech1_name_inchi_dct[rcts]),)
        for prds in prd_names:
            prd_ichs += ((mech1_name_inchi_dct[prds]),)
        mech1_reaction_ich_dct[(rct_ichs, prd_ichs)] = data

    mech2_reaction_ich_dct = {}
    for names, data in mech2_reaction_dct.items():
        [rct_names, prd_names] = names
        rct_ichs, prd_ichs = (), ()
        for rcts in rct_names:
            rct_ichs += ((mech2_name_inchi_dct[rcts]),)
        for prds in prd_names:
            prd_ichs += ((mech2_name_inchi_dct[prds]),)
        mech2_reaction_ich_dct[(rct_ichs, prd_ichs)] = data

    return mech1_reaction_ich_dct, mech2_reaction_ich_dct


def conv_ich_to_name_ktp_dct(ktp_ich_dct, csv_str):
    """ convert ktp dct from using ichs to using names
    """
    # Get dicts: dict[name] = inchi
    mech1_inchi_name_dct = mech_parser.spc_inchi_dct(csv_str)

    # Convert name dict to get: dict[inchi] = rxn_data
    ktp_name_dct = {}
    for ichs, params in ktp_ich_dct.items():
        [rct_ichs, prd_ichs] = ichs
        rct_names, prd_names = (), ()
        for rcts in rct_ichs:
            rct_names += ((mech1_inchi_name_dct[rcts]),)
        for prds in prd_ichs:
            prd_names += ((mech1_inchi_name_dct[prds]),)
        ktp_name_dct[(rct_names, prd_names)] = params

    return ktp_name_dct


def get_em_param_dct(rxn_param_dct):
    
    em_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        em_param = False
        if params[6] == '+M'
            em_param = True
        
        em_param_dct[rxn] = em_param
 
    return em_param_dct
    

# def spc_name_from_inchi(mech1_csv_str, mech2_csv_str, ich_pair):
#     """ uses dict[inchi]=name dicts to get
#         the mechanism name for a given InChI string
#     """
#    mech1_inchi_dct = chemkin_io.parser.mechanism.spc_inchi_dct(mech1_csv_str)
#    mech2_inchi_dct = chemkin_io.parser.mechanism.spc_inchi_dct(mech2_csv_str)
#
#     if ich in mech1_inchi_dct:
#        mech_name = mech1_inchi_dct[ich]
#     else:
#         mech_name = mech2_inchi_dct[ich]
#
#     return mech_name
