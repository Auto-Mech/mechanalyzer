"""
  Take data dictionaries from mechanisms and combine them under a common index
"""

import itertools
import numpy as np
import ratefit
from ioformat import remove_whitespace
from ioformat import phycon
from chemkin_io.parser import mechanism as mech_parser
from chemkin_io.parser.mechanism import reaction_units
from mechanalyzer.calculator import thermo
from mechanalyzer.calculator import rates


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


def mechanism_rates(mech1_ktp_dct, mech2_ktp_dct, temps,
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

    total_ktp_dct = {}

    print('\nFinding Reaction Matches to Combine Dictionaries')
    # Build full rates dictionary with common index
    # First loop through mech1: add common and mech1-unique species
    for mech1_name, mech1_ktp in mech1_ktp_dct.items():

        print('---------')
        print('- M1 RCTS', '+'.join(mech1_name[0]))
        print('- M1 PRDS', '+'.join(mech1_name[1]))
        # Check what (if/any) combination of mech2 matches with mech1
        mech2_name_match, reverse_rates = _assess_reaction_match(
            mech1_name, mech2_ktp_dct)

        # Calculate reaction rates, reverse if needed
        if mech2_name_match:
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
def _assess_reaction_match(mech1_names, mech2_dct):
    """ assess whether the reaction should be flipped
    """

    # Get all possible orderings of the reactants and products for mech1
    [mech1_rcts, mech1_prds] = mech1_names
    mech1_rct_perm = list(itertools.permutations(mech1_rcts, len(mech1_rcts)))
    mech1_prd_perm = list(itertools.permutations(mech1_prds, len(mech1_prds)))
    # print('\n\nTEST')
    # print('m1_rct_perm', mech1_rct_perm)
    # print('m1_prd_perm', mech1_prd_perm)

    mech2_key = ()
    flip_rxn = None
    for rxn in mech2_dct:
        [mech2_rcts, mech2_prds] = rxn
        # print('m2', mech2_rcts, mech2_prds)
        if mech2_rcts in mech1_rct_perm and mech2_prds in mech1_prd_perm:
            flip_rxn = False
            mech2_key = rxn
            break
        if mech2_rcts in mech1_prd_perm and mech2_prds in mech1_rct_perm:
            mech2_key = rxn
            flip_rxn = True
            break

    if mech2_key:
        print('\n- M2 RCTS', '+'.join(mech2_key[0]))
        print('- M2 PRDS', '+'.join(mech2_key[1]))
    else:
        print('\n- NO M2 MATCH')

    return mech2_key, flip_rxn


def _reverse_reaction_rates(mech_dct, thermo_dct, rxn, temps):
    """ For a given reaction, use the thermochemistry values of its
        constituent species to calculate the equilibrium constant
        and reverse the rate constants.

        :param mech_dct:
        :type mech_dct: dict[pressure: rates]
        :param thermo_dct: thermochemical values of all species in mechanism
        :type thermo_dct: dict[spc name: [thermo vals]]
        :param rxn: reactant-product pair for the reaction
        :type rxn: tuple(tuple(str), tuple(str))
        :return: rev_ktp_dct: reversed rates of the reaction
        :rtype: dict[pressure: reversed rates]
    """

    [rct_idxs, prd_idxs] = rxn
    k_equils = _calculate_equilibrium_constant(
        thermo_dct, rct_idxs, prd_idxs, temps)

    ktp_dct = mech_dct[rxn]
    rev_ktp_dct = {}
    for pressure, rate_ks in ktp_dct.items():

        # Calculate density to handle units, if needed
        if len(rct_idxs) > 1 and len(prd_idxs) == 1:
            densities = ratefit.calc.p_to_m(1.0, temps)
            rate_ks *= densities
            # print('flip1')
        elif len(rct_idxs) == 1 and len(prd_idxs) > 1:
            densities = ratefit.calc.p_to_m(1.0, temps)
            rate_ks /= densities
            # print('flip2')

        # Calculate the reverse rates with K_equil
        rev_rates = []
        for forw_k, k_equil in zip(rate_ks, k_equils):
            rev_rates.append(forw_k / k_equil)

        # Add reversed rates to dict
        rev_ktp_dct[pressure] = rev_rates

    return rev_ktp_dct


def _calculate_equilibrium_constant(thermo_dct, rct_idxs, prd_idxs, temps):
    """ Calculate the equilibrium constant for a given reaction at
        a set of temperatures using constiutent species' thermochemistry.

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
