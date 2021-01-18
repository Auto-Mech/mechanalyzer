"""
Create plots which compare thermochemical and kinetic values for two mechanisms
derived from the NASA polynomials and rate constants, respectively. Currently
assymed that the thermochemical and kinetic parameters will be read from
CHEMKIN-formatted mechanism files.
"""

import os
import numpy as np
from ioformat import remove_whitespace
import chemkin_io
from _format import read_file
from mechanalyzer import calculator
from mechanalyzer import plotter


def sm_rates(mech_str, t_ref, temps, pressures, dir_prefix='rate_plots',
             ignore_reverse=False, remove_bad_fits=False):
    """ Calculate rates for a single mechanism
    """
    # Parse mechanism file and cacluate the rates
    mech1_reaction_block = remove_whitespace(
        mech_parser.reaction_block(mech_str))
    mech1_units = chemkin_io.parser.mechanism.reaction_units(mech_str)
    mech1_ktp_dct = calculator.rates.mechanism(
        mech1_reaction_block, mech1_units, t_ref, temps, pressures,
        ignore_reverse=ignore_reverse, remove_bad_fits=remove_bad_fits)

    # Build the plots
    plotter.sm_rates.build(mech1_ktp_dct, temps, dir_prefix=dir_prefix)


def rates(mech1_str, mech2_str, t_ref, temps, pressures,
          mech1_csv_str=None, mech2_csv_str=None,
          dct_idx='name', dir_prefix='rate_plots',
          mech_labels=None, ignore_reverse=False):
    """ Read the reaction sections of two CHEMKIN files and
        produce plots of the rate constants together
    """

    # Build dictionaries containing the thermo data strings for each species
    if dct_idx == 'name':
        # _, mech2_thermo_dct = calculator.combine.build_thermo_name_dcts(
        #     mech1_str, mech2_str, temps)
        mech1_ktp_dct, mech2_ktp_dct = calculator.combine.build_reaction_name_dcts(
            mech1_str, mech2_str,
            t_ref, temps, pressures,
            remove_bad_fits=True)
    elif dct_idx == 'inchi':
        _, mech2_thermo_dct = calculator.combine.build_thermo_inchi_dcts(
            mech1_str, mech2_str, mech1_csv_str, mech2_csv_str, temps)
        mech1_ktp_dct, mech2_ktp_dct = calculator.combine.build_reaction_inchi_dcts(
            mech1_str, mech2_str, mech1_csv_str, mech2_csv_str,
            t_ref, temps, pressures,
            remove_bad_fits=True)

    # print('\n\nmech2 thermo')
    # for spc in mech2_thermo_dct:
    #     print(spc)
    #     print(mech2_thermo_dct[spc])

    print('\npre flip mech1')
    for rxn in mech1_ktp_dct:
        print(rxn)
        print(mech1_ktp_dct[rxn])
    # print('\npre flip mech2')
    # for rxn in mech2_ktp_dct:
    #     print(rxn)
    #     print(mech2_ktp_dct[rxn])

    # efe
    # print('\n\n\n')
    # ini_mech1_names = list(mech1_ktp_dct.keys())
    # mech2_names = list(mech2_ktp_dct.keys())
    # mech1_names = []
    # for i in range(len(mech2_names)):
    #     if i+1 <= len(ini_mech1_names):
    #         mech1_names.append(ini_mech1_names[i])
    #     else:
    #         mech1_names.append('AAA')

    # for m1, m2 in zip(mech1_names, mech2_names):
    #     print(m1, '   ', m2)
    # import sys
    # sys.exit()

    # Combine the two thermo dictionaries together under a common index
    # mech2_thermo_dct=None
    ktp_dct = calculator.combine.mechanism_rates(
        mech1_ktp_dct, mech2_ktp_dct, temps,
        mech2_thermo_dct=mech2_thermo_dct,
        ignore_reverse=ignore_reverse
    )

    # import sys
    # sys.exit()

    # Build list of CHEMKIN mechanism names for plotting
    if dct_idx != 'name':
        ktp_dct = calculator.combine.conv_ich_to_name_ktp_dct(
            ktp_dct, mech1_csv_str)
    names = list(ktp_dct.keys())

    print('mech dcts')
    for rxn, dcts in ktp_dct.items():
        print(rxn)
        print('mech1')
        print(dcts['mech1'])
        print('mech2')
        print(dcts['mech2'])

    # Plot the data from both mechanisms for each species
    plotter.rates.build(
        ktp_dct, temps, dir_prefix=dir_prefix,
        names=names, mech_labels=mech_labels)


def thermo(mech1_str, mech2_str, temps,
           mech1_csv_str=None, mech2_csv_str=None,
           dct_idx='name', dir_prefix='therm_plots'):
    """ Read the thermo sections of two CHEMKIN files and
        produce plots of the thermochemical parameters together
    """

    # Build dictionaries containing the thermo data strings for each species
    if dct_idx == 'name':
        mech1_thermo_dct, mech2_thermo_dct = calculator.combine.build_thermo_name_dcts(
            mech1_str, mech2_str, temps)
    elif dct_idx == 'inchi':
        mech1_thermo_dct, mech2_thermo_dct = calculator.combine.build_thermo_inchi_dcts(
            mech1_str, mech2_str, mech1_csv_str, mech2_csv_str, temps)

    # Combine the two thermo dictionaries together under a common index
    thermo_vals_dct = calculator.combine.mechanism_thermo(
        mech1_thermo_dct, mech2_thermo_dct)

    # Build list of CHEMKIN mechanism names for plotting
    if dct_idx == 'name':
        names = list(thermo_vals_dct.keys())
    elif dct_idx == 'inchi':
        names = [calculator.combine.spc_name_from_inchi(mech1_csv_str, mech2_csv_str, key)
                 for key in thermo_vals_dct]

    # Plot the data from both mechanisms for each species
    plotter.thermo.build(
        thermo_vals_dct, temps, dir_prefix=dir_prefix, names=names)
