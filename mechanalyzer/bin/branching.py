""" Calculate branching ratios for mechanisms
"""

import os
import pandas as pd
import mechanalyzer
from _format import read_file
from _format import format_rxn
from _format import chk_rxn


# Set printing options
pd.set_option('display.max_columns', 5)
pd.set_option('display.max_rows', 500)


def calc_multimech_rates(temps, pressures, t_ref, ckin_files, path,
                         rtyp='ignore', allow_rcts=()):
    """ Calculate rate constants across multiple mechanisms
    """

    # Loop over each pressure
    for pressure in pressures:

        # Initialize dct for each pressure
        full_ratek_dct = {}
        full_branch_dct = {}

        # Loop over each of the ckin files
        for ckin_idx, ckin_file in enumerate(ckin_files):

            # Initialize branching dict for this level
            lvl_str = '{0:4s}'.format('lvl'+str(ckin_idx+1))
            full_ratek_dct[lvl_str] = {}
            full_branch_dct[lvl_str] = {}

            # Read mechanism files
            mech_str = read_file(os.path.join(path, ckin_file))

            # Build the reactions total and branching dcts
            rxn_block = chemkin_io.parser.mechanism.reaction_block(mech_str)
            rxn_units = chemkin_io.parser.mechanism.reaction_units(mech_str)
            ratek_dct = mechanalyzer.calculator.rates.mechanism(
                rxn_block, rxn_units, t_ref, temps, pressures=[pressure],
                ignore_reverse=True, remove_bad_fits=True)
            branch_dct, _ = mechanalyzer.calculator.rates.branching_fractions(
                ratek_dct, [pressure])

            # Add to the overall dcts
            for rxn in branch_dct:
                # Only grab certain reactions
                if chk_rxn(rxn, rtyp, allow_rcts=allow_rcts):
                    # Build dictionary for data frame
                    bfs = branch_dct[rxn][pressure]
                    if bfs is not None:
                        rxn_str = format_rxn(rxn)
                        full_branch_dct[lvl_str].update({rxn_str: bfs[0]})

            # Add to the overall dcts
            for rxn in ratek_dct:
                # Only grab certain reactions
                if chk_rxn(rxn, rtyp, allow_rcts=allow_rcts):
                    # Build dictionary for data frame
                    rxn_str = format_rxn(rxn)
                    if pressure in ratek_dct[rxn]:
                        full_ratek_dct[lvl_str].update(
                            {rxn_str: ratek_dct[rxn][pressure][0]})

        # Print the table
        print('\n\nPressure: {}'.format(pressure))

        # Build Pandas data frame
        ktp_df = pd.DataFrame(full_ratek_dct)
        branch_df = pd.DataFrame(full_branch_dct)

        print(ktp_df)
        print('\n')
        print(branch_df)


def calc_multimech_rates2(temps, pressures, t_ref, ckin_files, path,
                          rtyp='ignore', allow_rcts=()):
    """ Calculate rate constants across multiple mechanisms
    """

    # Loop over each pressure
    for ckin_idx, ckin_file in enumerate(ckin_files):

        # Initialize branching dict for this level
        lvl_str = '{0:4s}'.format('lvl'+str(ckin_idx+1))
        full_ratek_dct = {}
        full_branch_dct = {}

        # Loop over each of the ckin files
        for pressure in pressures:

            # Initialize dct for each pressure
            full_ratek_dct[pressure] = {}
            full_branch_dct[pressure] = {}

            # Read mechanism files
            mech_str = read_file(os.path.join(path, ckin_file))

            # Build the reactions total and branching dcts
            rxn_block = chemkin_io.parser.mechanism.reaction_block(mech_str)
            rxn_units = chemkin_io.parser.mechanism.reaction_units(mech_str)
            ratek_dct = mechanalyzer.calculator.rates.mechanism(
                rxn_block, rxn_units, t_ref, temps, pressures=[pressure],
                ignore_reverse=True, remove_bad_fits=True)
            branch_dct, _ = mechanalyzer.calculator.rates.branching_fractions(
                ratek_dct, [pressure])

            # Add to the overall dcts
            for rxn in branch_dct:
                # Only grab certain reactions
                if chk_rxn(rxn, rtyp, allow_rcts=allow_rcts):
                    # Build dictionary for data frame
                    bfs = branch_dct[rxn][pressure]
                    if bfs is not None:
                        rxn_str = format_rxn(rxn)
                        full_branch_dct[pressure].update({rxn_str: bfs[0]})

            # Add to the overall dcts
            for rxn in ratek_dct:
                # Only grab certain reactions
                if chk_rxn(rxn, rtyp, allow_rcts=allow_rcts):
                    # Build dictionary for data frame
                    rxn_str = format_rxn(rxn)
                    if pressure in ratek_dct[rxn]:
                        full_ratek_dct[pressure].update(
                            {rxn_str: ratek_dct[rxn][pressure][0]})

        # Print the table
        print('\n\nlevel: {}'.format(lvl_str))

        # Build Pandas data frame
        ktp_df = pd.DataFrame(full_ratek_dct)
        branch_df = pd.DataFrame(full_branch_dct)

        print(ktp_df)
        print('\n')
        print(branch_df)
