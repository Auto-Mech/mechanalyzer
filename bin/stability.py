"""
  Determine the stability of radicals
"""

import os
import pandas as pd
import chemkin_io
from _format import read_file
from _format import format_rxn
from _format import formatp


# def read_radical_stabilities1():
#     """ Read the reaction sections of a CHEMKIN files and
#         calculate the rate constants
#     """
#
#     # Loop over each pressure
#     for pressure in PRESSURES:
#
#         # Initialize dct for each pressure
#         stab_dct = {}
#
#         # Loop over each of the ckin files
#         for ckin_idx, ckin_file in enumerate(CKIN_FILES):
#
#             # Read mechanism files
#             mech_str = _read_file(os.path.join(PATH, ckin_file))
#
#             # Build the reactions data dct
#             rxn_block = chemkin_io.parser.mechanism.reaction_block(mech_str)
#             rxn_dct = chemkin_io.parser.reaction.data_dct(rxn_block)
#
#             # Loop over the reaction dct
#             temp = None
#             for rxn, dstr in rxn_dct.items():
#
#                 # Unpack to get the reactants and products, and radical
#                 [reacs, prods] = rxn
#                 radical = prods[0]
#
#                 # ts, grab from addns
#                 if len(reacs) == 2 and len(prods) == 1:
#
#                     # Get the radical name
#                     radical = prods[0]
#
#                     # Get the max temps
#                     inf_dct = chemkin_io.parser.reaction.ratek_fit_info(dstr)
#                     temp = inf_dct[pressure]['temps'][1]
#
#                     # Add temperature to the stability dct
#                     if radical not in stab_dct:
#                         stab_dct[radical] = {ckin_idx: temp}
#                     else:
#                         stab_dct[radical].update({ckin_idx: temp})
#
#         # Print the table
#         print('\n\nPressure: {}'.format(pressure))
#
#         # radicals.sort(key=lambda x: x[1])
#         for radical in stab_dct:
#
#             # Print the header string
#             header_str = '{0:<25s}'.format('Radical')
#             for idx, _ in enumerate(CKIN_FILES):
#                 header_str += '{0:>6s}'.format('lvl'+str(idx+1))
#             print(header_str)
#
#             # Print the temperatures for the radical
#             rad_str = '{0:<25s}'.format(radical)
#             for idx, _ in enumerate(CKIN_FILES):
#                 if idx in stab_dct[radical]:
#                     rad_str += '{0:>6d}'.format(
#                         stab_dct[radical][idx])
#                 else:
#                     rad_str += '{0:>6s}'.format('')
#             rad_str += '\n'
#             print(rad_str)


def read_radical_stabilities2(pressures, ckin_files, path):
    """ Read the reaction sections of a CHEMKIN files and
        calculate the rate constants
    """

    # Loop over each pressure
    for ckin_idx, ckin_file in enumerate(ckin_files):

        # Print the table
        print('\n\nLevel: {}'.format(str(ckin_idx+1)))

        # Initialize branching dict for this level
        stab_dct = {}

        # Read mechanism files
        mech_str = read_file(os.path.join(path, ckin_file))

        # Build the reactions data dct
        rxn_block = chemkin_io.parser.mechanism.reaction_block(mech_str)
        rxn_dct = chemkin_io.parser.reaction.data_dct(
            rxn_block, remove_bad_fits=True)

        # Loop over each of the ckin files
        for pressure in pressures:

            # Initialize dct for each pressure
            pstr = '{}'.format(formatp(pressure))
            stab_dct[pstr] = {}

            # Loop over the reaction dct
            temp = None
            for rxn, dstr in rxn_dct.items():

                # Unpack to get the reactants and products, and radical
                [reacs, prods] = rxn
                radical = prods[0]

                # ts, grab from addns
                if len(reacs) == 2 and len(prods) == 1:
                    # Get the radical name
                    radical = prods[0]

                    # Get the max temps
                    inf_dct = chemkin_io.parser.reaction.ratek_fit_info(dstr)
                    if pressure in inf_dct:
                        temp = inf_dct[pressure]['temps'][1]
                    else:
                        print('No pressure in inf dct for', format_rxn(rxn))
                        continue

                    # Build dictionary for data frame
                    if radical in stab_dct[pstr]:
                        # update temp in dct if new temp is lower
                        if temp < stab_dct[pstr][radical]:
                            print(radical, temp, stab_dct[pstr][radical])
                            stab_dct[pstr].update({radical: temp})
                    else:
                        stab_dct[pstr].update({radical: temp})

        # Build Pandas data frame
        branch_df = pd.DataFrame(stab_dct)
        print(branch_df)
