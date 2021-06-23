""" Library of functions to handle various procedures
"""

import itertools
import numpy
import mess_io
import ioformat
# import autorun
from phydat import phycon
from ratefit.fit._read import gen_reaction_pairs
from ratefit.fit._read import read_rates


# Main calling function
# def well_lumped_input_file(script_str, run_dir, mess_inp_str,
#                            aux_dct, input_name,
#                            lump_pressure, lump_temp):
def well_lumped_input_file(inp_str, out_str, aux_str, log_str,
                           lump_pressure, lump_temp):
    """ Run MESS to get the wells and then parse the aux file for wells...
    """

    # Run MESS with input with no lumping specified
    # output_strs = autorun.mess.direct(
    #     script_str, run_dir, mess_inp_str,
    #     aux_dct=aux_dct,
    #     input_name=input_name)
    # out_str,aux_str,log_str=output_strs[0], output_strs[1], output_strs[2]

    # Get the requisite information from analyzing the output
    well_enes_dct = well_energies(out_str, log_str, lump_pressure)
    well_lump_str = well_lumping_scheme(
        aux_str, lump_pressure, lump_temp)

    # Write a new string containing the parsed information
    well_extend_str = _format_well_extension_inp(
        inp_str, well_enes_dct, well_lump_str)

    return well_extend_str


# Build the lumping scheme used for the well extension
def well_lumping_scheme(mess_aux_str, pressure, temp):
    """ Parse lumped wells from aux output; write into string for new input
    """

    well_lump_lst = mess_io.reader.merged_wells(mess_aux_str, pressure, temp)
    well_lump_str = mess_io.writer.well_lump_scheme(well_lump_lst)

    return well_lump_str


# Get the energies for definining the well extension cap
def well_energies(mess_out_str, mess_log_str, pressure):
    """ Analyze a well lumping algorithm.
    """

    # Get the temps where each well exists
    well_enes = {}
    well_rxns = _get_well_reactions(mess_out_str)
    for well, (lab_i, lab_j) in well_rxns.items():

        print('Obtaining information for well {} at P={}'.format(
            well, pressure))

        # Read the rate constants out of the mess outputs
        print('\nReading k(T,P)s from MESS output...')
        ktp_dct, _ = read_rates(mess_out_str, None, lab_i, lab_j)
        max_temp = _max_temp_well_exists(ktp_dct, pressure)

        print('\nMax T for energies is {}'.format(max_temp))

        # Read the average energy at the max temperature
        well_enes[well] = (
            mess_io.reader.well_average_energy(mess_log_str, well, max_temp)
        )

    return well_enes


def _get_well_reactions(mess_out_str):
    """ Get the reactions for each wells
    """

    rxn_pairs, _ = gen_reaction_pairs(mess_out_str, None)

    # Get the well labels from the reactions
    rgts = tuple(spc for spc in itertools.chain(*rxn_pairs) if 'W' in spc)
    wells = tuple(n for i, n in enumerate(rgts) if n not in rgts[:i])

    # Grab a reaction that contains the well
    well_rxns = {}
    for well in wells:
        for rxn in rxn_pairs:
            rct, _ = rxn
            if well == rct:
                well_rxns[well] = rxn
                break

    return well_rxns


def _max_temp_well_exists(ktp_dct, pressure):
    """ For a given reaction and pressure, find a max temperature

        Assumes a filtered ktp dct.
    """

    max_temp = 600.0  # fake value for now
    for _pressure, kt_lst in ktp_dct.items():
        if _pressure != 'high':
            if numpy.isclose(_pressure, pressure):
                max_temp = max(kt_lst[0])
                break

    return max_temp


# Handlies writing the new string
def _format_well_extension_inp(inp_str, well_enes_dct, well_lump_str):
    """ handles building new input will well lumping/extension info
    """

    # Reinitialize string
    new_inp_str = inp_str

    # Write string for each of the well enes
    for well, ene in well_enes_dct.items():
        # Find line for where well start, for-loop used for formatting diffs
        for line in inp_str.splitlines():
            if 'Well' in line and well in line:
                _search = line
                break
        _add = '  WellExtensionCap[kcal/mol]    {0:.2f}'.format(
            ene*phycon.EH2KCAL)
        new_inp_str = ioformat.add_line(
            string=new_inp_str, addline=_add,
            searchline=_search, position='after')

    # Write new strings with the lumped input
    well_extend_line = 'WellExtension\nExtensionCorrection    0.2'
    well_lump_line = ioformat.indent(well_lump_str, 2)
    new_inp_str = ioformat.add_line(
        string=new_inp_str, addline=well_extend_line,
        searchline='Model', position='before')
    new_inp_str = ioformat.add_line(
        string=new_inp_str, addline=well_lump_line,
        searchline='Model', position='after')

    return new_inp_str
