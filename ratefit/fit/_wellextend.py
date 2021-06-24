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
    """ Obtain the energies for each well at the given pressure.
    """

    mess_temps, _ = mess_io.reader.rates.temperatures(mess_out_str)
    max_run_temp = max(mess_temps)

    # Get the temps where each well exists
    well_enes = {}
    well_rxns = _get_well_reactions(mess_out_str)
    for well, rxn_lst in well_rxns.items():
        print('\n***********************************************\n')
        print('Obtaining information for well {} at P={}'.format(
            well, pressure))
        max_temp = -1.0
        for (lab_i, lab_j) in rxn_lst:

            rxn = '{}->{}'.format(lab_i, lab_j)

            # Read the rate constants out of the mess outputs
            print('\n-----------------------------------------------\n')
            print('Finding max T where k(T) exists for {}...'.format(rxn))
            ktp_dct, _ = read_rates(mess_out_str, None, lab_i, lab_j)
            rxn_temp = _max_temp_well_exists(ktp_dct, pressure, mess_temps)

            if rxn_temp > max_temp:
                max_temp = rxn_temp
                print('- New max temperature for well: {} K'.format(rxn_temp))

        # Determine if k(T) exist at highest T -> no Well cap exists
        if numpy.isclose(max_temp, max_run_temp):
            max_temp = None
            print('\nMax temperature at highest value from run.',
                  'No cap needed.')
        else:
            print('\nMax temperature for energies is {} K'.format(max_temp))

        # Read the average energy at the max temperature
        ene = (mess_io.reader.well_average_energy(mess_log_str, well, max_temp)
               if max_temp is not None else None)
        well_enes[well] = ene

    return well_enes


def _get_well_reactions(mess_out_str):
    """ Get the reactions for each wells
    """

    rxn_pairs, _ = gen_reaction_pairs(mess_out_str, None)

    # Get the well labels from the reactions
    rgts = tuple(spc for spc in itertools.chain(*rxn_pairs) if 'W' in spc)
    wells = tuple(n for i, n in enumerate(rgts) if n not in rgts[:i])

    # Grab reactions that contains the well as the reactant
    well_rxns = {}
    for well in wells:
        rxn_lst = ()
        for rxn in rxn_pairs:
            rct, _ = rxn
            if well == rct:
                rxn_lst += (rxn,)
        well_rxns[well] = rxn_lst

    return well_rxns


def _max_temp_well_exists(ktp_dct, pressure, mess_temps):
    """ For a given reaction and pressure, find a max temperature

        Assumes a filtered ktp dct.
    """

    max_temp = None
    for _pressure, kt_lst in ktp_dct.items():
        if _pressure != 'high':
            if numpy.isclose(_pressure, pressure):
                max_temp = max(kt_lst[0])
                break

    # default to the lowest temperature run if no rate constants found
    if max_temp is None:
        max_temp = min(mess_temps)
        print('\nNo k(T) values found for P = {} atm.'.format(pressure))
        print('T={}: minimum of all temps in output'.format(max_temp))
    else:
        print('\nT={}: max T where k(T) found'.format(max_temp))

    return max_temp


# Handlies writing the new string
def _format_well_extension_inp(inp_str, well_enes_dct, well_lump_str):
    """ handles building new input will well lumping/extension info
    """

    # Reinitialize string
    new_inp_str = inp_str

    # Write string for each of the well enes
    for well, ene in well_enes_dct.items():
        if ene is not None:
            # Find line for where well start, for-loop handle weird format
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
