""" Functions operating on Chemkin input files or strings
"""
import ioformat.pathtools as parser
from chemkin_io.parser import mechanism as parser_mech
from chemkin_io.parser import reaction as parser_rxn
from chemkin_io.parser import thermo as parser_thermo
from chemkin_io.parser import species as parser_spc
from mechanalyzer.calculator import rates as calc_rates
from mechanalyzer.calculator import thermo as calc_thermo


def load_rxn_ktp_dcts(mech_filenames, direc, temps, pressures):
    """ Read Chemkin mechanism files and calculate rates at the indicated
        pressures and temperatures. Return a list of rxn_ktp_dcts.

        :param mech_filenames: Chemkin mechanism filenames
        :type mech_filenames: list [filename1, filename2, ...]
        :param direc: directory with file(s) (all must be in same directory)
        :type direc: str 
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :param pressures: pressures at which to do calculations (atm)
        :type pressures: list [float]
        :return rxn_ktp_dcts: list of rxn_ktp_dcts
        :rtype: list of dcts [rxn_ktp_dct1, rxn_ktp_dct2, ...]
    """ 
    rxn_ktp_dcts = []
    for mech_filename in mech_filenames:
        rxn_ktp_dct = load_rxn_ktp_dct(mech_filename, direc, temps, pressures)
        rxn_ktp_dcts.append(rxn_ktp_dct)

    return rxn_ktp_dcts


def load_rxn_param_dcts(mech_filenames, direc):
    """ Read Chemkin-formatted mechanism files and return a list of
        rxn_param_dcts.

        :param mech_filenames: Chemkin mechanism filenames
        :type mech_filenames: list [filename1, filename2, ...]
        :param direc: directory with file
        :type direc: str 
        :return rxn_param_dcts: list of rxn_param_dcts 
        :rtype: list of dcts [rxn_param_dct1, rxn_param_dct2, ...]
    """ 
    rxn_param_dcts = []
    for mech_filename in mech_filenames:
        rxn_param_dct = load_rxn_ktp_dct(mech_filename, direc)
        rxn_param_dcts.append(rxn_param_dct)

    return rxn_param_dcts


def load_spc_therm_dcts(thermo_filenames, direc, temps):
    """ Reads Chemkin thermo files and calculates thermo at the indicated
        temperatures. Outputs a list of spc_therm_dcts.

        :param thermo_filename: filenames containing Chemkin thermo information
        :type thermo_filenames: list [filename1, filename2, ...]
        :param direc: directory with file(s) (all must be in same directory)
        :type direc: str
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :return spc_therm_dcts: list of spc_therm_dcts
        :rtype: list of dcts [spc_therm_dct1, spc_therm_dct2, ...]
    """
    spc_therm_dcts = []
    for thermo_filename in thermo_filenames:
        print('thermo filename test:', thermo_filename)
        spc_therm_dct = load_spc_therm_dct(thermo_filename, direc, temps)
        spc_therm_dcts.append(spc_therm_dct)

    return spc_therm_dcts


def load_rxn_ktp_dct(mech_filename, direc, temps, pressures):
    """ Read a Chemkin-formatted mechanism file and calculate rates at the indicated
        pressures and temperatures. Return a rxn_ktp_dct.

        :param mech_filename: Chemkin mechanism filename
        :type mech_filename: str
        :param direc: directory with file
        :type direc: str 
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :param pressures: pressures at which to do calculations (atm)
        :type pressures: list [float]
        :return rxn_ktp_dct: rxn_ktp_dct object
        :rtype: dct {rxn1: ktp_dct1, rxn2: ...} 
    """ 
    rxn_param_dct = load_rxn_param_dct(mech_filename, direc)
    rxn_ktp_dct = calc_rates.eval_rxn_param_dct(rxn_param_dct, pressures, temps)

    return rxn_ktp_dct


def load_rxn_param_dct(mech_filename, direc):
    """ Read a Chemkin-formatted mechanism file and return a rxn_param_dct.

        :param mech_filename: Chemkin mechanism filename
        :type mech_filename: str
        :param direc: directory with file
        :type direc: str 
        :return rxn_param_dct: rxn_param_dct object
        :rtype: dct {rxn1: param_tuple1, rxn2: ...} 
    """ 
    mech_str = parser.read_file(direc, mech_filename)
    rxn_param_dct = parse_rxn_param_dct(mech_str)

    return rxn_param_dct


def load_spc_therm_dct(thermo_filename, direc, temps):
    """ Reads a Chemkin thermo file and calculates thermo at the indicated
        temperatures. Outputs a spc_therm_dct.

        :param thermo_filename: filename containing Chemkin thermo information
        :type thermo_filename: str
        :param direc: directory with file
        :type direc: str
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :return spc_therm_dct: spc_therm_dct object
        :rtype: dct {spc1: therm_array1, spc2: ...}
    """
    mech_str = parser.read_file(direc, thermo_filename)
    spc_therm_dct = parse_spc_therm_dct(mech_str, temps)

    return spc_therm_dct


def parse_rxn_ktp_dct(mech_str, temps, pressures):
    """ Parses a raw Chemkin mechanism string and yields a rxn_ktp_dct

        :param mech_str: raw string from reading a Chemkin file
        :type mech_str: str
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :param pressures: pressures at which to do calculations (atm)
        :type pressures: list [float]
        :return rxn_ktp_dct: rxn_ktp_dct object
        :rtype: dct {rxn1: ktp_dct1, rxn2: ...} 
    """
    rxn_param_dct = parse_rxn_param_dct(mech_str)
    rxn_ktp_dct = calc_rates.eval_rxn_param_dct(rxn_param_dct, pressures, temps)

    return rxn_ktp_dct


def parse_rxn_param_dct(mech_str):
    """ Parses a raw Chemkin mechanism string and yields a rxn_param_dct

        :param mech_str: raw string from reading a Chemkin file
        :type mech_str: str
        :return rxn_param_dct: rxn_param_dct object
        :rtype: dct {rxn1: param_tuple1, rxn2: ...} 
    """
    ea_units, a_units = parser_mech.reaction_units(mech_str)
    rxn_block_str = parser_mech.reaction_block(mech_str)
    rxn_param_dct = parser_rxn.param_dct(rxn_block_str, ea_units, a_units)
     
    return rxn_param_dct


def parse_spc_therm_dct(mech_str, temps):
    """ Parses a raw Chemkin mechanism string and yields a spc_therm_dct

        *note: the input mech_str can be from reading the entire Chemkin file,
            not just the thermo portion

        :param mech_str: raw string from reading a Chemkin file
        :type mech_str: str
        :param temps: temperatures at which to do calculations (Kelvin)
        :type temps: list [float]
        :return rxn_param_dct: rxn_param_dct object
        :rtype: dct {rxn1: param_tuple1, rxn2: ...} 
    """
    thermo_block_str = parser_mech.thermo_block(mech_str)
    spc_nasa7_dct = parser_thermo.create_spc_nasa7_dct(thermo_block_str)
    spc_therm_dct = calc_thermo.create_spc_therm_dct(spc_nasa7_dct, temps)

    return spc_therm_dct


def parse_elem_tuple(mech_str):
    """ Parses a raw Chemkin mechanism string and yields an elem_tuple

        :param mech_str: raw string from reading a Chemkin file
        :type mech_str: str
        :return elem_tuple: a tuple of elements in the mechanism
        :rtype: tuple (elem1, elem2, ...)
    """
    el_block = parser_mech.element_block(mech_str)
    elem_tuple = parser_spc.names(el_block)

    return elem_tuple

