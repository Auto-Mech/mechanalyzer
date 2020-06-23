"""
Read the mechanism file
"""

import automol
import chemkin_io
from lib.amech_io.parser import ptt


ALLOWED_HEADERS = [
    'smiles', 'inchi', 'inchikey', 'mult', 'charge', 'sens'
]


# Build various dictionaries
def csv_dct(csv_str, values, key='name'):
    """ read the species file in a .csv format

        :param csv_str: string for species input file
        :type csv_str: str
        :param values: 
        :type values: lst(str)
        :param key: key to index the spc dictionary
        :type key: 
        :return: dict[key] = val
    """

    assert set(values) <= set(ALLOWED_HEADERS), (
        'Unallowed dict values requested.'
    )
    assert key in ('name', 'inchi'), (
        'Unallowed dict key requested.'
    )

    # Check if the csv file is formatted properly
    csv_proper = _check_csv_file(csv_str)

    if csv_proper:
        # Read the names from the file
        names = _read_csv_names(csv_str)

        # Read in the initial CSV file, 
        value_dcts = []
        for value in values:
            value_dct.append(
                chemkin_io.parser.mechanism.spc_name_dct(csv_str, value))

        # Build the final dictionary
        spc_dct = {}
        for name in names:
            spc_dct[name] = {}
            for value_dct in value_dcts:
                spc_dct[name][value] = value_dct.get(name, None)
    else:
        spc_dct = None

    return spc_dct


# def spc_name_dct(csv_str, entry):
#     """ Read the species.csv file and generate a dictionary that relates
#         structural information to the ChemKin mechanism name.
# 
#         :param csv_str: string of input csv file with species information
#         :type csv_str: str
#         :param entry: structural information that is desired
#         :type entry: str
#         :return spc_dct: all species with desired structural information
#         :rtype: dict[name:entry]
#     """
# 
#     data = _read_csv(csv_str)
# 
#     if entry == 'inchi':
#         spc_dct = _read_name_inchi(data)
#     elif entry == 'smiles':
#         spc_dct = _read_name_smiles(data)
#     elif entry == 'mult':
#         spc_dct = _read_name_mult(data)
#     elif entry == 'charge':
#         spc_dct = _read_name_charge(data)
#     elif entry == 'sens':
#         spc_dct = _read_name_sensitivity(data)
#     else:
#         raise NotImplementedError
# 
#     return spc_dct
# 
# 
# def spc_inchi_dct(csv_str):
#     """ Read the species.csv file and generate a dictionary that relates
#         ChemKin mechanism name to InChI string.
# 
#         :param csv_str: string of input csv file with species information
#         :type csv_str: str
#         :return spc_dct: all species with names and InChI strings
#         :rtype: dict[InChI: name]
#     """
# 
#     data = _read_csv(csv_str)
# 
#     spc_dct = {}
#     if hasattr(data, 'inchi'):
#         spc_dct = dict(zip(data.name, data.inchi))
#     elif hasattr(data, 'smiles'):
#         ichs = [_inchi(smiles) for smiles in data.smiles]
#         spc_dct = dict(zip(ichs, data.name))
#     else:
#         spc_dct = {}
# 
#     return spc_dct


def _read_csv_names(csv_str):
    """ Obtain a list of names from a name 
    """

    data = _read_csv(csv_str)

    if hasattr(data, 'name'):
        names = list(data.name)
    else:
        names = None

    return names


def _read_name_inchi(data):
    """ Build the species dictionary relating ChemKin name to InChI string.
        The InChI strings are read directly from the data object if available.
        Otherwise they are generated using the SMILES strings.

        :param data: information from input species.csv file
        :type data: pandas
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: InChI]
    """

    if hasattr(data, 'inchi'):
        spc_dct = dict(zip(data.name, data.inchi))
    elif hasattr(data, 'smiles'):
        print('No inchi column in csv file, getting inchi from SMILES')
        ichs = [_inchi(smiles) for smiles in data.smiles]
        spc_dct = dict(zip(data.name, ichs))
    else:
        spc_dct = {}
        print('No "inchi" or "SMILES" column in csv file')

    # Fill remaining inchi entries if inchi
    for i, name in enumerate(data.name):
        if str(spc_dct[name]) == 'nan':
            spc_dct[name] = _inchi(data.smiles[i])

    return spc_dct


def _read_name_smiles(data):
    """ Build the species dictionary relating ChemKin name to SMILES string.
        The SMILES strings are read directly from the data object if available.
        Otherwise they are generated using the InChI strings.

        :param data: information from input species.csv file
        :type data: pandas
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: SMILES]
    """

    spc_dct = {}
    if hasattr(data, 'smiles'):
        spc_dct = dict(zip(data.name, data.smiles))
    elif hasattr(data, 'inchi'):
        smiles = [_smiles(ich) for ich in data.inchi]
        spc_dct = dict(zip(data.name, smiles))
    else:
        spc_dct = {}
        print('No "SMILES" or "InChI" column in csv file')

    return spc_dct


def _read_name_mult(data):
    """ Build the species dictionary relating ChemKin name to multiplicity.

        :param data: information from input species.csv file
        :type data: pandas
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: multiplicity]
    """

    if hasattr(data, 'mult'):
        spc_dct = dict(zip(data.name, data.mult))
    else:
        spc_dct = {}
        print('No "mult" column in csv file')

    return spc_dct


def _read_name_charge(data):
    """ Build the species dictionary relating ChemKin name to charge.
        If the charge is missing for a given species, the dictionary
        element will be set to zero, assuming a neutral species.

        :param data: information from input species.csv file
        :type data: pandas
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: charge]
    """

    fill = 0
    if hasattr(data, 'charge'):
        spc_dct = dict(zip(data.name, data.charge))
    else:
        if fill is not None:
            spc_dct = dict(zip(data.name, [fill for name in data.name]))
        else:
            spc_dct = {}

    return spc_dct


def _read_name_sensitivity(data):
    """ Build the species dictionary relating ChemKin name to sensitivity.
        If the sensitivity is missing for a given species, the dictionary
        element will be set to zero.

        :param data: information from input species.csv file
        :type data: pandas
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: sensitivity]
    """

    fill = 0.
    if hasattr(data, 'sens'):
        spc_dct = dict(zip(data.name, data.sens))
    else:
        if fill is not None:
            spc_dct = dict(zip(data.name, [fill for name in data.name]))
        else:
            spc_dct = {}

    return spc_dct


def spc_inchi_dct(csv_str):
    """ Read the species.csv file and generate a dictionary that relates
        ChemKin mechanism name to InChI string.

        :param csv_str: string of input csv file with species information
        :type csv_str: str
        :return spc_dct: all species with names and InChI strings
        :rtype: dict[InChI: name]
    """

    data = _read_csv(csv_str)

    spc_dct = {}
    if hasattr(data, 'inchi'):
        spc_dct = dict(zip(data.name, data.inchi))
    elif hasattr(data, 'smiles'):
        ichs = [_inchi(smiles) for smiles in data.smiles]
        spc_dct = dict(zip(ichs, data.name))
    else:
        spc_dct = {}

    return spc_dct


def _read_csv(csv_str):
    """ Read the csv file and generate data using pandas.

        :param csv_str: string of input csv file with species information
        :type csv_str: str
        :return data: data for the species from csv file
        :rtype: pandas.data object?
    """

    # Read in csv file while removing whitespace and make all chars lowercase
    csv_file = StringIO(csv_str)
    data = pandas.read_csv(csv_file, comment='#', quotechar="'")

    # Parse CSV string into data columns
    # data.columns = data.columns.str.strip()
    data.columns = map(str.lower, data.columns)

    return data


def _check_csv(csv_str):
    """ Check if the csv is formatted properly
    """

    data = _read_csv(csv_str)
    if set(list(data)) <= set(ALLOWED_HEADERS):
        proper = True
    else:
        print('Unallowed Headers in Header file')
        proper = False

    return proper
