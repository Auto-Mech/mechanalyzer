"""
Read the csv file
"""


from io import StringIO
import pandas
from automol.smiles import inchi as _inchi
from automol.inchi import smiles as _smiles


# What columns are allowed in the CSV file
ALLOWED_HEADERS = (
    'name', 'smiles', 'inchi', 'inchikey', 'mult', 'charge', 'sens'
)
DEFAULT_HEADERS = (
    'smiles', 'inchi', 'inchikey', 'mult', 'charge', 'sens'
)


# Build various dictionaries
def csv_dct(csv_str, values=DEFAULT_HEADERS, key='name'):
    """ read the species file in a .csv format

        :param csv_str: string for species input file
        :type csv_str: str
        :param values: values to build dcts for
        :type values: lst(str)
        :param key: key to index the spc dictionary
        :type key: str
        :return: dict[key] = val
    """

    assert set(values) <= set(ALLOWED_HEADERS), (
        'Unallowed dict values requested.'
    )
    assert key in ('name', 'inchi'), (
        'Unallowed dict key requested.'
    )

    # Read the CSV file
    data = _read_csv(csv_str)

    # Check if the csv file is formatted properly
    csv_proper = _check_csv(data)

    if csv_proper:
        # Read the names from the file
        names = _read_csv_names(data)

        # Read in the initial CSV file,
        value_dcts = []
        for value in values:
            value_dcts.append(READERS[value](data, names))

        # Build the final dictionary
        spc_dct = {}
        for name in names:
            spc_dct[name] = {}
            for value_idx, value_dct in enumerate(value_dcts):
                spc_dct[name][values[value_idx]] = value_dct.get(name, None)
    else:
        spc_dct = None

    return spc_dct


def _read_csv_names(data):
    """ Obtain a list of names from a name
    """

    if hasattr(data, 'name'):
        names = list(data.name)
    else:
        names = None

    return names


def _read_csv_inchi(data, idxs):
    """ Build the species dictionary relating ChemKin name to InChI string.
        The InChI strings are read directly from the data object if available.
        Otherwise they are generated using the SMILES strings.

        :param data: information from input species.csv file
        :type data: pandas
        :param idxs: index for the dict to be created
        :type idxs: list(str)
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: InChI]
    """
    if hasattr(data, 'inchi'):
        spc_dct = dict(zip(idxs, data.inchi))
    elif hasattr(data, 'smiles'):
        print('No inchi column in csv file, getting inchi from SMILES')
        ichs = [_inchi(smiles) for smiles in data.smiles]
        spc_dct = dict(zip(idxs, ichs))
    else:
        spc_dct = {}
        print('No "inchi" or "SMILES" column in csv file')

    # Fill remaining inchi entries if inchi
    for i, name in enumerate(idxs):
        if str(spc_dct[name]) == 'nan':
            spc_dct[name] = _inchi(data.smiles[i])

    return spc_dct


def _read_csv_inchikey(data, idxs):
    """ Build the species dictionary relating ChemKin name to InChI key.
    """

    if hasattr(data, 'inchikey'):
        spc_dct = dict(zip(idxs, data.inchikey))
    else:
        spc_dct = dict(zip(idxs, ['' for _ in range(len(idxs))]))

    return spc_dct


def _read_csv_smiles(data, idxs):
    """ Build the species dictionary relating ChemKin name to SMILES string.
        The SMILES strings are read directly from the data object if available.
        Otherwise they are generated using the InChI strings.

        :param data: information from input species.csv file
        :type data: pandas
        :param idxs: index for the dict to be created
        :type idxs: list(str)
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: SMILES]
    """

    spc_dct = {}
    if hasattr(data, 'smiles'):
        spc_dct = dict(zip(idxs, data.smiles))
    elif hasattr(data, 'inchi'):
        smiles = [_smiles(ich) for ich in data.inchi]
        spc_dct = dict(zip(idxs, smiles))
    else:
        spc_dct = {}
        print('No "SMILES" or "InChI" column in csv file')

    return spc_dct


def _read_csv_mult(data, idxs):
    """ Build the species dictionary relating ChemKin name to multiplicity.

        :param data: information from input species.csv file
        :type data: pandas
        :param idxs: index for the dict to be created
        :type idxs: list(str)
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: multiplicity]
    """

    if hasattr(data, 'mult'):
        spc_dct = dict(zip(idxs, data.mult))
    else:
        spc_dct = {}
        print('No "mult" column in csv file')

    return spc_dct


def _read_csv_charge(data, idxs):
    """ Build the species dictionary relating ChemKin name to charge.
        If the charge is missing for a given species, the dictionary
        element will be set to zero, assuming a neutral species.

        :param data: information from input species.csv file
        :type data: pandas
        :param idxs: index for the dict to be created
        :type idxs: list(str)
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: charge]
    """

    fill = 0
    if hasattr(data, 'charge'):
        spc_dct = dict(zip(idxs, data.charge))
    else:
        if fill is not None:
            spc_dct = dict(zip(idxs, [fill for name in idxs]))
        else:
            spc_dct = {}

    return spc_dct


def _read_csv_sensitivity(data, idxs):
    """ Build the species dictionary relating ChemKin name to sensitivity.
        If the sensitivity is missing for a given species, the dictionary
        element will be set to zero.

        :param data: information from input species.csv file
        :type data: pandas
        :param idxs: index for the dict to be created
        :type idxs: list(str)
        :return spc_dct: output dictionary for all species
        :rtype spc_dct: dict[name: sensitivity]
    """

    fill = 0.
    if hasattr(data, 'sens'):
        spc_dct = dict(zip(idxs, data.sens))
    else:
        if fill is not None:
            spc_dct = dict(zip(idxs, [fill for name in idxs]))
        else:
            spc_dct = {}

    return spc_dct


def _read_csv_headers(data):
    """ Get the names of the CSV headers
    """
    return list(data.head)


def _check_csv(data):
    """ Check if the csv is formatted properly
    """

    headers = set(list(data.head()))

    req1 = {'name', 'smiles', 'mult'}
    req2 = {'name', 'inchi', 'mult'}
    if req1 <= headers or req2 <= headers:
        proper = True
    else:
        proper = False
        print('Required Headers (name, smiles/inchi/mult) missing.')

    if headers <= set(ALLOWED_HEADERS):
        proper = True
    else:
        print('Unallowed Headers in Header file.')
        proper = False

    # Check for repeat names
    if len(list(data.name)) > len(set(list(data.name))):
        proper = False

    # Check validity of inchi and multiplicity combinations (and chg?)
    # assert _is_valid_inchi_multiplicity(ich, mul)

    # Add check to see for equivalence of columns
    # if proper:
    #     num_rows = len(data.name)
    #     print(num_rows)
    #     print(data)
    #     for col in data.columns:
    #         print(col)
    #         print(len(col))
    #         if len(col) != num_rows:
    #             proper = False
    # print('proper3', proper)

    return proper


# Set dct for above readers
READERS = {
    'smiles': _read_csv_smiles,
    'inchi': _read_csv_inchi,
    'inchikey': _read_csv_inchikey,
    'mult': _read_csv_mult,
    'charge': _read_csv_charge,
    'sens': _read_csv_sensitivity
}


def read_csv_headers(csv_str):
    """ Get the names of the CSV headers
    """
    data = _read_csv(csv_str)
    return list(data.head())


def _read_csv(csv_str):
    """ Read the csv file and generate data using pandas.

        :param csv_str: string of input csv file with species information
        :type csv_str: str
        :return data: data for the species from csv file
        :rtype: pandas.data object?
    """

    # Read in csv file while removing whitespace and make all chars lowercase
    csv_file = StringIO(csv_str)
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")

    # Parse CSV string into data columns
    data.columns = data.columns.str.strip()
    data.columns = map(str.lower, data.columns)

    return data
