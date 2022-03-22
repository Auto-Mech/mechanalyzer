"""
Parses a spc.csv file containing species information to obtain a mech_spc_dct
"""

import csv
import copy
from ioformat import pathtools
from automol.smiles import inchi as smi_to_ich
from automol.inchi import smiles as ich_to_smi
from automol.inchi import formula as ich_to_fml
from automol.inchi import low_spin_multiplicity as _low_spin_mult
from automol.inchi import recalculate  # don't think I need this
from automol.inchi import is_complete
from automol.formula import from_string as str_to_fml
import rdkit.Chem as _rd_chem

ALLOWED_COLUMN_NAMES = (
    'name',
    'smiles',
    'inchi',
    'inchikey',
    'mult',
    'charge',
    'exc_flag',
    'fml',
    'sens'
)

TRIP_DCT = {
    'O': 'InChI=1S/O',
    'O2': 'InChI=1S/O2/c1-2',
    'CH2': 'InChI=1S/CH2/h1H2',
    'C3H2': 'InChI=1S/C3H2/c1-3-2/h1-2H',
    'H2CCC': 'InChI=1S/C3H2/c1-3-2/h1H2',
    'C6H6O-B1': 'InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-2,4-5H,3H2',
    'CH2CHN': 'InChI=1S/C2H3N/c1-2-3/h2H,1H2',
}


def load_mech_spc_dcts(filenames, path, quotechar="'", chk_ste=True,
                       chk_match=True):
    """ Obtains multiple mech_spc_dcts given a list of spc.csv filenames

        :param filenames: filenames of the spc.csv file to be read
        :type filename: list
        :param quotechar: character used to optionally ignore commas; " or '
        :type quotechar: str
        :return mech_spc_dcts: list of mech_spc_dcts
        :rtype: list
    """

    mech_spc_dcts = []
    for filename in filenames:
        print(f'Loading mech_spc_dct for the file {filename}...')
        mech_spc_dct = load_mech_spc_dct(filename, path,
                                             quotechar=quotechar,
                                             chk_ste=chk_ste,
                                             chk_match=chk_match)
        mech_spc_dcts.append(mech_spc_dct)

    return mech_spc_dcts


def load_mech_spc_dct(filename, path, quotechar="'", chk_ste=True,
                          chk_match=True):
    """ Obtains a single mech_spc_dct given a spc.csv filename

        :param filename: filename of the spc.csv file to be read
        :type filename: str
        :param quotechar: character used to optionally ignore commas; " or '
        :type quotechar: str
        :return mech_spc_dct: identifying information on species in a mech
        :rtype: dct {spc1: spc_dct1, spc2: ...}
    """

    file_str = pathtools.read_file(path, filename, remove_comments='!')
    mech_spc_dct = parse_mech_spc_dct(file_str, quotechar=quotechar,
                                    chk_ste=chk_ste, chk_match=chk_match)

    return mech_spc_dct


def parse_mech_spc_dct(file_str, quotechar="'", chk_ste=True, chk_match=True):
    """ Obtains a single mech_spc_dct given a string parsed from a spc.csv file

        :param file_str: the string that was read directly from the .csv file
        :type file_str: str
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return mech_spc_dct: identifying information on species in a mech
        :rtype: dct {spc1: spc_dct1, spc2: ...}
    """

    # Check for incorrect quotechar usage
    if quotechar == '"':
        assert "'" not in file_str, (
            f'A quotechar input of double quotation marks was selected, but at '
            f'least one instance of a single quotation mark exists in the '
            f'csv file. Use only double quotation marks.')
    elif quotechar == "'":
        assert '"' not in file_str, (
            f'A quotechar input of single quotation marks was selected, but at '
            f'least one instance of a double quotation mark exists in the '
            f'csv file. Use only single quotation marks.')
    else:
        raise NotImplementedError(
            f'The quotechar {quotechar} is not allowed. Options are either '
            f'a single or double quotation mark.')

    # Split into lines
    lines = file_str.split('\n')
    if lines[-1] == '':
        del lines[-1]  # gets rid of the annoying last line that gets added

    # Build the mech_spc_dct line by line
    mech_spc_dct = {}
    errors = False
    for idx, line in enumerate(lines):
        if idx == 0:
            col_headers, num_cols = parse_first_line(line, quotechar=quotechar)
        else:
            cols = parse_line(line, idx, num_cols, quotechar=quotechar)
            if cols is not None:
                spc, spc_dct = read_spc_dct(cols, col_headers)
                # Check that the species name was not already defined
                assert spc not in mech_spc_dct.keys(), (
                    f'The species name {spc} appears in the csv file more than'
                    f' once. The second time is on line {idx + 1}, {line}.')
                # Fill in the spc_dct and then add it to the mech_spc_dct
                spc_dct, error = fill_spc_dct(spc_dct, spc, chk_ste=chk_ste,
                                              chk_match=chk_match)
                mech_spc_dct[spc] = spc_dct
                if error:
                    errors = True
            else:
                print(f'Line {idx + 1} appears to be empty. Skipping...')

    # Find species with the same chemical identifiers but different names
    check_for_dups(mech_spc_dct)  # prints warnings

    assert not errors, (f'Errors while parsing the .csv file! See printouts')

    return mech_spc_dct


def read_spc_dct(cols, col_headers):
    """ Reads the parsed contents of a single line and converts to a spc_dct

        :param cols: a list of the entries in the line
        :type cols: list
        :param col_headers: a list of the entries in the first line
        :type col_headers: list
        :return spc: the mechanism name for the species
        :rtype: str
        :return spc_dct: identifying information for a single species
        :rtype: dct
    """

    # Build the spc_dct for this species, column by column
    spc_dct = {}
    for idx, col in enumerate(cols):
        col_header = col_headers[idx]
        if col_header == 'name':
            spc = col
        # If charge or mult or exc_flag, convert to int if not empty
        elif col_header in ('charge', 'mult', 'exc_flag') and col != '':
            spc_dct[col_header] = int(col)
        # If on formula, convert to dict
        elif col_header == 'fml':
            spc_dct[col_header] = str_to_fml(col)
        else:
            spc_dct[col_header] = col

    return spc, spc_dct


def fill_spc_dct(spc_dct, spc, chk_ste=True, chk_match=True):
    """ Fills in missing values in a spc_dct

        :param spc_dct: identifying information for a single species
        :type: dct
        :return full_spc_dct: beefed-up spc_dct
        :rtype: dct
    """

    full_spc_dct = copy.deepcopy(spc_dct)

    # Add SMILES or InChI if missing
    if full_spc_dct['inchi'] == '' and full_spc_dct['smiles'] == '':
        print(f'For {spc}, both inchi and smiles are empty')
        error = True
    elif 'inchi' not in full_spc_dct or full_spc_dct['inchi'] == '':
        smi = full_spc_dct['smiles']
        error = check_smi(smi, spc)
        full_spc_dct['inchi'] = smi_to_ich(smi)
    elif 'smiles' not in full_spc_dct:
        ich = full_spc_dct['inchi']
        error = check_ich(ich, spc, chk_ste=chk_ste)
        full_spc_dct['smiles'] = ich_to_smi(ich)
    else:  # if smiles and inchi were both already included
        if chk_match:
            error = check_smi_and_ich(
                full_spc_dct['smiles'], full_spc_dct['inchi'], spc,
                chk_ste=chk_ste)
        else:
            error = check_ich(full_spc_dct['inchi'], spc, chk_ste=chk_ste)
            if not error:  # if the inchi passed, check the smiles
                error = check_smi(full_spc_dct['smiles'], spc)

    # Add charge and exc_flag if missing; assume 0 for both
    if 'charge' not in full_spc_dct or full_spc_dct['charge'] == '':
        full_spc_dct['charge'] = 0
    if 'exc_flag' not in full_spc_dct or full_spc_dct['exc_flag'] == '':
        full_spc_dct['exc_flag'] = 0

    # Add multiplicity if missing by guessing
    if 'mult' not in full_spc_dct or full_spc_dct['mult'] == '':
        ich = full_spc_dct['inchi']
        smi = full_spc_dct['smiles']
        mult = _low_spin_mult(ich)
        if mult == 1 and ich in TRIP_DCT.values():
            sing = False
            if '(S)' in smi or '(S)' in ich or '(S)' in spc:
                sing = True
            elif 'singlet' in smi or 'singlet' in ich or 'singlet' in spc:
                sing = True
            if not sing:  # otherwise, just leave as 1
                mult = 3
        full_spc_dct['mult'] = mult

    # Add formula if missing
    if 'fml' not in full_spc_dct or full_spc_dct['fml'] == '':
        fml = ich_to_fml(full_spc_dct['inchi'])
        full_spc_dct['fml'] = fml

    return full_spc_dct, error


def parse_first_line(first_line, quotechar="'"):
    """ Parses the first line in the spc.csv file

        :param first_line: a string with the contents of the first line in the
            csv file
        :type first_line: str
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return headers: a list of the entries in the first line
        :rtype: list
        :return num_cols: the number of columns indicated by the first line
        :rtype: int
    """

    # Remove the UTF encoding if present
    if '\ufeff' in first_line:
        first_line = first_line.replace('\ufeff', '')

    # Parse the line
    headers = next(csv.reader([first_line], quotechar=quotechar))
    num_cols = len(headers)

    # Check for appropriateness of column headers
    for col_header in headers:
        assert col_header in ALLOWED_COLUMN_NAMES, (
            f'The column header {col_header} is not allowed. Options are'
            f' {ALLOWED_COLUMN_NAMES}.')
    assert 'name' in headers, (
        "The species name, 'name', must be included in the csv file header.")
    assert 'inchi' in headers or 'smiles' in headers, (
        "At least one of the following chemical identifiers must be included in"
        " the csv file header: 'inchi' or 'smiles'.")

    return headers, num_cols


def parse_line(line, idx, num_cols, quotechar="'"):
    """ Parses a line in the spc.csv file (other than the first line)

        :param line: a string with the contents of a single line in the csv file
        :type line: str
        :param idx: the index of the line
        :type idx: int
        :param num_cols: the number of columns in the csv file
        :type num_cols: int
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return cols: a list of the entries in the line
        :rtype: list
    """

    cols = next(csv.reader([line], quotechar=quotechar))
    if line == '':
        cols = None
    else:
        assert len(cols) == num_cols, (
            f'Line {idx + 1}, {line}, has {len(cols)} columns. The first line '
            f'in the csv file indicates {num_cols} columns are needed.')

    return cols


def check_ich(ich, spc, chk_ste=True):

    error = False
    mol = _rd_chem.MolFromInchi(ich)
    if mol is None:
        print(f"The spc '{spc}' has an invalid InChI, '{ich}'")
        error = True

    # If indicated, check for stereochemical completeness of inchi
    if chk_ste:
        mol2 = is_complete(ich)
        if mol2 is None:
            print(f"The spc '{spc}' has a stereo-invalid InChI, '{ich}'")
            error = True

    return error


def check_smi(smi, spc):

    error = False
    mol = _rd_chem.MolFromSmiles(smi)
    if mol is None:
        print(f"The spc '{spc}' has an invalid SMILES string, '{smi}'")
        error = True

    return error


def check_smi_and_ich(smi, ich, spc, chk_ste=True):

    # Check the smiles and inchi for correctness
    error1 = check_smi(smi, spc)
    error2 = check_ich(ich, spc, chk_ste=chk_ste)

    # Recalculate inchi from smiles and see if it's the same
    # Note: only do if both initial tests passed and if not checking stereo
    if not any([error1, error2]) and not chk_ste:
        recalc_ich = smi_to_ich(smi)
        if ich != recalc_ich:
            print(f"The spc '{spc}' has non-matching SMILES and InChI strings")
            error3 = True
        else:
            error3 = False
    else:
        error3 = False

    error = any([error1, error2, error3])

    return error


def check_for_dups(mech_spc_dct):
    """ Checks a mech_spc_dct for species that are identical except in name
    """

    def are_spcs_same(spc_dct1, spc_dct2):
        """ Checks if two species are chemically identical
        """

        ich1 = spc_dct1['inchi']
        mlt1 = spc_dct1['mult']
        chg1 = spc_dct1['charge']
        exc1 = spc_dct1['exc_flag']

        ich2 = spc_dct2['inchi']
        mlt2 = spc_dct2['mult']
        chg2 = spc_dct2['charge']
        exc2 = spc_dct2['exc_flag']

        if ich1 == ich2 and mlt1 == mlt2 and chg1 == chg2 and exc1 == exc2:
            are_same = True
        else:
            are_same = False

        return are_same

    spcs = tuple(mech_spc_dct.keys())  # convert to tuple for indexing
    spc_dcts = tuple(mech_spc_dct.values())
    for outer_idx, outer_spc in enumerate(spcs):
        outer_dct = spc_dcts[outer_idx]
        for inner_idx in range(outer_idx + 1, len(spcs)):
            inner_spc = spcs[inner_idx]
            inner_dct = spc_dcts[inner_idx]
            if are_spcs_same(outer_dct, inner_dct):
                print(f'{outer_spc} and {inner_spc} are chemical twins!')

