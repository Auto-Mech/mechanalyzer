""" Computes the Heat of Formation at 0 K for a given species
"""

import os
import csv
from phydat import phycon


# Path  the database files (stored in the thermo src directory)
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def calc_hform_0k(hzero_mol, hzero_basis, basis, coeff, ref_set):
    """ calculates the heat-of-formation at 0 K
    """

    dhzero = hzero_mol * phycon.EH2KCAL
    # ioprinter.debug_message(
    #     'ABS of molecule in kcal: {:.5f}'.format(dhzero))
    for i, spc in enumerate(basis):
        _ts = True
        if isinstance(spc, str):
            _ts = False
        h_basis = reference_enthalpy(spc, ref_set, 0, _ts)
        if h_basis is None:
            h_basis = 0.0
            # ioprinter.warning_message(
            #     'No Heat of Formation in database for {} {} ', spc, ref_set)
        if _ts:
            h_basis *= phycon.KJ2KCAL
        dhzero += coeff[i] * h_basis
        dhzero -= coeff[i] * hzero_basis[i] * phycon.EH2KCAL
        # ioprinter.debug_message('Contribution from:', spc)
        # ioprinter.debug_message(
        #     'HF0K in kcal: {:g} * {:.5f}'.format(coeff[i], h_basis))
        # ioprinter.debug_message(
        #     'ABS in kcal: {:g} * {:.5f}'.format(
        #         coeff[i], hzero_basis[i] * phycon.EH2KCAL))

    return dhzero


def reference_enthalpy(species, ref, temp, transt=False):
    """ gets a reference value
    """

    # Set path and name to thermo database file
    thermodb_file = _thermo_database(temp, transt)

    # Find the energy value for the given species and enery type
    h_species = None
    if transt:
        rcts, prds = species
        rct_str = '+'.join(rcts)
        prd_str = '+'.join(prds)
        species = '='.join([rct_str, prd_str])
    with open(thermodb_file, 'r') as db_file:
        reader = csv.DictReader(db_file)
        for row in reader:
            if row['inchi'] == species:
                val = row[ref]
                if val == '':
                    val = None
                h_species = float(val)
    assert h_species is not None, (
        'Could not find heat of formation for {} '.format(species)
        )

    return h_species


def _thermo_database(temp, transt):
    """ Set the the name of the thermo dabase to read.
    """

    if transt:
        thermodb_name = 'tsthermodb_{}K.csv'.format(str(int(temp)))
    else:
        thermodb_name = 'Hfdb_{}K.csv'.format(str(int(temp)))
    thermodb_file = os.path.join(SRC_PATH, 'thermdb', thermodb_name)

    return thermodb_file
