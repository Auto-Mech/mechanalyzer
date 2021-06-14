""" Computes the Heat of Formation at 0 K for a given species
"""

import os
import csv
from phydat import phycon


# Path  the database files (stored in the thermo src directory)
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def calc_hform_0k(spc_h0, basis_h0,
                  basis_ichs, basis_coeffs, ref_set, rxn=False):
    """ Calculates the change in the heat-of-formation for a species or
        transition state at 0 K. Returns value in Hartrees.

        Add TEMP VALUE FOR GENERALIZATION

        :param spc_h0: 0 K Heat-of-Formation for species or transition state
        :type spc_h0: float
        :param basis_ichs: InChI strings for basis species
        :type basis_ichs: tuple(str)
        :param basis_coeffs: basis-set coefficients corresponding to basis
        :type basis_coeffs: tuple(float)
        :param ref_set: database set to read value from
        :type ref_set: str
        :param rxn: parameter to calculate value for species or reaction
        :type rxn: bool
    """

    dhzero = spc_h0
    for i, ich in enumerate(basis_ichs):

        # Read reference enthalpies for basis molecule
        ref_h0 = reference_enthalpy(ich, ref_set, 0, rxn=rxn)
        ref_h0 = ref_h0 if ref_h0 is not None else 0.0  # break loop?

        # Add basis and reference energies to overall va
        dhzero += basis_coeffs[i] * (ref_h0 - basis_h0[i])

    return dhzero


def reference_enthalpy(ich_lookup, ref_set, temp, rxn=False):
    """ Reads a reference enthalpy for a species or transition state
        from either the ANL0 or ATcT data set. Value returned in Hartrees.

        TS InChI format: rct1_ich+rct2_ich=prd1_ich+prd2_ich

        :param ich_lookup: InChI string(s) of species/TS to lookup
        :type ich_lookup: str
        :param ref_set: database set to read value from
        :type ref_set: str
        :param temp: temperature to obtain enthalpy value for
        :type temp: int
        :param rxn: parameter to read value for species or reaction
        :type rxn: bool
    """

    # Set path and name to thermo database file
    thermodb_file = _thermo_database(temp, rxn=rxn)

    # Find the energy value for the given species and enery type
    hf_val = None
    with open(thermodb_file, 'r') as db_file:
        reader = csv.DictReader(db_file)
        for row in reader:
            if row['inchi'] == ich_lookup:
                val = row[ref_set]
                if val == '':
                    val = None
                hf_val = float(val)

    # Convert units if val found, else print error message
    if hf_val is not None:
        if not rxn:
            hf_val *= phycon.KCAL2EH
        else:
            hf_val *= phycon.KJ2EH
    else:
        print('No Heat for Formation exists:')
        print('SPC:{} Set:{} Temp:{}K'.format(ich_lookup, ref_set, temp))

    return hf_val


def format_reaction_inchi(rxn_ichs):
    """ Format the InChI strings of the reaction into a string to
        look up a reaction enthalpy in the database.

        :param rxn_ichs: InChI strings of reactants and products
        :type rxn_ichs: tuple(tuple(str), tuple(str))
        :rtype: str
    """

    rct_ichs, prd_ichs = rxn_ichs
    rct_str, prd_str = '+'.join(rct_ichs), '+'.join(prd_ichs)
    species = '='.join([rct_str, prd_str])

    return species


def _thermo_database(temp, rxn=False):
    """ Set the path to appropriate database to read based on the
        desired temperature and whether the user intends to read a value
        for a species or transition state.

        :param temp: temperature to obtain enthalpy value for
        :type temp: int
        :param rxn: parameter to look for species or reaction
        :type rxn: bool
        :rtype: str
    """

    if rxn:
        thermodb_name = 'tsthermodb_{}K.csv'.format(str(int(temp)))
    else:
        thermodb_name = 'Hfdb_{}K.csv'.format(str(int(temp)))

    return os.path.join(SRC_PATH, 'thermdb', thermodb_name)
