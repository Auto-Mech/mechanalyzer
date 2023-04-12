""" Tests the strip_ste module
"""

import os
import tempfile
import numpy
from mechanalyzer.parser import new_spc as spc_parser
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.builder import strip_ste
from mechanalyzer.calculator.rates import check_p_t
from chemkin_io.writer import mechanism
from ioformat import pathtools

# Set Paths to test/data directory and output directory
DAT_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
TMP_OUT = tempfile.mkdtemp()

# Filenames
SPC_CSV = 'strip_ste.csv'
CKIN = 'strip_ste.ckin'

# Define temps and pressures (for refitting rates)
TEMPS_LST = [numpy.linspace(300, 2500, 12)]
PRESSURES = [1, 10]
TEMPS_LST = check_p_t(TEMPS_LST, PRESSURES)  # checks temps list and reforms

# Load things
MECH_SPC_DCT = spc_parser.load_mech_spc_dct(
    SPC_CSV, DAT_PATH, chk_ste=False, chk_match=False)
RXN_PARAM_DCT = ckin_parser.load_rxn_param_dct(CKIN, DAT_PATH)


def test_renaming():
    """ Tests the renaming of the stereo species since this depends on the
        name generation done in builder/_names
    """

    mech_spc_dct_strpd_ich = {
        'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3': {
            'smiles': 'CC=CC', 'inchi': 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3',
            'inchikey': 'IAQRGUVFOMOMEM-UHFFFAOYSA-N', 'mult': 1,
            'charge': 0, 'exc_flag': 0, 'fml': {'C': 4, 'H': 8}},
        'InChI=1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3': {
            'smiles': 'CC1OC1C',
            'inchi': 'InChI=1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3',
            'inchikey': 'PQXKWPLDPFFDJP-UHFFFAOYSA-N',
            'mult': 1, 'charge': 0, 'exc_flag': 0,
            'fml': {'C': 4, 'H': 8, 'O': 1}}
    }

    re_mech_spc_dct, _ = strip_ste.regenerate_names(mech_spc_dct_strpd_ich, {})
    assert set(list(re_mech_spc_dct.keys())) == set(
        ['C4H8OCETH-RvEsWv', 'C4H8ALK-808sWv']), ('Something has changed about the'
        f' name from inchi generator! {re_mech_spc_dct.keys()}')


def test_main():
    """ Tests the main strip_ste function
    """

    re_rxn_param_dct_comb, re_mech_spc_dct_comb = strip_ste.main(
        RXN_PARAM_DCT, MECH_SPC_DCT, TEMPS_LST, PRESSURES)

    # Check values
    # Reaction 1
    ref_vals1 = [5e12, 0, 0]
    vals1 = re_rxn_param_dct_comb[
        (('C4H8ALK-808sWv', 'HO2'), ('C4H8OCETH-RvEsWv', 'OH'), (None,))].arr
    assert numpy.allclose(ref_vals1, vals1)
    # Reaction 2
    ref_vals2 = [3e12, 0, 0]
    vals2 = re_rxn_param_dct_comb[
        (('C4H8ALK-808sWv', 'HO2'), ('H', 'OH'), (None,))].arr
    assert numpy.allclose(ref_vals2, vals2)
    # Reaction 3
    ref_vals3 = [3.5e12, 0, 0]
    vals3 = re_rxn_param_dct_comb[
        (('H', 'HO2'), ('C4H8OCETH-RvEsWv', 'OH'), (None,))].arr
    assert numpy.allclose(ref_vals3, vals3)
    # Reaction 4
    ref_vals4 = [1e12, 0, 0]
    vals4 = re_rxn_param_dct_comb[
        (('O', 'O'), ('O2',), (None,))].arr
    assert numpy.allclose(ref_vals4, vals4)
    # Reaction 5
    ref_vals5 = [1e12, 0, 0]
    vals5 = re_rxn_param_dct_comb[
        (('H', 'OH'), ('H2O',), (None,))].arr
    assert numpy.allclose(ref_vals5, vals5)

    # Just writing this in case someone wants to go look at it
    mech_str = mechanism.write_chemkin_file(
        rxn_param_dct=re_rxn_param_dct_comb,
        mech_spc_dct=re_mech_spc_dct_comb)
    pathtools.write_file(mech_str, DAT_PATH, 'strip_ste.out')


if __name__ == '__main__':
    test_renaming()
    test_main()
