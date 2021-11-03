""" species dict function tests
"""

import os
import ioformat
import mechanalyzer


# Set paths
PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

# Set information for running tests
SMILES = ('CC', 'CC(O)Cl')

HEADERS = ('smiles', 'inchi', 'mult', 'charge')


def test__csv_io():
    """ test mechanalyzer.parser.spc.dct
        test mechanalyzer.parser.spc.dct_to_csv_str
    """

    ref_csv_str = ioformat.pathtools.read_file(DAT_PATH, 'spc3.csv')

    spc_dct = mechanalyzer.parser.spc.build_spc_dct(ref_csv_str, 'csv')
    csv_str = mechanalyzer.parser.spc.csv_string(spc_dct, HEADERS)

    assert ref_csv_str == csv_str


# builders
def test__spc_dct_build():
    """ test mechanalyzer.parser.spc.spc_dct_from_smiles
    """

    ref_spc_dct = {
        'C2H6': {
            'smiles': 'CC',
            'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3',
            'inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'C2H5ClO': {
            'smiles': 'C[C@@H](O)Cl',
            'inchi': 'InChI=1S/C2H5ClO/c1-2(3)4/h2,4H,1H3/t2-/m1/s1',
            'inchikey': 'KJESGYZFVCIMDE-UWTATZPHSA-N',
            'charge': 0, 'mult': 1,
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)}
    }

    spc_dct = mechanalyzer.parser.spc.spc_dct_from_smiles(
        SMILES, stereo=True)
    assert set(ref_spc_dct.keys()) == set(spc_dct.keys())
    for name, ref_dct in ref_spc_dct.items():
        assert set(ref_dct.keys()) == set(spc_dct[name].keys())


# modify/add functionality
def test__mod_spc_dct_atomcount():
    """ test mechanalyzer.parser.spc.reorder_by_atomcount
    """

    ref_spc_dct = {
        'O2': {'inchi': 'InChI=1S/O2/c1-2'},
        'H2': {'inchi': 'InChI=1S/H2/h1H'},
        'OH': {'inchi': 'InChI=1S/HO/h1H'},
        'CH4': {'inchi': 'InChI=1S/CH4/h1H4'},
        'CH3OH': {'inchi': 'InChI=1S/CH4O/c1-2/h2H,1H3'},
        'CH3SH': {'inchi': 'InChI=1S/CH4S/c1-2/h2H,1H3'},
        'C2H6': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
        'C2H5NH2': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
        'C2H5OH': {'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'},
        'C3H8': {'inchi': 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'}
    }

    spc_dct = {
        'C2H6': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
        'O2': {'inchi': 'InChI=1S/O2/c1-2'},
        'C3H8': {'inchi': 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'},
        'CH3OH': {'inchi': 'InChI=1S/CH4O/c1-2/h2H,1H3'},
        'C2H5OH': {'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'},
        'CH4': {'inchi': 'InChI=1S/CH4/h1H4'},
        'H2': {'inchi': 'InChI=1S/H2/h1H'},
        'CH3SH': {'inchi': 'InChI=1S/CH4S/c1-2/h2H,1H3'},
        'OH': {'inchi': 'InChI=1S/HO/h1H'},
        'C2H5NH2': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
    }

    spc_dct = mechanalyzer.parser.spc.reorder_by_atomcount(ref_spc_dct)
    assert ref_spc_dct == spc_dct


def test__mod_spc_dct_hof_basis():
    """ test mechanalyzer.parser.spc.add_heat_of_formation_basis
    """

    spc_dct = {
        'C4H9OH': {'inchi': 'InChI=1S/C4H10O/c1-2-3-4-5/h5H,2-4H2,1H3'},
        'CH4': {'inchi': 'InChI=1S/CH4/h1H4'},
        'H2O': {'inchi': 'InChI=1S/H2O/h1H2'}
    }

    ref_spc_dct = {
        'C4H9OH': {
            'inchi': 'InChI=1S/C4H10O/c1-2-3-4-5/h5H,2-4H2,1H3'},
        'CH4': {
            'inchi': 'InChI=1S/CH4/h1H4'},
        'H2O': {
            'inchi': 'InChI=1S/H2O/h1H2'},
        'cbh0_[HH]': {
            'smiles': '[HH]',
            'inchi': 'InChI=1S/H2/h1H',
            'inchikey': 'UFHFLCQGNIYNRP-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'cbh1_CC': {
            'smiles': 'CC',
            'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3',
            'inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'cbh1_CO': {
            'smiles': 'CO',
            'inchi': 'InChI=1S/CH4O/c1-2/h2H,1H3',
            'inchikey': 'OKKJLVBELUTLKV-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'cbh2_CCC': {
            'smiles': 'CCC',
            'inchi': 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3',
            'inchikey': 'ATUOYWHBWRKTHZ-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.5235987755982988,
            'hbond_cutoffs': (4.55, 1.92)},
        'cbh2_CCO': {
            'smiles': 'CCO',
            'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
            'inchikey': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
            'charge': 0, 'mult': 1,
            'mc_nsamp': (True, 3, 1, 3, 100, 12),
            'hind_inc': 0.523598775598298,
            'hbond_cutoffs': (4.55, 1.92)}
    }

    spc_dct = mechanalyzer.parser.spc.add_heat_of_formation_basis(
        spc_dct, ref_schemes=('cbh0', 'cbh1', 'cbh2'))
    spc_dct2 = mechanalyzer.parser.spc.add_heat_of_formation_basis(
        spc_dct, ref_schemes=('cbh0', 'cbh1', 'cbh2'))

    assert set(ref_spc_dct.keys()) == set(spc_dct.keys())
    assert set(ref_spc_dct.keys()) == set(spc_dct2.keys())
    for name, ref_dct in ref_spc_dct.items():
        assert set(ref_dct.keys()) == set(spc_dct[name].keys())
        assert set(ref_dct.keys()) == set(spc_dct2[name].keys())


def test__mod_spc_dct_stereo():
    """ test mechanalyzer.parser.spc.stereochemical_spc_dct
    """

    spc_dct = {
        'CC(O)Cl': {'inchi': 'InChI=1S/C2H5ClO/c1-2(3)4/h2,4H,1H3'},
        'CC=CC': {'inchi': 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+'},
        'OC=CN': {'inchi': 'InChI=1S/C2H5NO/c3-1-2-4/h1-2,4H,3H2'},
        'CC': {'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3'},
    }

    ref_spc_dct = {
        'CC(O)Cl': {
            'inchi': 'InChI=1S/C2H5ClO/c1-2(3)4/h2,4H,1H3/t2-/m0/s1',
            'inchikey': 'KJESGYZFVCIMDE-REOHCLBHSA-N'
        },
        'CC=CC': {
            'inchi': 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+',
            'inchikey': 'IAQRGUVFOMOMEM-ONEGZZNKSA-N'
        },
        'OC=CN': {
            'inchi': 'InChI=1S/C2H5NO/c3-1-2-4/h1-2,4H,3H2/b2-1-',
            'inchikey': 'UEVZOFFYIPNZNW-UPHRSURJSA-N'
        },
        'CC': {
            'inchi': 'InChI=1S/C2H6/c1-2/h1-2H3',
            'inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N'
        }
    }

    spc_dct = mechanalyzer.parser.spc.stereochemical_spc_dct(
        ref_spc_dct, all_stereo=False)

    assert ref_spc_dct == spc_dct
