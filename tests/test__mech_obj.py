"""
Test the handling of objects
"""

from mechanalyzer import par
from mechanalyzer import spc
from mechanalyzer import thy


# DATA FOR SPECIES
INCHI = 'InChI=1S/H2O/h1H2'
CHARGE = 0
MULT = 1

# SPECIES DCT
SPC_DCT = {
    par.SPC.INCHI: INCHI,
    par.SPC.CHARGE: CHARGE,
    par.SPC.MULT: MULT
}


# DATA FOR THEORY
PROGRAM = 'psi4'
METHOD = 'ccsd(t)'
BASIS = 'cc-pvtz'
ORB_LABEL = 'RU'
MOD_ORB_LABEL = 'R'


def test__spc_info():
    """ test mechanalyzer._obj.par
        test mechanalyzer._obj.spc
    """

    # Build reference spc_info object
    ref_spc_info = (INCHI, CHARGE, MULT)

    # Build the spc info objects from sources
    spc_info1 = spc.from_data(INCHI, CHARGE, MULT)
    spc_info2 = spc.from_dct(SPC_DCT)
    assert spc_info1 == ref_spc_info
    assert spc_info2 == ref_spc_info

    # Read info from the objects
    assert spc.value(ref_spc_info, par.SPC.INCHI) == INCHI
    assert spc.value(ref_spc_info, par.SPC.CHARGE) == CHARGE
    assert spc.value(ref_spc_info, par.SPC.MULT) == MULT


def test__thy_info():
    """ test mechanalyzer._obj.par
        test mechanalyzer._obj.spc
    """

    # Build reference spc_info object
    ref_thy_info = (PROGRAM, METHOD, BASIS, ORB_LABEL)
    ref_mod_thy_info = (PROGRAM, METHOD, BASIS, MOD_ORB_LABEL)

    # Build the thy info objects from sources
    thy_info = thy.from_data(PROGRAM, METHOD, BASIS, ORB_LABEL)
    mod_thy_info = thy.from_data(PROGRAM, METHOD, BASIS, MOD_ORB_LABEL)
    assert thy_info == ref_thy_info
    assert mod_thy_info == ref_mod_thy_info

    # Read info from the objects
    assert thy.value(ref_thy_info, par.THY.PROGRAM) == PROGRAM
    assert thy.value(ref_thy_info, par.THY.METHOD) == METHOD
    assert thy.value(ref_thy_info, par.THY.BASIS) == BASIS
    assert thy.value(ref_thy_info, par.THY.ORB_LABEL) == ORB_LABEL

    assert thy.value(ref_mod_thy_info, par.THY.PROGRAM) == PROGRAM
    assert thy.value(ref_mod_thy_info, par.THY.METHOD) == METHOD
    assert thy.value(ref_mod_thy_info, par.THY.BASIS) == BASIS
    assert thy.value(ref_mod_thy_info, par.THY.ORB_LABEL) == MOD_ORB_LABEL


if __name__ == '__main__':
    test__spc_info()
    test__thy_info()
