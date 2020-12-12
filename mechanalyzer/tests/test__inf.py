""" test chemkin_io.calculator.thermo
"""

from mechanalyzer import par
from mechanalyzer import inf


INCHI1 = 'InChI=1S/C2H4/c1-2/h1-2H2'
CHARGE1 = 0
MULT1 = 1

INCHI2 = 'InChI=1S/H'
CHARGE2 = 0
MULT2 = 2

INCHI3 = 'InChI=1S/C2H5/c1-2/h1H2,2H3'
CHARGE3 = 0
MULT3 = 2

SPC_DCT = {
    'C2H4': {
        par.SPC.INCHI: INCHI1,
        par.SPC.CHARGE: CHARGE1,
        par.SPC.MULT: MULT1
    },
    'H': {
        par.SPC.INCHI: INCHI2,
        par.SPC.CHARGE: CHARGE2,
        par.SPC.MULT: MULT2
    },
    'C2H5': {
        par.SPC.INCHI: INCHI3,
        par.SPC.CHARGE: CHARGE3,
        par.SPC.MULT: MULT3
    },
}


PROG1 = 'gaussian09'
METHOD1 = 'wb97xd'
BASIS1 = '6-31G*'
ORB_RESTRICT1 = 'RU'

PROG2 = 'molpro2025'
METHOD2 = 'ccsd(t)'
BASIS2 = 'cc-pvtz'
ORB_RESTRICT2 = 'RR'

THY_DCT = {
    par.THY.PROGRAM: PROG2,
    par.THY.METHOD: METHOD2,
    par.THY.BASIS: BASIS2,
    par.THY.ORB_RESTRICT: ORB_RESTRICT2
}

REACS = ['C2H4', 'H']
PRODS = ['C2H5']


def test__inf():
    """ test mechanalyzer.inf
    """

    # SPC INFO OBJECTS

    # Build spc_info objects
    spc_info1 = inf.spc.from_data(INCHI1, CHARGE1, MULT1)
    spc_info2 = inf.spc.from_dct(SPC_DCT['H'])
    spc_info3 = inf.spc.combine(spc_info1, spc_info2, mval='max')

    assert spc_info1 == (INCHI1, CHARGE1, MULT1)
    assert spc_info2 == (INCHI2, CHARGE2, MULT2)
    assert spc_info3 == ((INCHI1, INCHI2), 0, 2)

    # Get a value
    ich = inf.spc.value(spc_info1, par.SPC.INCHI)

    assert ich == INCHI1

    # THEORY INFO OBJECTS

    # Build thy_info objects
    thy_info1 = inf.thy.from_data(PROG1, METHOD1, BASIS1, ORB_RESTRICT1)
    thy_info2 = inf.thy.from_dct(THY_DCT)

    assert thy_info1 == (PROG1, METHOD1, BASIS1, ORB_RESTRICT1)
    assert thy_info2 == (PROG2, METHOD2, BASIS2, ORB_RESTRICT2)

    # Modify the orb label of the thy_info object
    mod_thy_info1 = inf.thy.modify_orb_label(thy_info1, spc_info2)

    assert mod_thy_info1 == (PROG1, METHOD1, BASIS1, 'U')

    # Get a value
    basis = inf.thy.value(thy_info1, par.THY.BASIS)

    assert basis == BASIS1

    # REACTION INFO OBJECTS
    rxn_info = inf.rxn.from_dct(REACS, PRODS, SPC_DCT, rxn_mul='low')
    assert rxn_info == (
        (('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/H'),
         ('InChI=1S/C2H5/c1-2/h1H2,2H3',)),
        ((0, 0),
         (0,)),
        ((1, 2), (2,)),
        2
    )

    sort_rxn_info = inf.rxn.sort(rxn_info)
    assert sort_rxn_info == (
        (('InChI=1S/C2H5/c1-2/h1H2,2H3',),
         ('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/H')),
        ((0,),
         (0, 0)),
        ((2,),
         (1, 2)),
        2
    )

    rxn_ichs = inf.rxn.value(rxn_info, par.SPC.INCHI)
    assert rxn_ichs == (
        ('InChI=1S/C2H4/c1-2/h1-2H2', 'InChI=1S/H'),
        ('InChI=1S/C2H5/c1-2/h1H2,2H3',)
    )


if __name__ == '__main__':
    test__inf()
