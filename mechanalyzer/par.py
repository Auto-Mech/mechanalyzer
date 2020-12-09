"""
  Set various type parameters used in automol
"""


class SPC():
    """ Params for species
    """
    NAME = 'name'
    INCHI = 'inchi'
    INCHI_KEY = 'inchikey'
    SMILES = 'smiles'
    CHARGE = 'charge'
    MULT = 'mult'
    # MULT = 'multiplicity'
    SENS = 'sensitivity'
    ELECLVLS = 'electronic_levels'


class RXN():
    """ Params for reactions
    """
    REACS = 'reacs'
    PRODS = 'prods'


class TS():
    """ Params for TSs
    """
    FRM_BND_KEYS = 'frm_bnd_keys'
    BRK_BND_KEYS = 'brk_bnd_keys'


class THY():
    """ Params for thy
    """
    PROGRAM = 'program'
    METHOD = 'method'
    BASIS = 'basis'
    ORB_RESTRICT = 'orb_restrict'  # fix
