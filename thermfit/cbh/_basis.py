# Basis builders
def get_basis(ich):
    """ get a basis
    """

    formula = automol.inchi.formula_string(ich)
    atm_dict = automol.formula.from_string(formula)
    return select_basis(atm_dict)


def get_basic(ich):
    """ get basis for basic scheme
    """
    formula_dct = automol.inchi.formula(ich)
    spc_bas = select_basis(formula_dct)
    if len(spc_bas) == 1 and ich == spc_bas[0]:
        clist = [1]
    else:
        clist = calc_coefficients(spc_bas, formula_dct)
    return spc_bas, clist
