
# basis.py
def get_reduced_basis(basis_ich, species_formula):
    """
    Form a matrix for a given basis and atomlist
    INPUT:
    input_basis     - ich strings for set of reference molecules
    atomlist  - list of atoms (all atoms that appear
                in basis should be in atomlist)
    OUTPUT:
    mat       - matrix (length of basis by length of atomlist)
                (square if done right)
    """

    # Get the basis formulae list
    basis_formulae = [automol.inchi.formula_string(spc) for spc in basis_ich]

    reduced_basis = []
    for i, basis_formula in enumerate(basis_formulae):
        basis_atom_dict = automol.formula.from_string(basis_formula)
        flag = True
        for key, _ in basis_atom_dict.items():
            if key not in species_formula:
                flag = False

        if flag:
            reduced_basis.append(basis_ich[i])

    return reduced_basis
