""" BASIC SCHEME (JUST MOVE TO SPECIES?)
"""

import numpy
import automol.inchi
import automol.geom


# Basic Basis List builders
def basic_ts_basis(zrxn):
    """ Determine a basis for relative enthalpy calculations for a
        transition states by looping over the reactants and building
        a list of a simple species in a inear combination which
        reproduces the number of atoms in the stoichiometry of the species.

        Basis represented by list of InChI strings and array of coefficients.

        :param zrxn: reaction object oriented to Z-Matrix
        :type zrxn: automol.reac.Reaction object
        :rtype: (tuple(str), numpy.ndarray)
    """

    # Just use reactants
    rxn_ichs = automol.reac.reaction_inchis(zrxn)
    rct_ichs, _ = rxn_ichs

    basis, coeff_lst = [], []
    for ich in rct_ichs:
        spc_bas_i, coeff_bas_i = basic_spc_basis(ich)
        for bas_i, c_bas_i in zip(spc_bas_i, coeff_bas_i):
            if bas_i not in basis:
                # Add basis and coefficients to list
                basis.append(bas_i)
                coeff_lst.append(c_bas_i)
            else:
                # Add coefficient value to existing coefficient value
                for j, bas_j in enumerate(basis):
                    if bas_i == bas_j:
                        coeff_lst[j] += c_bas_i

    return (tuple(basis), numpy.array(coeff_lst))


def basic_spc_basis(ich):
    """ Determine a basis for relative enthalpy calculations for species
        by building a list of a simple species in a inear combination which
        reproduces the number of atoms in the stoichiometry of the species.

        Basis represented by list of InChI strings and array of coefficients.

        :param ich: InChI string for spc
        :type ich: str
        :rtype: (tuple(str), numpy.ndarray)
    """

    # Get a list of all the atom types in the molecule
    fml = automol.inchi.formula(ich)
    symbs = tuple(fml.keys())

    # Create list of inchi keys corresponding to basis species
    basis = ()
    # H2
    basis += ('InChI=1S/H2/h1H',)
    # CH4
    if 'C' in symbs:
        basis += ('InChI=1S/CH4/h1H4',)
    # H2O
    if 'O' in symbs:
        basis += ('InChI=1S/H2O/h1H2',)
    # NH3
    if 'N' in symbs:
        basis += ('InChI=1S/H3N/h1H3',)
    # Cl2
    if 'Cl' in symbs:
        basis += ('InChI=1S/ClH/h1H',)
        # basis += ('InChI=1S/Cl2/c1-2',)
    # SO2
    if 'S' in symbs:
        basis += ('InChI=1S/O2S/c1-3-2',)
        if 'O' not in symbs:
            basis += ('InChI=1S/H2O/h1H2',)

    # Generate the coefficients for the basis
    if len(basis) == 1 and ich == basis[0]:
        coeff_lst = (1.0,)
    else:
        coeff_lst = _coefficients(basis, fml)

    return basis, coeff_lst


def _coefficients(basis, spc_fml):
    """ Form a matrix consisting of the coefficients for the basis
        species to balance out the atoms.

        :param basis: InChI strings for basis species
        :type basis: tuple(str)
        :param spc_fml: stoichiometric formula of species basis corresponds to
        :type spc_fml: dict[str:int]
        :rtype: numpy.ndarray
    """

    # Initialize an natoms x natoms matrix
    nbasis = len(basis)
    basis_mat = numpy.zeros((nbasis, nbasis))

    # Get the basis formulae list
    basis_fml_str = [automol.inchi.formula_string(spc) for spc in basis]
    for spc in basis_fml_str:
        basis_atom_dict = automol.formula.from_string(spc)
        for atom in basis_atom_dict:
            if atom not in spc_fml:
                spc_fml[atom] = 0

    # Set the elements of the matrix
    for i, spc in enumerate(basis_fml_str):
        basis_atom_dict = automol.formula.from_string(spc)
        basis_vals = []
        for key in spc_fml.keys():
            if key in basis_atom_dict:
                basis_vals.append(basis_atom_dict[key])
            else:
                basis_vals.append(0)
        basis_mat[i] = basis_vals

    #  Transpose
    basis_mat = basis_mat.T

    # Form stoichometry vector vector
    stoich_vec = numpy.zeros(len(spc_fml))
    for i, key in enumerate(spc_fml.keys()):
        stoich_vec[i] = spc_fml[key]

    # Solve C = B^-1 S
    basis_mat = numpy.linalg.inv(basis_mat)
    coeff_vec = numpy.dot(basis_mat, stoich_vec)

    return coeff_vec


# unused basis reduction code
# def get_reduced_basis(basis_ich, species_formula):
#     """
#     Form a matrix for a given basis and atomlist
#     INPUT:
#     input_basis     - ich strings for set of reference molecules
#     atomlist  - list of atoms (all atoms that appear
#                 in basis should be in atomlist)
#     OUTPUT:
#     mat       - matrix (length of basis by length of atomlist)
#                 (square if done right)
#     """
#
#     # Get the basis formulae list
#     basis_formulae = [automol.inchi.formula_string(spc) for spc in basis_ich]
#
#     reduced_basis = []
#     for i, basis_formula in enumerate(basis_formulae):
#         basis_atom_dict = automol.formula.from_string(basis_formula)
#         flag = True
#         for key, _ in basis_atom_dict.items():
#             if key not in species_formula:
#                 flag = False
#
#         if flag:
#             reduced_basis.append(basis_ich[i])
#
#     return reduced_basis
