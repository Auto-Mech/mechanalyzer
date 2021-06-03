tandardized mechanism objects

Basic objects
spc
description: spc_name
type: str
format: spc
rxn
description: reaction name
type: tuple
format: ((rct1, rct2, ...), (prd1, prd2, ...), (third_bod1, third_bod2, ...))
each entry (e.g., rct1) is a species (see the spc entry)
the third bodies are confusing in that they have a ‘+’ in front and may also have ‘()’ enclosing them
can be a generic third body instead of a specific species: ‘+M’ or ‘(+M)’
spc_ident_dct
description: species and their accompanying chemically unique descriptions 
type: dct
format: {spc1: ident_dct1, spc2: ...}
ident_dct
description: chemically unique description of a spc
type: dct
format: {‘smiles’: SMILES, ‘inchi’: InChI, ‘inchikey’: InChI_key, ‘mult’: multiplicity, ‘charge’: charge, ‘sens’: sensitivity, ‘fml’: fml_dct}
SMILES is a str
InChI is a str
InChI_key is a str
multiplicity is an int
charge is an int
sensitivity is a float
fml_dct is a dct describing the chemical formula of a species
for example, for formaldehyde, the fml_dct would be {‘C’: 1, ‘H’: 2, ‘O’: 1}


Objects describing rate expressions

rxn_param_dct
description: reactions and their accompanying rate expressions
type: dct
.. code_block: python
    {rxn1: (param_tup1, param_tup2, ...), rxn2: ...}

.. example:: .. code_block: python
the different param_tuples pertain to different duplicate expressions
example of double Arrhenius: {((‘H’, ‘O2’), (‘OH’, ‘O’), (None,)): ([1e15, 0, 15000], None, None, None, None, None), ([1e10, 0, 5000], None, None, None, None, None)}
example of PLOG with some (but not all) pressures with duplicate fits: {(('H', 'O2'), ('OH', 'O'), (None,)): (([1E+15, 0.00, 25000], None, None, None, {0.1: [1E+15, 0.00, 25000], 1.0: [1E+16, 0.00, 25000], 10.0: [1E+17, 0.00, 25000], 100.0: [1E+18, 0.00, 25000]}, None), ([1E+15, 0.00, 25000], None, None, None, {0.1: [1E+15, 0.00, 25000], 1.0: [1E+16, 0.00, 25000], 100.0: [1E+18, 0.00, 25000]}, None),)}
param_tup 
description: rate expression for a reaction
type: tuple
format: (highp _params, lowp_ params, troe_params, cheb_dct, plog_dct, collider_dct)
highp_params
description: Arrhenius parameters for the high-pressure limit
type: list
format: [A, n, Ea]
units:
A is on a molar basis
n is relative to a reference temp of 1 Kelvin
Ea is in cal/mol
note: highp_params should only ever contain a single Arrhenius expression (i.e.
lowp_params
description: Arrhenius parameters for the low-pressure limit
only used for Lindemann and Troe expressions
type: list
format: [A, n, Ea]
units:
A is on a molar basis
n is relative to a reference temp of 1 Kelvin
Ea is in cal/mol
note: lowp_params should only ever contain a single Arrhenius expression
troe_params
description: Troe parameters
type: list
format: [alpha, T***, T*, T**]
units: 
alpha is dimensionless
T***, T*, and T** are in Kelvin
T** is optional; it can either be omitted from the array or specified as None
cheb_dct
description: Chebyshev parameters
type: dct
format: {'t_limits': [tmin, tmax], 'p_limits': [pmin, pmax], 'alpha_elm': cheb_coeffs, ‘a_units’: units of the output rate coefficient
units: 
tmin and tmax are in Kelvin
pmin and pmax are in atmospheres
the units of the output rate constant are given by the ‘a_units’ value, which is a str that can be either ‘moles’ or ‘molecules’
cheb_coeffs 
description: Chebyshev polynomial coefficients
type: Numpy array 
shape (N, M), where N is the number of basis functions along the temperature axis and M is the number of basis functions along the pressure axis
Note: N, M is the same order that these parameters are defined in the Chemkin CHEB command
For N=2, M=3, this Numpy array would look like [[a,b,c], [d,e,f]]
format: Numpy_array[[coeff1, coeff2, ...], [...], ...]
units: the units of the output rate constant are given by the ‘a_units’ value, which is a str that can be either ‘moles’ or ‘molecules’ (see the cheb_dct entry)
plog_dct
description: PLOG parameters
type: dct
format: {pressure1: highp _params1, pressure2: ...}
units: pressures are in atmospheres
note: for pressures with more than one Arrhenius expression, duplicates are described by multiple param_tuples (see the rxn_param_dct entry)


Objects describing rate constants

rxn_ktp_dct
description: rate constants for each reaction as a function of T and P
type: dct
format: {rxn1: ktp_dct1, rxn2: ...}
ktp_dct
description: rate constants as a function of T and P
type: dct
format: {pressure1: (temps1, kts1), pressure2: ...}
units: 
pressures are in atmospheres
see the temps and kts entries
temps
description: temperatures
type: Numpy array 
shape (N,), where N is the number of temps
format: Numpy _array[temp1, temp2, ...]
units: Kelvin
kts
description: rate constants as a function of temperature
type: Numpy array 
shape (N,), where N is the number of rate constants
format: Numpy _array[kt1, kt2, ...]
units: mol, cm, s basis, with values determined by molecularity of the reaction


Objects describing thermodynamic polynomial fits

spc_nasa7_dct
type: dct
{spc1: nasa7_params1, spc2: ...}
nasa7_params
type: list
need to work on this...



Objects describing thermodynamic property values

spc_thermo_dct
type: dct
{spc1: thermo_array1, spc2: ...}
thermo_array 
type: lst
list of lists
[temps, h, cp, s, g]
each item is a 1xN numpy array 
is this correct? might need to change; check mechanalyzer/calculator/thermo.py


Objects describing comparisons of mechanisms

aligned_rxn_ktp_dct
type: dct
.. code_block: python
    {rxn1: [ktp_dct1, ktp_dct2, ...], rxn2: ...}
aligned_rxn_ratio_dct
type: dct
.. code_block: python
    {rxn1: [ratio_dct1, ratio_dct2, ...], rxn2: ...}
aligned_spc_thermo_dct
type: dct
.. code_block: python
    {spc1: [thermo_array1, thermo_array2, ...] , spc2: ...}
aligned_spc_diff_dct
type: dct
.. code_block: python
    {spc1: [diff_array1, diff_array2, ...] , spc2: ...}
ratio_dct:
type: dct 
similar structure to a ktp_dct, except give a ratio of k(T,P) values relative to another ktp_dct
{pressure1: (temps1, ratios1), pressure2: ...}



