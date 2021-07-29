""" Script to generate files with reactions
"""

import automol
import chemkin_io
import mechanalyzer
import thermfit


# Get a list of InChI strings for the reactants; make spc_dct for code
RCT_SMIS = [
    ('C4H7O2(1)', 'O=CCCC[O]'),
    ('C4H7O2(2)', 'CC(=O)CC[O]'),
    ('C4H7O2(3)', 'CCC(=O)C[O]'),
    ('C4H7O2(4)', 'O=CCC([O])C'),
    ('C4H7O2(5)', 'CC(=O)C([O])C')
]

MECH_SPC_DCT = {
    smi[0]: thermfit.create_spec(automol.smiles.inchi(smi[1]))
    for smi in RCT_SMIS}
MECH_RXN_DCT = {}

# Set reaction type
RSERIES = (
    ('all', (automol.par.ReactionClass.Typ.HYDROGEN_MIGRATION,)),
)

# Set headers for csv string
HEADERS = ('smiles', 'inchi', 'mult', 'charge')

# Generate the reactions that are specified
spc_dct, rxn_dct = mechanalyzer.builder.rxn.build_mechanism(
    MECH_SPC_DCT, MECH_RXN_DCT, rxn_series=RSERIES)

# Write the dictionaries to strings
csv_str = mechanalyzer.parser.spc.csv_str(spc_dct, HEADERS)
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=None,
    spc_dct=spc_dct,
    spc_nasa7_dct=None,
    rxn_param_dct=rxn_dct,
    comments=None)

with open('newspc.csv', 'w') as fobj:
    fobj.write(csv_str)
with open('newrxn.dat', 'w') as fobj:
    fobj.write(mech_str)
