""" Script to generate files with reactions
"""

import automol
from automol.par import ReactionClass
import chemkin_io
import mechanalyzer
import thermfit


# Get a list of InChI strings for the reactants; make spc_dct for code
MECH_SPC_DCT = {
    'C4H10': thermfit.create_spec(automol.smiles.inchi('CCCC')),
    'O2': thermfit.create_spec(automol.smiles.inchi('O=O')),
    'OH': thermfit.create_spec(automol.smiles.inchi('[OH]')),
    'H': thermfit.create_spec(automol.smiles.inchi('[H]'))
}
MECH_RXN_DCT = {}

# Set reaction type
RSERIES = (
    # Series 1
    ('all',
     ('H', 'OH', 'O2'),
     (ReactionClass.Typ.HYDROGEN_ABSTRACTION,)),
    # Series 2
    ('radicals',
     ('O2',),
     (ReactionClass.Typ.ADDITION,
      ReactionClass.Typ.HYDROGEN_MIGRATION,
      ReactionClass.Typ.BETA_SCISSION)),
    # Series 3
    ('all',
     (),
     (ReactionClass.Typ.ELIMINATION,
      ReactionClass.Typ.HYDROGEN_MIGRATION)),
    # Series 4
    ('radicals',
     ('O2',),
     (ReactionClass.Typ.ADDITION,
      ReactionClass.Typ.BETA_SCISSION,
      ReactionClass.Typ.RING_FORM_SCISSION)),
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
