import numpy as np
TEMPS = np.array([500, 1000, 1500])
KTS = np.array([1e10, 1e11, 1e12])
h = np.array([-5000, -4000, -3000])
cp = np.array([-5000, -4000, -3000])
s = np.array([-5000, -4000, -3000])

SPC_IDENT_DCT1 = {
    'H': {'smiles': '', 'inchi': 'InChI=1S/H', 'inchikey': '', 'mult': 2, 'charge': 0, 'sens': 0, 'fml': {'H': 1}},
    'OH': {'smiles': '', 'inchi': 'InChI=1S/HO/h1H', 'inchikey': '', 'mult': 2, 'charge': 0, 'sens': 0, 'fml': {'H': 1, 'O': 1}},
    'O': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 3, 'charge': 0, 'sens': 0, 'fml': {'O': 1}},
    'H2': {'smiles': '', 'inchi': 'InChI=1S/H2/h1H', 'inchikey': '', 'mult': 1, 'charge': 0, 'sens': 0, 'fml': {'H': 2}},
    'O2': {'smiles': '', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': '', 'mult': 1, 'charge': 0, 'sens': 0, 'fml': {'O': 2}},
    'O(S)': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 1, 'charge': 0, 'sens': 0, 'fml': {'O': 1}},
}

SPC_IDENT_DCT2 = {
    'HV': {'smiles': '', 'inchi': 'InChI=1S/H', 'inchikey': '', 'mult': 2, 'charge': 0, 'sens': 0, 'fml': {'H': 1}},
    'OHV': {'smiles': '', 'inchi': 'InChI=1S/HO/h1H', 'inchikey': '', 'mult': 2, 'charge': 0, 'sens': 0, 'fml': {'H': 1, 'O': 1}},
    'OV': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 3, 'charge': 0, 'sens': 0, 'fml': {'O': 1}},
    'O2V': {'smiles': '', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': '', 'mult': 1, 'charge': 0, 'sens': 0, 'fml': {'O': 2}},
    'H2V': {'smiles': '', 'inchi': 'InChI=1S/H2/h1H', 'inchikey': '', 'mult': 1, 'charge': 0, 'sens': 0, 'fml': {'H': 2}},
    'HO2V': {'smiles': '', 'inchi': 'InChI=1S/HO2/c1-2/h1H', 'inchikey': '', 'mult': 2, 'charge': 0, 'sens': 0, 'fml': {'H': 1, 'O': 2}},
    'O(S)V': {'smiles': '', 'inchi': 'InChI=1S/O', 'inchikey': '', 'mult': 1, 'charge': 0, 'sens': 0, 'fml': {'O': 1}},
}

SPC_THERMO_DCT1 = {
    'H': [TEMPS, h, cp, s, np.array([38112.076194796915, 22159.224380360094, 4906.62142234519])],
    'OH': [TEMPS, h, cp, s, np.array([-13437.71616429529, -38594.34206874827, -65664.834890985])],
    'O': [TEMPS, h, cp, s, np.array([40013.92863513222, 18462.12185893287, -4399.717127561729])],
    'H2': [TEMPS, h, cp, s, np.array([-16011.013763631481, -34787.09040495749, -55449.572706248])],
    'O2': [TEMPS, h, cp, s, np.array([-24919.542351752494, -52791.48721222688, -82818.01126124])],
    'O(S)': [TEMPS, h, cp, s, np.array([41013.92863513222, 19462.12185893287, -4899.7171261729])],
}

SPC_THERMO_DCT2 = {
    'HV': [TEMPS, h, cp, s, np.array([38112.076194796915, 22159.224380360094, 4906.621422734519])],
    'OHV': [TEMPS, h, cp, s, np.array([-13437.71616429529, -38594.34206874827, -65664.80890985])],
    'OV': [TEMPS, h, cp, s, np.array([40013.92863513222, 18462.12185893287, -4399.717127561729])],
    'O2V': [TEMPS, h, cp, s, np.array([-24919.542351752494, -52791.48721222688, -82818.01121114])],
    'H2V': [TEMPS, h, cp, s, np.array([-16011.013763631481, -34787.09040495749, -55449.57270738])],
    'HO2V': [TEMPS, h, cp, s, np.array([-24930.759856779678, -56555.863471982615, -91105.70559])],
    'O(S)V': [TEMPS, h, cp, s, np.array([41013.92863513222, 19462.12185893287, -4899.712756729])],
}

RXN_KTP_DCT1 = {
    (('H2', 'O'), ('OH', 'H'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H', 'O2'), ('OH', 'O'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H2', 'O'), ('OH', 'OH'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H', 'O'), ('OH',), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H', 'O'), ('OH',), ('(+M)',)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H', 'O'), ('OH',), ('+O(S)',)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H2', 'O(S)'), ('OH', 'O'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
}

RXN_KTP_DCT2 = {
    (('H2V', 'OV'), ('OHV', 'HV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('OV', 'OHV'), ('O2V', 'HV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H2V', 'O2V'), ('HO2V', 'HV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('OHV',), ('HV', 'OV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('OHV',), ('HV', 'OV'), ('(+M)',)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('OHV',), ('HV', 'OV'), ('+O(S)V',)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
    (('H2V', 'O(S)V'), ('OV', 'OHV'), (None,)): {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
}

