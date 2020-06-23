"""
Read the mechanism file
"""

import automol
import chemkin_io
from lib.amech_io.parser import ptt


# Build a spc dct containing all info species input file
def build_spc_dct(job_path, spc_type, check_stereo=False):
    """ Get a dictionary of all the input species
        indexed by InChi string
    """
    spc_csv_str = ptt.read_inp_str(job_path, CSV_INP, remove_comments=False)
    if spc_type == 'csv':
        spc_dct = csv_dct(spc_csv_str, check_stereo=check_stereo)
    else:
        raise NotImplementedError

    return mod_spc_dct


# SPC list from reactions
def build_queue(rxn_lst):
    """ Build spc queue from the reaction lst for the drivers
        :return spc_queue: all the species and corresponding models in rxn
        :rtype: list[(species, model),...]
    """

    if 'all' in rxn_lst:
        # First check if rxn_lst is a bunch of species
        spc_queue = rxn_lst['all']['species']
    else:
        # Build the list from expanding the reacs and prods
        spc_queue = []
        for rxn in rxn_lst:
            model = rxn['model']
            spc_queue.extend(((reac, model) for reac in rxn['reacs']))
            spc_queue.extend(((prod, model) for prod in rxn['prods']))

    return spc_queue


# Write new files
def write_stereo_csv(spc_str, path):
    """ read the species file in a .csv format and write a new one
        that has stero information
    """

    # Read in the initial CSV information (deal with mult stereo)
    smi_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'smiles')
    ich_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'inchi')
    mul_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'mult')
    chg_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'charge')
    sens_dct = chemkin_io.parser.mechanism.spc_name_dct(spc_str, 'sens')

    # Build the string for the new file
    spc_str = 'name,SMILES,InChI,IchKey,mult,charge,sens \n'
    for name in ich_dct:
        ich = ich_dct[name]
        smi = smi_dct[name]
        mul = mul_dct[name]
        chg = chg_dct[name]
        sens = sens_dct[name]

        # Generate ichs
        ichs_wstereo = _generate_inchi_stereo(ich)
        for idx, ich_wstereo in enumerate(ich_stereos):
            
            # Augment name if needed
            if idx == 0:
                sname = name
            else:
                sname = name + '-{}'.format(str(idx+2))

            # Generate hash key from InChI
            hashkey = automol.inchi.inchi_key(ich_wstereo)

            # Write to the string 
            spc_str += (
                '{0},\'{1}\',\'{2}\',\'{3}\',{4},{5},{6} \n'.format(
                sname, smi, ich_wstereo, hashkey, mul, chg, sens
            )

    # Write the file
    stereo_file = os.path.join(path, 'species_stereo.csv')
    with open(stereo_path, 'w') as stereo_csv_file:
        stereo_csv_file.write(spc_str)
