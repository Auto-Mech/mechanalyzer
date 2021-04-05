"""
  Code to expand the mechanism using stereochemistry
"""


def add_stereo(mech_info, spc_str):
    """ Take a mechanism with reactions and species and
        expand it to include stereoselective reactions.
    """

    # Build an initial dict and expand it to all possible isomerizations
    init_dct = csv_dct(spc_str, values=headers, key='name')

    # 


