"""
 Util stuff, probably move
"""

def write_file(path, outname):
    """ Write a file to path
    """

    basis_file = os.path.join(path, outname)
    with open(basis_file, 'w') as file_obj:
        file_obj.write(spc_str)
