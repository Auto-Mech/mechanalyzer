""" Format and check things for assessing mechanisms
    and then produce data.
"""

import os

# Format strings


def format_rxn(rxn):
    """ format reaction into a string
    """
    return '{} = {}'.format('+'.join(rxn[0]), '+'.join(rxn[1]))


def formatp(pressure):
    """ format pressure into a string
    """
    if pressure in ('High', 'high'):
        pressure_str = '{0:>10s}'.format(pressure)
    else:
        pressure_str = '{0:>10.3f}'.format(pressure)
    return pressure_str


# Check reaction types
def chk_rxn(rxn, typ, allow_rcts=()):
    """ Check if a reaction meets criteria to
    """
    if typ == 'abstraction':
        nreacs, nprods = 2, 2
    elif typ == 'addition':
        nreacs, nprods = 2, 1
    elif typ == 'isomerization':
        nreacs, nprods = 1, 1
    elif typ == 'decomposition':
        nreacs, nprods = 1, 2
    elif typ == 'ignore':
        chk = True

    if typ != 'ignore':
        [reacs, prods] = rxn
        if len(reacs) == nreacs and len(prods) == nprods:
            if allow_rcts:
                chk = bool(any(rct in reacs for rct in allow_rcts))
            else:
                chk = True
        else:
            chk = False

    return chk


# Deal with files
def read_file(file_name):
    """ Read a file
    """
    print('Read mechanism file:', file_name, '\n')
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


def combine_mech_files(lvl_dirs):
    """ Put mechanism from various stoichiometries and put them in
        a single mechanism file formatted for parsing by analysis
        functions.
    """

    for lvldir in lvl_dirs:

        # Go to path
        os.chdir(lvldir)

        # Build paths
        full_ckin = os.path.join(lvldir, 'full.ckin')
        ckin_files = [os.path.join(lvldir, cdir) for cdir in os.listdir('.')]

        # Write full ckin file using individual files
        with open(full_ckin, 'a') as full_file:
            full_file.write('REACTIONS\n\n')
            for ckin_file in ckin_files:
                with open(ckin_file, 'r') as ckin_file:
                    ckin_str = ckin_file.read()
                full_file.write(ckin_str)
            full_file.write('END')

        # Go back to above dir
        os.chdir('../')
