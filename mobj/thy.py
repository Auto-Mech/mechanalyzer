"""
  Handles data objects
"""

THY_PROPS = [ 
    mechanalyzer.THY.PROGRAM,
    mechanalyzer.THY.METHOD
    mechanalyzer.THY.BASIS,
    mechanalyzer.THY.ORB_RESTRICT
]


# Constructors
def from_data(program, method, basis, orb_restrict):
    """ thy info data structure from necessary data 
    """
   
    # assert to put in a data for a real species
    return (program, method, basis, orb_restrict)


def from_dct(dct):
    """ thy info data structure from 
    """

    assert set(THY_PROPS) <= set(dct), (
        'Properties {} not in dict'.format(
            ' '.join(THY_PROPS))
    )

    inf_obj = tuple()
    for prop in THY_PROPS:
        inf_obj += (prop,)

    return inf_dct


def modify_orb_label(thy_info, mult):
    """ Convert the orbital restriction label to one for a species
        using the multiplicity.
    """

    mod_thy_info = thy_info[:3]
    orb_label = value(thy_info, mechanalyzer.THY.ORB_RESTRICT)
    mod_thy_info += (_mod_orbital_label(orb_label, mult),)

    return mod_thy_info


def _mod_orbital_label(orb_label, mult):
    """ orbital restriction logical
    """

    if orb_label == 'RR':
        mod_orb_label = 'R'
    elif orb_label == 'UU':
        mod_orb_label = 'U'
    elif orb_label == 'RU':
        if mult == 1:
            mod_orb_label = 'R'
        else:
            mod_orb_label = 'U'

    return mod_orb_label


# Getters
def value(inf_obj, val):
    """ obtain a value
    """

    assert val in THY_PROPS, (
        'Desired value {} not in thy info object'.format(
            val)
    )

    return inf_obj[THY_PROPS.index(val)]
