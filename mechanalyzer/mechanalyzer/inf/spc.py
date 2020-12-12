"""
  Handles data objects
"""

from mechanalyzer import par


SPC_PROPS = [
    par.SPC.INCHI,
    par.SPC.CHARGE,
    par.SPC.MULT
]


def from_data(inchi, charge, mult):
    """ spc info data structure from necessary data
    """

    # assert to put in a data for a real species
    return (inchi, charge, mult)


def from_dct(dct):
    """ spc info data structure from
    """

    assert set(SPC_PROPS) <= set(dct), (
        'Properties {} not in dict'.format(
            ' '.join(SPC_PROPS))
    )

    inf_obj = tuple()
    for prop in SPC_PROPS:
        inf_obj += (dct[prop],)

    return inf_obj


def value(inf_obj, val):
    """ obtain a value
    """

    assert val in SPC_PROPS, (
        'Desired value {} not in spc info object'.format(
            val)
    )

    return inf_obj[SPC_PROPS.index(val)]


def combine(inf_obj1, inf_obj2, mval='max'):
    """ Created a spc info for two species where charge and mult
        have been combined
    """

    assert mval in ('max', 'min')

    ich = (inf_obj1[0], inf_obj2[0])
    chg = inf_obj1[1] + inf_obj2[1]
    if mval == 'max':
        mult = max([inf_obj1[2], inf_obj2[2]])
    else:
        mult = min([inf_obj1[2], inf_obj2[2]])

    return (ich, chg, mult)
