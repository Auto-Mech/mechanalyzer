"""
  Handles data objects
"""

SPC_PROPS = [
    mechanalyzer.SPC.INCHI,
    mechanalyzer.SPC.CHARGE,
    mechanalyzer.SPC.MULT
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
        inf_obj += (prop,)

    return inf_obj


def value(inf_obj, val):
    """ obtain a value
    """

    assert val in SPC_PROPS, (
        'Desired value {} not in spc info object'.format(
            val)
    )

    return inf_obj[SPC_PROPS.index(val)]
