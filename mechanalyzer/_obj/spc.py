"""
  Handles data objects
"""

from mechanalyzer._obj import par


SPC_PROPS = (
    par.SPC.INCHI,
    par.SPC.CHARGE,
    par.SPC.MULT
)


def from_data(inchi, charge, mult):
    """ spc info data structure from necessary data
    """

    assert isinstance(charge, int), (
        'charge={} is not an integer'.format(charge)
    )
    assert isinstance(mult, int), (
        'multiplicity={} is not an integer'.format(mult)
    )

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
