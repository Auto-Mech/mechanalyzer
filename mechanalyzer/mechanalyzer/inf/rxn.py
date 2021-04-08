"""
  Reaction info objects
"""

import autofile
from mechanalyzer import par
from mechanalyzer.inf import spc


RXN_PROPS = [
    par.SPC.INCHI,
    par.SPC.CHARGE,
    par.SPC.MULT,
    par.SPC.TSMULT
]


def from_dct(reacs, prods, spc_dct, rxn_mul='low'):
    """ Build a reaction info object using a species dictionary and names
       
        Add the names to the info object?
        Have a way to get a dict object {(rctnames, prdnames) = rxn_info}
    """

    # Build the tuples of the reacs+prods infos
    rxn_ichs, rxn_chgs, rxn_muls = tuple(), tuple(), tuple()
    for side in (reacs, prods):
        ichs, chgs, muls = tuple(), tuple(), tuple()
        for rct in side:
            spc_info = spc.from_dct(spc_dct[rct])
            ichs += (spc.value(spc_info, par.SPC.INCHI),)
            chgs += (spc.value(spc_info, par.SPC.CHARGE),)
            muls += (spc.value(spc_info, par.SPC.MULT),)
        rxn_ichs += (ichs,)
        rxn_chgs += (chgs,)
        rxn_muls += (muls,)

    # Determine the multiplicity of the full reaction
    fake_rxn_info = (rxn_ichs, rxn_chgs, rxn_muls, ())
    ts_mul = ts_mult(fake_rxn_info, rxn_mul=rxn_mul)

    rxn_info = (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)

    return rxn_info


def value(inf_obj, val):
    """ obtain a value
    """

    assert val in RXN_PROPS, (
        'Desired value {} not in rxn info object'.format(val)
    )

    return inf_obj[RXN_PROPS.index(val)]


# write ts mult and charge getter functions
def ts_info(inf_obj):
    """ Build a spc info object for the transisiton state for the reaction
    """

    _chg = ts_chg(inf_obj)
    _mul = value(inf_obj, par.SPC.TSMULT)  # wrong

    return ('', _chg, _mul)


def rgts_info(inf_obj):
    """ obtain all of the info of the rgt info

        get list of spc info objects from a rxn info
    """

    _rgts_info = ()
    for rgt in ('reacs', 'prods'):
        _rgts_info += (rgt_info(inf_obj, rgt),)

    return _rgts_info


def rgt_info(inf_obj, rgt):
    """ obtain all of the info of the rgt info

        get list of spc info objects from a rxn info
    """

    assert rgt in ('reacs', 'prods')

    rxn_ichs, rxn_chgs, rxn_muls, _ = inf_obj
    if rgt == 'reac':
        rgt_ichs, rgt_chgs, rgt_muls = rxn_ichs[0], rxn_chgs[0], rxn_muls[0]
    else:
        rgt_ichs, rgt_chgs, rgt_muls = rxn_ichs[1], rxn_chgs[1], rxn_muls[1]

    _rgt_info = ()
    for rgt_ich, rgt_chg, rgt_mul in zip(rgt_ichs, rgt_chgs, rgt_muls):
        _rgt_info += (spc.from_data(rgt_ich, rgt_chg, rgt_mul),)

    return _rgt_info


def replace(inf_obj, val):
    """ Replace some value of the info object
    """
    raise NotImplementedError


def sort(inf_obj, scheme='autofile'):
    """ Resort the reacs and prods based on some scheme,
        currently just autofile
    """

    rxn_ichs, rxn_chgs, rxn_muls, ts_mul = inf_obj

    if scheme == 'autofile':
        rxn_ichs, rxn_chgs, rxn_muls = autofile.schema.sort_together(
            rxn_ichs, rxn_chgs, rxn_muls)

    sort_inf_obj = (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)

    return sort_inf_obj


def reverse(inf_obj):
    """ Reverse the reactants and products in the info object
    """

    rxn_ichs, rxn_chgs, rxn_muls, ts_mul = inf_obj

    rxn_ichs = (rxn_ichs[1], rxn_ichs[0])
    rxn_chgs = (rxn_chgs[1], rxn_chgs[0])
    rxn_muls = (rxn_muls[1], rxn_muls[0])

    rev_inf_obj = (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)

    return rev_inf_obj


def ts_chg(inf_obj):
    """ Build the transition charge
    """

    rxn_chgs = value(inf_obj, par.SPC.CHARGE)

    _chg = 0
    for rct_chg in rxn_chgs[0]:
        _chg += rct_chg

    return _chg


def ts_mult(inf_obj, rxn_mul='low'):
    """ evaluate the ts multiplicity from the multiplicities
        of the reactants and products
    """

    rxn_muls = value(inf_obj, par.SPC.MULT)

    nrcts, nprds = len(rxn_muls[0]), len(rxn_muls[1])

    # Set the multiplicities
    rct_spin_sum, prd_spin_sum = 0, 0
    rct_muls, prd_muls = [], []
    if nrcts == 1 and nprds == 1:
        _mul = max(rxn_muls[0][0], rxn_muls[1][0])
    else:
        for rct_mul in rxn_muls[0]:
            rct_spin_sum += (rct_mul - 1.)/2.
            rct_muls.append(rct_mul)
        for prd_mul in rxn_muls[1]:
            prd_spin_sum += (prd_mul - 1.)/2.
            prd_muls.append(prd_mul)

        if rxn_mul == 'low':
            _mul = min(rct_spin_sum, prd_spin_sum)
        else:
            _mul = max(rct_spin_sum, prd_spin_sum)
        _mul = int(round(2*_mul + 1))

    return _mul
