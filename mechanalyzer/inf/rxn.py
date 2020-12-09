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
    'rxnmul'
]


def from_dct(reacs, prods, spc_dct, rxn_mul='low', sort=False):
    """ prepare rxn info and reverse the reactants and products
        if reaction is endothermic
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
    _, ts_mul_low, ts_mul_high = rxn_chg_mult(rxn_muls, rxn_chgs)
    ts_mul = ts_mul_low if rxn_mul == 'low' else ts_mul_high

    rxn_info = (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)
    if sort:
        rxn_info = sort(rxn_info)

    return rxn_info


def value(inf_obj, val):
    """ obtain a value
    """

    assert val in RXN_PROPS, (
        'Desired value {} not in rxn info object'.format(val)
    )

    return inf_obj[RXN_PROPS.index(val)]


def sort(inf_obj):
    """ Resort the reacs and prods based on some scheme,
        currently just autofile
    """

    rxn_ichs, rxn_chgs, rxn_muls, ts_mul = inf_obj

    rxn_ichs, rxn_chgs, rxn_muls = autofile.schema.sort_together(
        rxn_ichs, rxn_chgs, rxn_muls)

    sort_inf_obj = (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)

    return sort_inf_obj


def rxn_chg_mult(rxn_muls, rxn_chgs):
    """ evaluate the ts multiplicity from the multiplicities
        of the reactants and products
    """
    nrcts, nprds = len(rxn_muls[0]), len(rxn_muls[1])

    # Set the multiplicities
    rct_spin_sum, prd_spin_sum = 0, 0
    rct_muls, prd_muls = [], []
    if nrcts == 1 and nprds == 1:
        ts_mul_low = max(rxn_muls[0][0], rxn_muls[1][0])
        ts_mul_high = ts_mul_low
    else:
        for rct_mul in rxn_muls[0]:
            rct_spin_sum += (rct_mul - 1.)/2.
            rct_muls.append(rct_mul)
        for prd_mul in rxn_muls[1]:
            prd_spin_sum += (prd_mul - 1.)/2.
            prd_muls.append(prd_mul)
        ts_mul_low = min(rct_spin_sum, prd_spin_sum)
        ts_mul_low = int(round(2*ts_mul_low + 1))
        ts_mul_high = max(rct_spin_sum, prd_spin_sum)
        ts_mul_high = int(round(2*ts_mul_high + 1))

    # Set the charges
    ts_chg = 0
    for rct_chg in rxn_chgs[0]:
        ts_chg += rct_chg

    return ts_chg, ts_mul_low, ts_mul_high
