""" functions operating on the reactions block string
"""


import itertools
import operator
import numpy as np
import ratefit
from ioformat import phycon
from chemkin_io.parser import reaction as rxn_parser


def mechanism(rxn_block, rxn_units, t_ref, temps, pressures, collider=None,
              ignore_reverse=True, remove_bad_fits=False):
    """ Parses the all the reactions data string in the reaction block
        in a mechanism file for their fitting parameters and
        uses them to calculate rate constants [k(T,P)]s.

        :param rxn_block: string for reaction block from the mechanism input
        :type rxn_block: str
        :param rxn_units: units for parameters specifies
        :type rxn_units: str
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :param collider: bath gas molecule collider
        :type collider: str
        :param ignore_reverse: don't include any reverse reactions
        :type ignore_reverse: bool
        :return: branch_dct: branching fractions for all reactions in mechanism
        :rtype: dict[reaction: branch_ktp_dict]
        :return: total_rate_dct: total k(T,P)s for all reactants in mechanism
        :rtype: dict[reactants: total_ktp_dict]
    """

    # reaction_data_strings = rxn_parser.data_strings(rxn_block)
    reaction_data_dct = rxn_parser.data_dct(
        rxn_block, remove_bad_fits=remove_bad_fits)
    # for rxn in reaction_data_dct:
    #    print('ckin calc rate', rxn)
    mech_dct = {}
    # for rxn, dstrs in reaction_data_dct.items():
    for rxn, dstr in reaction_data_dct.items():
        ktp_dct = reaction(dstr, rxn_units,
                           t_ref, temps, pressures=pressures)
        if rxn not in mech_dct:
            mech_dct[rxn] = ktp_dct
        else:
            mech_dct[rxn] = _add_rates(mech_dct[rxn], ktp_dct)
        # if rxn not in mech_dct and rxn_rev not in mech_dct:
        # rct_names = rxn_parser.reactant_names(dstr)
        # prd_names = rxn_parser.product_names(dstr)
        # rxn = (rct_names, prd_names)
        # rxn_rev = (prd_names, rct_names)
        # elif rxn not in mech_dct and rxn_rev in mech_dct:
        #     if not ignore_reverse:
        #         # have to flip
        #         new_ktp_dct = reaction(dstr, rxn_units,
        #                                t_ref, temps, pressures=pressures)
        #         mech_dct[rxn] = _add_rates(mech_dct[rxn], new_ktp_dct)
        # elif rxn in mech_dct:
        #     new_ktp_dct = reaction(dstr, rxn_units,
        #                            t_ref, temps, pressures=pressures)
        #     mech_dct[rxn] = _add_rates(mech_dct[rxn], new_ktp_dct)
        # if rxn not in mech_dct and rxn_rev not in mech_dct:
        #     mech_dct[rxn] = reaction(dstr, rxn_units,
        #                              t_ref, temps, pressures=pressures)
        # elif rxn in mech_dct and rxn_rev not in mech_dct:
        #     new_ktp_dct = reaction(dstr, rxn_units,
        #                            t_ref, temps, pressures=pressures)
        #     mech_dct[rxn] = _add_rates(mech_dct[rxn], new_ktp_dct)
        # elif rxn not in mech_dct and rxn_rev in mech_dct:
        #     new_ktp_dct = reaction(dstr, rxn_units,
        #                            t_ref, temps, pressures=pressures)
        #     mech_dct[rxn_rev] = _add_rates(mech_dct[rxn_rev], new_ktp_dct)

    # ktp_dct = {}
    # for name, rxn_dstr in rxn_dct.items():
    #     ktp_dct[name] = reaction(
    #         rxn_dstr, rxn_units, t_ref, temps, pressures=pressures)

    return mech_dct


def branching_fractions(mech_dct, pressures):
    """ Parses the all the reactions data string in the reaction block
        in a mechanism file for their fitting parameters and
        uses them to calculate rate constants [k(T,P)]s.
        These rate constants are then used to calculate the
        branching fractions for all the unique reactants in the mechanism.

        :param mech_dct: mechanism dct
        :type mech_dct: dict
        :return: branch_dct: branching fractions for all reactions in mechanism
        :rtype: dict[reaction: branch_ktp_dict]
        :return: total_rate_dct: total k(T,P)s for all reactants in mechanism
        :rtype: dict[reactants: total_ktp_dict]
    """

    # Obtain groups of rxns which share common reactants
    rcts, rct_grps = [], []
    rxns = sorted(mech_dct.keys())
    for key, group in itertools.groupby(rxns, operator.itemgetter(0)):
        rcts.append(key)
        rct_grps.append(list(x for x in group))

    # Build a dct where the rate constants have been combined
    total_rate_dct = {}
    for rct, rct_grp in zip(rcts, rct_grps):
        if len(rct_grp) > 1:
            # Initialize empty dct for all the pressures for below sum to work
            total_rate_dct[rct] = dict(
                zip(pressures, [None for _ in range(len(pressures))]))
            # Sum over all the rates for each reaction, at each pressure
            for pressure in pressures:
                if all(pressure in mech_dct[grp] for grp in rct_grp):
                    total_rate_dct[rct][pressure] = sum(
                        (mech_dct[grp][pressure] for grp in rct_grp))
                else:
                    total_rate_dct[rct][pressure] = None

    # Now get a dct of the branching ration
    branch_dct = {}
    for rxn, rate_dct in mech_dct.items():
        if rxn[0] in total_rate_dct:
            # Initialize empty dct for all the pressures for below sum to work
            branch_dct[rxn] = dict(
                zip(pressures, [None for _ in range(len(pressures))]))
            # Calc ratio: rate / total rate for each reaction, at each pressure
            for pressure in pressures:
                if total_rate_dct[rxn[0]][pressure] is not None:
                    branch_dct[rxn][pressure] = (
                        rate_dct[pressure] / total_rate_dct[rxn[0]][pressure]
                    )
                else:
                    branch_dct[rxn][pressure] = None

    return branch_dct, total_rate_dct


def reaction(rxn_dstr, rxn_units, t_ref, temps, pressures=None, collider=None):
    """ Parses the data string for a reaction for fitting parameters and
        uses those parameters to calculate rate constants at input temps
        and pressures [k(T,P)]s.

        :param rxn_dstr: string for reaction containing fitting parameters
        :type rxn_dstr: str
        :param rxn_units: units for parameters specifies
        :type rxn_units: str
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :param collider: bath gas molecule collider
        :type collider: str
        :return ktp_dct: k(T,P)s at all temps and pressures
        :rtype: dict[pressure: temps]
    """

    # Accepts a params dictionary
    # Read the parameters from the reactions string
    print('dstr\n', rxn_dstr)
    highp_params = rxn_parser.high_p_parameters(rxn_dstr)
    lowp_params = rxn_parser.low_p_parameters(rxn_dstr)
    troe_params = rxn_parser.troe_parameters(rxn_dstr)
    chebyshev_params = rxn_parser.chebyshev_parameters(rxn_dstr)
    plog_params = rxn_parser.plog_parameters(rxn_dstr)
    collid_dct = rxn_parser.collider_enhance_factors(rxn_dstr)

    # Determine if any pdep params at all are found
    any_pdep = any(params is not None
                   for params in (lowp_params, troe_params,
                                  chebyshev_params, plog_params))

    # First check the pressure region that is being specified
    pressure_region = rxn_parser.pressure_region_specification(rxn_dstr)

    # Set the collider efficiency
    collid_factor = collid_dct.get(collider, 1.0)

    print('any_pdep', any_pdep)
    print('p region', pressure_region)

    # Calculate high_pressure rates
    highp_ks = _arrhenius(highp_params, temps, t_ref, rxn_units)
    ktp_dct = {}
    if 'high' in pressures:
        if not any_pdep and pressure_region == 'indep':
            if not rxn_parser.are_highp_fake(highp_params):
                ktp_dct['high'] = highp_ks

    # Get a pdep list of pressures
    pdep_pressures = [pressure for pressure in pressures
                      if pressure != 'high']

    # Calculate pressure-dependent rate constants based on discovered params
    # Either linear Pressure dependence, if specified or using
    # Either (1) Plog, (2) Chebyshev, (3) Lindemann, or (4) Troe
    # Update units if necessary

    pdep_dct = {}
    if pressure_region == 'lowp':
        pdep_dct = ratefit.calc.lowp_limit(
            highp_ks, temps, pdep_pressures, collid_factor=collid_factor)
    else:
        if plog_params is not None:
            pdep_dct = _plog(plog_params, temps, pdep_pressures,
                             t_ref, rxn_units)

        elif chebyshev_params is not None:
            pdep_dct = _chebyshev(chebyshev_params, temps, pdep_pressures)

        elif lowp_params is not None:
            lowp_ks = _arrhenius(lowp_params, temps, t_ref, rxn_units)
            if troe_params is not None:
                pdep_dct = _troe(troe_params, highp_ks, lowp_ks,
                                 temps, pdep_pressures,
                                 collid_factor=collid_factor)
            else:
                pdep_dct = ratefit.calc.lindemann(
                    highp_ks, lowp_ks, temps, pdep_pressures,
                    collid_factor=collid_factor)

    # Build the rate constants dictionary with the pdep dict
    if pdep_dct:
        ktp_dct.update(pdep_dct)

    return ktp_dct


def _add_rates(ktp_dct1, ktp_dct2):
    """ Adds the rates of two dictionaries together.

        :param ktp_dct1: k(T,P)s at all temps and pressures for mechanism 1
        :type ktp_dct1: dict[pressure: temps]
        :param ktp_dct2: k(T,P)s at all temps and pressures for mechanism 2
        :type ktp_dct2: dict[pressure: temps]
        :return ktp_dct1: combined k(T,P)s at all temps and pressures
        :rtype: dict[pressure: temps]
    """

    for pressure in ktp_dct1:
        ktp_dct1[pressure] += ktp_dct2[pressure]

    return ktp_dct1


# Rate calculators
def _arrhenius(arr_params, temps, t_ref, rxn_units):
    """ Calculates rate constants [k(T)]s with the Arrhenius expression
        using the parameters parsed from the reaction string.

        :param arr_params: Arrhenius fitting parameters from string
        :type arr_params: list(float)
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param rxn_units: units for parameters specifies
        :type rxn_units: str
        :return kts: k(T)s at all temps and pressures
        :rtype: numpy.ndarray
    """

    arr_params = _update_params_units(arr_params, rxn_units)
    kts = ratefit.calc.arrhenius(arr_params, t_ref, temps)

    return kts


def _plog(plog_params, temps, pressures, t_ref, rxn_units):
    """ Calculates rate constants [k(T,P)]s with the PLOG expression
        using the parameters parsed from the reaction string.

        :param plog_params: PLOG fitting parameters from string
        :type plog_params: dict[pressure: params]
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param rxn_units: units for parameters specifies
        :type rxn_units: str
        :return ktp_dct: k(T,P)s at all temps and pressures
        :rtype: dict[pressure: temps]
    """

    for pressure, params in plog_params.items():
        plog_params[pressure] = _update_params_units(params, rxn_units)
    ktp_dct = ratefit.calc.plog(plog_params, t_ref, temps, pressures)

    return ktp_dct


def _chebyshev(chebyshev_params, temps, pressures):
    """ Calculates rate constants [k(T,P)]s with the Chebyshev expression
        using the parameters parsed from the reaction string.

        :param chebyshev_params: Chebyshev fitting parameters from string
        :type chebyshev_params: dict[param: val]
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :return ktp_dct: k(T,P)s at all temps and pressures
        :rtype: dict[pressure: temps]
    """

    [tmin, tmax] = chebyshev_params['t_limits']
    [pmin, pmax] = chebyshev_params['p_limits']
    [arows, acols] = chebyshev_params['alpha_dim']
    alpha = np.array(chebyshev_params['alpha_elm'])
    assert alpha.shape == (arows, acols)
    ktp_dct = ratefit.calc.chebyshev(
        alpha, tmin, tmax, pmin, pmax, temps, pressures)

    return ktp_dct


def _troe(troe_params, highp_ks, lowp_ks, temps, pressures, collid_factor=1.0):
    """ Calculates rate constants [k(T,P)]s with the Troe expression
        using the parameters parsed from the reaction string.

        :param troe_params: Troe fitting parameters from string
        :type troe_params: list(float)
        :param highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param lowp_ks: k(T)s determined at low-pressure
        :type lowp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :return ktp_dct: k(T,P)s at all temps and pressures
        :rtype: dict[pressure: temps]
    """

    if len(troe_params) == 3:
        ts2 = None
    elif len(troe_params) == 4:
        ts2 = troe_params[3]
    ktp_dct = ratefit.calc.troe(
        highp_ks, lowp_ks, temps, pressures,
        troe_params[0], troe_params[1], troe_params[2], ts2=ts2,
        collid_factor=collid_factor)

    return ktp_dct


def _update_params_units(params, rxn_units):
    """  Check the units of the Arrhenius fitting parameters
         in the reaction string according to the units given in
         the mechanism file.

        :param params: Arrhenius fitting parameters
        :type params: list(float)
        :param rxn_units: units for parameters specifies
        :type rxn_units: str
        :return params: Arrhenius fitting parameters (A->mol, Ea->kcal.mol-1)
        :rtype: list(float)
    """

    # Determine converstion factor for A parameter units
    if rxn_units[1] == 'molecules':
        a_conv_factor = phycon.NAVO
    else:
        a_conv_factor = 1.0

    # Determine converstion factor for Ea parameter units
    ea_units = rxn_units[0]
    if ea_units == 'cal/mole':
        ea_conv_factor = phycon.CAL2KCAL
    elif ea_units == 'joules/mole':
        ea_conv_factor = phycon.J2KCAL
    elif ea_units == 'kjoules/mole':
        ea_conv_factor = phycon.KJ2KCAL
    elif ea_units == 'kelvin':
        ea_conv_factor = phycon.KEL2KCAL
    else:
        ea_conv_factor = 1.0

    # Set the units of the parameters using the conversion factors
    if params is not None:
        params[0][0] *= a_conv_factor
        params[0][2] *= ea_conv_factor
        if len(params) > 1:
            params[1][0] *= a_conv_factor
            params[1][2] *= ea_conv_factor

    return params
