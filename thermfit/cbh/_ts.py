""" Computes the Heat of Formation at 0 K for a given species
"""

import numpy
import automol.inchi
import automol.graph
import automol.formula
import automol.reac
from automol.par import ReactionClass
from thermfit.cbh import _tsgra as tsutil
from thermfit.cbh import _util as util
from thermfit.cbh._spc import species_basis


# (Reaction Type, IsRadRad)
CBH_TS_CLASSES = [
    (ReactionClass.Typ.HYDROGEN_ABSTRACTION, False),
    (ReactionClass.Typ.HYDROGEN_ABSTRACTION, True),
    (ReactionClass.Typ.ADDITION, False),
    (ReactionClass.Typ.ELIMINATION, False),
    (ReactionClass.Typ.BETA_SCISSION, False)
]


# Main CBH functions to call
def ts_basis(zrxn, scheme):
    """ Get the basis for the appropriate CBH scheme

        :param zrxn: reaction object oriented to Z-Matrix
        :type zrxn: automol.reac.Reaction object
        :param scheme: CBH Scheme used to generate basis
        :type scheme: str
    """
    if scheme == 'basic':
        frag_lst, coeff_lst = basic_ts_basis(zrxn, spc_scheme=scheme)
    else:
        frag_lst, coeff_lst = cbh_basis(zrxn, scheme)

    return frag_lst, coeff_lst


def basic_ts_basis(zrxn, spc_scheme):
    """ Determine a basis for relative enthalpy calculations for a
        transition states by looping over the reactants and building
        a list of a simple species in a inear combination which
        reproduces the number of atoms in the stoichiometry of the species.

        Basis represented by list of InChI strings and array of coefficients.

        :param zrxn: reaction object oriented to Z-Matrix
        :type zrxn: automol.reac.Reaction object
        :param spc_scheme: CBH Scheme for species fragments
        :type spc_scheme: str
        :rtype: (tuple(str), numpy.ndarray)
    """

    # Just use reactants
    rxn_ichs = automol.reac.reaction_inchis(zrxn)
    rct_ichs, _ = rxn_ichs

    basis, coeff_lst = [], []
    for ich in rct_ichs:
        spc_bas_i, coeff_bas_i = species_basis(ich, spc_scheme, balance=True)
        for bas_i, c_bas_i in zip(spc_bas_i, coeff_bas_i):
            if bas_i not in basis:
                # Add basis and coefficients to list
                basis.append(bas_i)
                coeff_lst.append(c_bas_i)
            else:
                # Add coefficient value to existing coefficient value
                for j, bas_j in enumerate(basis):
                    if bas_i == bas_j:
                        coeff_lst[j] += c_bas_i

    return (tuple(basis), numpy.array(coeff_lst))


def cbh_basis(zrxn, scheme):
    """ Determine basis species required for relative enthalpy calculations
        for transition state structures for various CBH-n schemes.

        :param zrxn: reaction object oriented to Z-Matrix
        :type zrxn: automol.reac.Reaction object
        :param scheme: CBH scheme to determine basis for
        :type scheme: str
        :rtype: (tuple(str), list(float))
    """

    zrxn = automol.reac.without_dummy_atoms(zrxn)

    rxnclass = automol.reac.reaction_class(zrxn)
    radrad = automol.reac.is_radical_radical(zrxn)

    frm_bnd_keys = automol.reac.forming_bond_keys(zrxn)
    brk_bnd_keys = automol.reac.breaking_bond_keys(zrxn)
    frm_key1, frm_key2 = tsutil.split_bnd_keys(frm_bnd_keys)
    brk_key1, brk_key2 = tsutil.split_bnd_keys(brk_bnd_keys)
    gra = zrxn.forward_ts_graph
    # Set up graph and reaction site information
    site = None
    site2 = None

    #  Elimination missing the forming double bond
    if rxnclass == ReactionClass.Typ.ELIMINATION:
        # brk_key, brk_key2 = _elimination_find_brk_bnds(gra, frm_key1)
        frm_key2 = tsutil.elimination_second_forming_bond(
            gra, brk_key1, brk_key2)

    #  Addition is missing the 2nd order bond in the graph
    elif rxnclass == ReactionClass.Typ.ADDITION:
        gra, brk_key1 = tsutil.add_appropriate_pi_bonds(gra)
        if not brk_key1:
            gra = tsutil.remove_frm_bnd(gra, brk_key1, frm_key1)
            gra, brk_key1 = tsutil.add_appropriate_pi_bonds(gra)

    # The first set of forming and breaking bonds makes the first reaction site
    if frm_key1 and brk_key1 and rxnclass != ReactionClass.Typ.ELIMINATION:
        site = [
            tsutil.xor(frm_key1, brk_key1),
            tsutil.intersec(frm_key1, brk_key1),
            tsutil.xor(brk_key1, frm_key1)]
    #  eliminations are one large reaction site that we split into
    # site1 and site2 for convieninece
    if rxnclass == ReactionClass.Typ.ELIMINATION:
        try:
            site = [
                tsutil.xor(frm_key1, brk_key1),
                tsutil.intersec(frm_key1, brk_key1),
                tsutil.xor(brk_key1, frm_key1)]
            site2 = [
                tsutil.xor(frm_key2, brk_key2),
                tsutil.intersec(frm_key2, brk_key2),
                tsutil.xor(brk_key2, frm_key2)]
        except:  # noqa: E722
            site = [
                tsutil.xor(frm_key1, brk_key2),
                tsutil.intersec(frm_key1, brk_key2),
                tsutil.xor(brk_key2, frm_key1)]
            site2 = [
                tsutil.xor(frm_key2, brk_key2),
                tsutil.intersec(frm_key2, brk_key2),
                tsutil.xor(brk_key2, frm_key2)]

    elif rxnclass == ReactionClass.Typ.BETA_SCISSION:
        rad_atm = list(automol.graph.sing_res_dom_radical_atom_keys(gra))[0]
        adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
        site = [rad_atm, None, None]
        for atm in brk_key1:
            if rad_atm in adj_atms[atm]:
                site[1] = atm
            else:
                site[2] = atm

    #  radical radical hydrogen abstraction needs a second site
    #  where the pi bond is formed
    elif rxnclass == ReactionClass.Typ.HYDROGEN_ABSTRACTION and radrad:
        rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
        adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)
        atmc, atmd = frm_key1
        if atmc not in rad_atms:
            atmd, atmc = atmc, atmd
        for atma in brk_key1:
            if atma != atmd:
                for atmb in adj_atms[atma]:
                    if atmb in rad_atms:
                        frm_key2 = frozenset({atma, atmb})
                        site2 = [atmb, atma, atmd]

    fclasses = (
        (ReactionClass.Typ.HYDROGEN_ABSTRACTION, False),
        (ReactionClass.Typ.HYDROGEN_MIGRATION, False),
        (ReactionClass.Typ.BETA_SCISSION, False),
        (ReactionClass.Typ.ADDITION, False)
    )
    if rxnclass == ReactionClass.Typ.HYDROGEN_ABSTRACTION and radrad:
        if scheme == 'cbh0':
            frags = cbhzed_radradabs(gra, site, site2)
        elif scheme == 'cbh1':
            frags = cbhone_radradabs(gra, site, site2)
    elif any(cls[0] == rxnclass for cls in fclasses):
        if scheme == 'cbh0':
            frags = cbhzed_habs(gra, site)
        elif scheme == 'cbh1':
            frags = cbhone_habs(gra, site)
    elif rxnclass == ReactionClass.Typ.ELIMINATION:
        if scheme == 'cbh0':
            frags = cbhzed_elim(gra, site, site2)
        elif scheme == 'cbh1':
            frags = cbhone_elim(gra, site, site2)
    else:
        raise NotImplementedError
    # Split the transformed graphs into a list of inchis
    fraglist = []
    clist = []
    for frag_gra in frags.values():
        if 'exp_gra' in frag_gra:
            fraglist.append(automol.graph.inchi(frag_gra['exp_gra']))
            clist.append(frag_gra['coeff'])
        else:
            # if rxnclass == ReactionClass.Typ.HYDROGEN_ABSTRACTION and radrad:
            #     fraglist.append(
            #         tsutil.split_radradabs_gras(frag_gra['ts_gra']))
            # else:
            #     fraglist.append(
            #         tsutil.split_gras(frag_gra['ts_gra']))
            fraglist.append(
                tsutil.split_gras(frag_gra['ts_gra']))
            clist.append(frag_gra['coeff'])
    return fraglist, clist


# Individual CBH-n calculators
def cbhzed_radradabs(gra, site1, site2, bal=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    _, atms, bnd_ords, atm_vals, adj_atms = tsutil.ts_graph(gra, site1, site2)
    # Determine CBHzed fragments
    frags = {}
    for atm in atm_vals:
        grai = (atms.copy(), bnd_ords.copy())
        if (atms[atm][0] != 'H' or atm in site1 or atm in site2):
            if atm in [site1[1], site1[2], site2[0], site2[2]]:
                continue
            coeff = 1.0
            if not bal:
                if atm in site1 + site2:
                    nonhyd_adj_atms1 = tsutil.remove_hyd_from_adj_atms(
                        atms, adj_atms[site1[0]], site2,
                        other_adj=adj_atms[site2[0]])
                    nonhyd_adj_atms2 = tsutil.remove_hyd_from_adj_atms(
                        atms, adj_atms[site2[0]], site1,
                        other_adj=adj_atms[site1[0]])
                    nonhyd_adj_atms1 = tuple(adj for adj in nonhyd_adj_atms1
                                             if adj not in site1)
                    nonhyd_adj_atms1 = tuple(adj for adj in nonhyd_adj_atms2
                                             if adj not in site2)
                    coeff = (
                        util.branch_point(
                            nonhyd_adj_atms1, nonhyd_adj_atms2) *
                        util.terminal_moiety(
                            nonhyd_adj_atms1, nonhyd_adj_atms2)
                    )
                else:
                    nonhyd_adj_atms = tsutil.remove_hyd_from_adj_atms(
                        atms, adj_atms[atm])
                    coeff = (
                        util.branch_point(nonhyd_adj_atms) *
                        util.terminal_moiety(nonhyd_adj_atms)
                    )
            if atm == site1[0]:
                key = 'ts_gra'
                extended_site = [*site1, *site2]
            else:
                key = 'exp_gra'
                extended_site = [atm]
            for site_atm in extended_site:
                for atm_x in adj_atms[site_atm]:
                    if atm_x not in extended_site and atms[atm_x][0] != 'H':
                        grai = cleave_group_and_saturate(
                            grai, site_atm, atm_x)
            grai = automol.graph.explicit(grai)
            frags = _add_frag_to_frags(key, coeff, grai, frags)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags_ts(gra, frags)
    frags = tsutil.simplify_gra_frags(frags)

    return frags


def cbhone_radradabs(gra, site1, site2, bal=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    _, atms, bnd_ords, atm_vals, adj_atms = tsutil.ts_graph(gra, site1, site2)

    # Determine CBHzed fragments
    frags = {}
    for bnd in bnd_ords:
        atma, atmb = bnd
        extended_site = False
        grai = (atms.copy(), bnd_ords.copy())
        if atma not in site1 + site2 or atmb not in site1 + site2:
            coeff = 1.0
            if atmb in site1 + site2:
                atma, atmb = atmb, atma
            if (atma in site1 or atma in site2) and (atms[atmb][0] != 'H'):
                key = 'ts_gra'
                extended_site = [*site1, *site2]
            elif atms[atma][0] != 'H' and atms[atmb][0] != 'H':
                key = 'exp_gra'
                extended_site = [atma, atmb]
            if extended_site:
                for site_atm in extended_site:
                    for atm_x in adj_atms[site_atm]:
                        if atm_x not in extended_site and atms[atm_x][0] != 'H':
                            grai = cleave_group_and_saturate(
                                grai, site_atm, atm_x)
                grai = automol.graph.explicit(grai)
                frags = _add_frag_to_frags(key, coeff, grai, frags)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            if not frags:
                zedfrags = cbhzed_radradabs(gra, site1, site2, bal=True)
            else:
                zedfrags = cbhzed_radradabs(gra, site1, site2, bal=False)
            newfrags = frags.copy()
            for zedfrags_dct in zedfrags.values():
                newname = None
                repeat = False
                if 'exp_gra' in zedfrags_dct:
                    key = 'exp_gra'
                    for onename in newfrags:
                        if 'exp_gra' in newfrags[onename]:
                            if automol.graph.full_isomorphism(
                                   newfrags[onename][key],
                                   zedfrags_dct[key]):
                                newname = onename
                                repeat = True
                else:
                    key = 'ts_gra'
                if not repeat:
                    newname = len(newfrags.keys())
                    newfrags[newname] = {}
                    newfrags[newname][key] = zedfrags_dct[key]
                if not frags:
                    util.add2dic(
                        newfrags[newname], 'coeff', zedfrags_dct['coeff'])
                else:
                    util.add2dic(
                        newfrags[newname], 'coeff', -zedfrags_dct['coeff'])
            frags = newfrags
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags_ts(gra, frags)
    frags = tsutil.simplify_gra_frags(frags)

    return frags


def cbhzed_elim(gra, site1, site2, bal=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    _, atms, bnd_ords, atm_vals, adj_atms = tsutil.ts_graph(gra, site1, site2)

    # Determine CBHzed fragments
    frags = {}
    if not site1[0] == site2[2]:
        site2, site1 = site1, site2
    for atm in atm_vals:
        grai = (atms.copy(), bnd_ords.copy())
        if (atms[atm][0] != 'H' or atm in site1 + site2):
            if atm in [site1[1], site1[2], site2[0], site2[1]]:
                # Dont overcount reactions site
                # Ignore atm if its any of these because we'll count
                # all these through site1[0] and site2[0]
                continue
            coeff = 1.0
            if not bal:
                # The coefficient changes for branch and terminal points if
                # cbhzed is being used to balance cbhx then coeffs are normal
                coeff = _coeff_for_elim_sites(
                    atm, atms, adj_atms, site1, site2)
            if atm == site1[0]:
                key = 'ts_gra'
                extended_site = [*site1, *site2]
            else:
                key = 'exp_gra'
                extended_site = [atm]
            for site_atm in extended_site:
                for atm_x in adj_atms[site_atm]:
                    if atm_x not in extended_site and atms[atm_x][0] != 'H':
                        grai = cleave_group_and_saturate(
                            grai, site_atm, atm_x)

            grai = automol.graph.explicit(grai)
            frags = _add_frag_to_frags(key, coeff, grai, frags)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags_ts(gra, frags)
    frags = tsutil.simplify_gra_frags(frags)

    return frags


def cbhzed_habs(gra, site, bal=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    _, atms, bnd_ords, atm_vals, adj_atms = tsutil.ts_graph(gra, site)
    # Determine CBHzed fragments
    frags = {}
    for atm in atm_vals:
        grai = (atms.copy(), bnd_ords.copy())
        if (atms[atm][0] != 'H' or atm in site):
            if atm in (site[1], site[2]):
                continue
            coeff = 1.0
            if not bal:
                coeff = _coeff_for_ts_sites(
                    atm, atms, adj_atms, site)
            if atm == site[0]:
                key = 'ts_gra'
                extended_site = [*site]
            else:
                key = 'exp_gra'
                extended_site = [atm]
            for site_atm in extended_site:
                for atm_x in adj_atms[site_atm]:
                    if atm_x not in extended_site and atms[atm_x][0] != 'H':
                        grai = cleave_group_and_saturate(
                            grai, site_atm, atm_x)
            grai = automol.graph.explicit(grai)
            frags = _add_frag_to_frags(key, coeff, grai, frags)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags_ts(gra, frags)
    frags = tsutil.simplify_gra_frags(frags)

    return frags


def cbhone_elim(gra, site1, site2, bal=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    _, atms, bnd_ords, atm_vals, adj_atms = tsutil.ts_graph(gra, site1, site2)
    if not site1[0] == site2[2]:
        site2, site1 = site1, site2
    # Determine CBHone fragments
    frags = {}
    for bnd in bnd_ords:
        grai = (atms.copy(), bnd_ords.copy())
        extended_site = None
        atma, atmb = bnd
        if atma not in site1 + site2 or atmb not in site1 + site2:
            coeff = 1.0
            if atmb in site1 + site2:
                atma, atmb = atmb, atma
            if (atma in site1 or atma in site2) and (atms[atmb][0] != 'H'):
                key = 'ts_gra'
                extended_site = [*site1, *site2, atmb]
            elif atms[atma][0] != 'H' and atms[atmb][0] != 'H':
                key = 'exp_gra'
                extended_site = [atma, atmb]

            if extended_site is not None:
                for site_atm in extended_site:
                    for atm_x in adj_atms[site_atm]:
                        if atm_x not in extended_site and atms[atm_x][0] != 'H':
                            grai = cleave_group_and_saturate(
                                grai, site_atm, atm_x)
                grai = automol.graph.explicit(grai)
                frags = _add_frag_to_frags(key, coeff, grai, frags)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            if not frags:
                zedfrags = cbhzed_elim(gra, site1, site2, bal=True)
            else:
                zedfrags = cbhzed_elim(gra, site1, site2, bal=False)
            newfrags = frags.copy()
            for zedfrags_dct in zedfrags.values():
                newname = None
                repeat = False
                if 'exp_gra' in zedfrags_dct:
                    key = 'exp_gra'
                    for onename in newfrags:
                        if 'exp_gra' in newfrags[onename]:
                            if automol.graph.full_isomorphism(
                                    newfrags[onename][key],
                                    zedfrags_dct[key]):
                                newname = onename
                                repeat = True
                else:
                    key = 'ts_gra'
                if not repeat:
                    newname = len(newfrags.keys())
                    newfrags[newname] = {}
                    newfrags[newname][key] = zedfrags_dct[key]
                if not frags:
                    util.add2dic(
                        newfrags[newname], 'coeff',  zedfrags_dct['coeff'])
                else:
                    util.add2dic(
                        newfrags[newname], 'coeff',  -zedfrags_dct['coeff'])
            frags = newfrags
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags_ts(gra, frags)
    frags = tsutil.simplify_gra_frags(frags)
    return frags


def cbhone_habs(gra, site, bal=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    _, atms, bnd_ords, atm_vals, adj_atms = tsutil.ts_graph(gra, site)

    # Determine CBHone fragments
    frags = {}

    for bnd in bnd_ords:
        atma, atmb = bnd
        grai = (atms.copy(), bnd_ords.copy())
        extended_site = None
        if ((atms[atma][0] != 'H' or atma in site) and
           (atms[atmb][0] != 'H' or atmb in site)):
            if atma in site and atmb in site:
                continue
            coeff = 1.0
            if atmb in site:
                atmb, atma = atma, atmb
            if atma in site:
                key = 'ts_gra'
                extended_site = [*site, atmb]
            elif atms[atma][0] != 'H' and atms[atmb][0] != 'H':
                key = 'exp_gra'
                extended_site = [atma, atmb]

            if extended_site is not None:
                for site_atm in extended_site:
                    for atm_x in adj_atms[site_atm]:
                        if atm_x not in extended_site and atms[atm_x][0] != 'H':
                            grai = cleave_group_and_saturate(
                                grai, site_atm, atm_x)

                grai = automol.graph.explicit(grai)
                frags = _add_frag_to_frags(key, coeff, grai, frags)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            if not frags:
                zedfrags = cbhzed_habs(gra, site, bal=True)
            else:
                zedfrags = cbhzed_habs(gra, site, bal=False)
            newfrags = frags.copy()
            for zedfrags_dct in zedfrags.values():
                newname = None
                repeat = False
                if 'exp_gra' in zedfrags_dct:
                    key = 'exp_gra'
                    for onename in newfrags:
                        if 'exp_gra' in newfrags[onename]:
                            if automol.graph.full_isomorphism(
                               newfrags[onename][key],
                               zedfrags_dct[key]):
                                newname = onename
                                repeat = True
                else:
                    key = 'ts_gra'
                if not repeat:
                    newname = len(newfrags.keys())
                    newfrags[newname] = {}
                    newfrags[newname][key] = zedfrags_dct[key]
                if not frags:
                    util.add2dic(
                        newfrags[newname], 'coeff', zedfrags_dct['coeff'])
                else:
                    util.add2dic(
                        newfrags[newname], 'coeff', -zedfrags_dct['coeff'])
            frags = newfrags
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags_ts(gra, frags)
    frags = tsutil.simplify_gra_frags(frags)

    return frags


def _add_frag_to_frags(key, coeff, grai, frags):
    newname = None
    repeat = False
    for name, frags_dct in frags.items():
        if key in frags_dct:
            if key == 'exp_gra':
                if automol.graph.full_isomorphism(
                       frags_dct[key], grai):
                    newname = name
                    repeat = True
            else:
                if frags_dct[key] == grai:
                    newname = name
                    repeat = True
    if not repeat:
        newname = len(frags.keys())
        frags[newname] = {}
        frags[newname][key] = grai
    util.add2dic(frags[newname], 'coeff', coeff)
    return frags


def cleave_group_and_saturate(gra, atmi, atmj):
    """
    Fragments the bnd between atmi-atmj and adds hydrogens
    to atmi for each bond order that was removed. Returns
    the graph that atmi is in
    INPUT:
    gra -- graph
    atmi -- atm we are saturating
    atmj -- atm we are cleaving
    OUTPUT
    gra -- new graph
    """
    # Graphical info about molecule

    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
    order = list(bnd_ords[frozenset({atmi, atmj})])[0]
    for _ in range(order):
        gra = automol.graph.add_bonded_atom(gra, 'H', atmi)
    gra = automol.graph.remove_bonds(gra, (frozenset({atmi, atmj}),))
    gras = automol.graph.connected_components(gra)
    new_gra = None
    for gra_k in gras:
        atmsk = automol.graph.atom_keys(gra_k)
        if atmi in atmsk:
            new_gra = gra_k
    return new_gra


def _coeff_for_elim_sites(atm, atms, adj_atms, site1, site2):
    if atm in site1:
        nonhyd_adj_atms1 = tsutil.remove_hyd_from_adj_atms(
            atms, adj_atms[site1[1]], othersite=site2)
        nonhyd_adj_atms2 = tsutil.remove_hyd_from_adj_atms(
            atms, adj_atms[site2[2]], othersite=site1)
        nonhyd_adj_atms1 = tuple(adj for adj in nonhyd_adj_atms1
                                 if adj not in site1)
        nonhyd_adj_atms2 = tuple(adj for adj in nonhyd_adj_atms2
                                 if adj not in site2)
        coeff = (
            util.branch_point(
                nonhyd_adj_atms1, nonhyd_adj_atms2) *
            util.terminal_moiety(
                nonhyd_adj_atms1, nonhyd_adj_atms2)
        )
    else:
        nonhyd_adj_atms = tsutil.remove_hyd_from_adj_atms(
            atms, adj_atms[atm])
        coeff = (
            util.branch_point(nonhyd_adj_atms) *
            util.terminal_moiety(nonhyd_adj_atms)
        )
    return coeff


def _coeff_for_ts_sites(atm, atms, adj_atms, site):
    if atm in site:
        nonhyd_adj_atms1 = tsutil.remove_hyd_from_adj_atms(
            atms, adj_atms[site[0]], site,
            other_adj=adj_atms[site[2]])
        nonhyd_adj_atms2 = tsutil.remove_hyd_from_adj_atms(
            atms, adj_atms[site[2]],
            site, other_adj=adj_atms[site[0]])
        nonhyd_adj_atms3 = []
        for adj in adj_atms[site[0]]:
            if adj in adj_atms[site[2]]:
                nonhyd_adj_atms3 = tsutil.remove_hyd_from_adj_atms(
                    atms, adj_atms[adj], othersite=site)
        coeff = (
            util.branch_point(
                nonhyd_adj_atms1, nonhyd_adj_atms2,
                nonhyd_adj_atms3) *
            util.terminal_moiety(
                nonhyd_adj_atms1, nonhyd_adj_atms2,
                nonhyd_adj_atms3, endisterm=False)
        )
    else:
        nonhyd_adj_atms = tsutil.remove_hyd_from_adj_atms(
            atms, adj_atms[atm])
        coeff = (
            util.branch_point(nonhyd_adj_atms) *
            util.terminal_moiety(nonhyd_adj_atms)
        )
    return coeff
