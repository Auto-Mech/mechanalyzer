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


CBH_TS_CLASSES = [
    'hydrogen abstraction',
    'beta scission',
    'elimination high',
    'radical radical hydrogen abstraction high',
    'addition high'
]


# Main CBH functions to call
def ts_basis(zrxn, scheme, spc_scheme=None):
    """ Get the basis for the appropriate CBH scheme

        :param zrxn: reaction object oriented to Z-Matrix
        :type zrxn: automol.reac.Reaction object
        :param scheme: CBH Scheme used to generate basis
        :type scheme: str
        :param spc_scheme: CBH Scheme for species fragments used in `basic` TS scheme
        :type spc_scheme: str

        spc_scheme either basic, cbh0, cbh1
    """

    if scheme == 'basic':
        # For a cbh_m scheme, spc = cbh_m
        # For a cbh_m_n scheme, spc = cbh_n
        if spc_scheme is None:
            spc_scheme = scheme
            if '_' in spc_scheme:
                spc_scheme = 'cbh' + spc_scheme.split('_')[1]
        frag_lst, coeff_lst = basic_ts_basis(zrxn, spc_scheme=spc_scheme)
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

    gra = tsutil.remove_frm_bnd(gra, brk_key1, frm_key1)
    gra = tsutil.remove_frm_bnd(gra, brk_key2, frm_key2)

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
        except:
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
    elif 'radical radical hyd' in rxnclass:
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

    # Run the appropriate function to transform the graph so that
    # breaking bonds are order N.6 and forming are N.4 and the valences
    # of the atoms involved are appropriately altered
    fclasses = ('hydrogen abstraction', 'beta scission',
                'hydrogen migration', 'addition')
    if 'radical radical hyd' in rxnclass:
        if scheme == 'cbh0':
            frags = cbhzed_radradabs(gra, site, site2)
        elif scheme == 'cbh1':
            frags = cbhone_radradabs(gra, site, site2)
    elif any(cls in rxnclass for cls in fclasses):
        if scheme == 'cbh0':
            frags = cbhzed_habs(gra, site)
        elif scheme == 'cbh1':
            frags = cbhone_habs(gra, site)
    elif 'elimination' in rxnclass:
        if scheme == 'cbh0':
            frags = cbhzed_elim(gra, site, site2)
        elif scheme == 'cbh1':
            frags = cbhone_elim(gra, site, site2)
    else:
        raise NotImplementedError

    # Split the transformed graphs into a list of inchis
    fraglist = []
    clist = []
    for frag in frags:
        if 'exp_gra' in frags[frag]:
            # Remove dummy atoms from graph, broke for H-ABS (KBM)
            # egra = automol.graph.without_dummy_atoms(frags[frag]['exp_gra'])
            # if egra != ({}, {}):  # Check if graph is empty
            #     fraglist.append(automol.graph.inchi(frags[frag]['exp_gra']))
            #     clist.append(frags[frag]['coeff'])
            fraglist.append(automol.graph.inchi(frags[frag]['exp_gra']))
            clist.append(frags[frag]['coeff'])
        else:
            if 'beta' in rxnclass:
                fraglist.append(
                    tsutil.split_beta_gras(frags[frag]['ts_gra']))
            elif 'elim' in rxnclass:
                fraglist.append(
                    tsutil.split_elim_gras(frags[frag]['ts_gra']))
            elif 'radical radical hyd' in rxnclass:
                fraglist.append(
                    tsutil.split_radradabs_gras(frags[frag]['ts_gra']))
            else:
                fraglist.append(
                    tsutil.split_gras(frags[frag]['ts_gra']))
            clist.append(frags[frag]['coeff'])

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
                key1 = [site1[0], site1[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic = {0: (atms[atm][0], atm_vals[atm]-bnd_ord1, None)}
                bnd_dct = {frozenset({0, 1}): (bnd_ord1, None)}
                key1 = [site1[2], site1[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0] + 0.6
                atm_dic[1] = (
                    atms[site1[1]][0],
                    atm_vals[site1[1]]-bnd_ord1-bnd_ord2, None)
                atm_dic[2] = (
                    atms[site1[2]][0],
                    atm_vals[site1[2]]-bnd_ord2, None)
                bnd_dct[frozenset({1, 2})] = (bnd_ord2, None)
                key1 = [site2[0], site2[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic[3] = (
                    atms[site2[0]][0], atm_vals[site2[0]]-bnd_ord1, None)
                atm_dic[2] = (atms[site1[2]][0], atm_dic[2][1]-bnd_ord1, None)
                bnd_dct[frozenset({3, 2})] = (bnd_ord1, None)
            else:
                bnd_dct = {}
                atm_dic = {0: (atms[atm][0], int(atm_vals[atm]), None)}
            grai = (atm_dic, bnd_dct)
            try:
                grai = automol.graph.explicit(grai)
                key = 'exp_gra'
            except:
                key = 'ts_gra'
            newname = None
            repeat = False
            for name in frags:
                if key in frags[name]:
                    if automol.graph.full_isomorphism(frags[name][key], grai):
                        newname = name
                        repeat = True
            if not repeat:
                newname = len(frags.keys())
                frags[newname] = {}
                frags[newname][key] = grai
            util.add2dic(frags[newname], 'coeff', coeff)
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
    _, atms, bnd_ords, atm_vals, _ = tsutil.ts_graph(gra, site1, site2)

    # Determine CBHzed fragments
    frags = {}
    for bnd in bnd_ords:
        atma, atmb = bnd
        unique_bond = False
        if atma not in site1 + site2 or atmb not in site1 + site2:
            coeff = 1.0
            if atmb in site1 + site2:
                atma, atmb = atmb, atma
            if (atma in site1 or atma in site2) and (atms[atmb][0] != 'H'):
                key1 = [site1[0], site1[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic = {
                    0: (atms[site1[0]][0], atm_vals[site1[0]]-bnd_ord1, None)}
                bnd_dct = {frozenset({0, 1}): (bnd_ord1, None)}
                key1 = [site1[2], site1[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0] + 0.6
                atm_dic[1] = (
                    atms[site1[1]][0],
                    atm_vals[site1[1]]-bnd_ord1-bnd_ord2, None)
                atm_dic[2] = (
                    atms[site1[2]][0], atm_vals[site1[2]]-bnd_ord2, None)
                bnd_dct[frozenset({1, 2})] = (bnd_ord2, None)
                key1 = [site2[0], site2[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic[3] = (
                    atms[site2[0]][0], atm_vals[site2[0]]-bnd_ord1, None)
                atm_dic[2] = (atms[site1[2]][0], atm_dic[2][1]-bnd_ord1, None)
                bnd_dct[frozenset({3, 2})] = (bnd_ord1, None)
                key1 = [atma, atmb]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic[4] = (atms[atmb][0], atm_vals[atmb]-bnd_ord1, None)
                if atma == site1[0]:
                    atm_dic[0] = (
                        atm_dic[0][0], atm_dic[0][1] - bnd_ord1, None)
                    bnd_dct[frozenset({0, 4})] = (bnd_ord1, None)
                elif atma == site2[0]:
                    atm_dic[3] = (
                        atm_dic[3][0], atm_dic[3][1] - bnd_ord1, None)
                    bnd_dct[frozenset({3, 4})] = (bnd_ord1, None)
                elif atma == site1[1]:
                    atm_dic[1] = (
                        atm_dic[1][0], atm_dic[1][1] - bnd_ord1, None)
                    bnd_dct[frozenset({1, 4})] = (bnd_ord1, None)
                elif atma == site1[2]:
                    atm_dic[2] = (
                        atm_dic[2][0], atm_dic[2][1] - bnd_ord1, None)
                    bnd_dct[frozenset({2, 4})] = (bnd_ord1, None)
                unique_bond = True
            elif atms[atma][0] != 'H' and atms[atmb][0] != 'H':
                bnd_ord = list(bnd_ords[bnd])[0]
                atm_dic = {
                    0: (atms[atma][0], int(atm_vals[atma])-bnd_ord, None)}
                atm_dic[1] = (atms[atmb][0], int(atm_vals[atmb])-bnd_ord, None)
                bnd_dct = {frozenset({0, 1}): (bnd_ord, None)}
                unique_bond = True
            if unique_bond:
                grai = (atm_dic, bnd_dct)
                try:
                    grai = automol.graph.explicit(grai)
                    key = 'exp_gra'
                except:
                    key = 'ts_gra'
                newname = None
                repeat = False
                for name in frags:
                    if key in frags[name]:
                        if key == 'exp_gra':
                            if automol.graph.full_isomorphism(
                                    frags[name][key], grai):
                                newname = name
                                repeat = True
                        else:
                            if frags[name][key] == grai:
                                newname = name
                                repeat = True
                if not repeat:
                    newname = len(frags.keys())
                    frags[newname] = {}
                    frags[newname][key] = grai
                util.add2dic(frags[newname], 'coeff', coeff)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            zedfrags = cbhzed_radradabs(gra, site1, site2, bal=False)
            newfrags = frags.copy()
            for zedname in zedfrags:
                newname = None
                repeat = False
                if 'exp_gra' in zedfrags[zedname]:
                    key = 'exp_gra'
                    for onename in newfrags:
                        if 'exp_gra' in newfrags[onename]:
                            if automol.graph.full_isomorphism(
                                   newfrags[onename][key],
                                   zedfrags[zedname][key]):
                                newname = onename
                                repeat = True
                else:
                    key = 'ts_gra'
                if not repeat:
                    newname = len(newfrags.keys())
                    newfrags[newname] = {}
                    newfrags[newname][key] = zedfrags[zedname][key]
                util.add2dic(
                    newfrags[newname], 'coeff', -zedfrags[zedname]['coeff'])
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
        if (atms[atm][0] != 'H' or atm in site1 or atm in site2):
            if atm not in [site1[1], site1[2], site2[0], site2[1]]:
                continue
            coeff = 1.0
            if not bal:
                if atm in site1:
                    nonhyd_adj_atms1 = tsutil.remove_hyd_from_adj_atms(
                        atms, adj_atms[site1[1]], othersite=site2)
                    nonhyd_adj_atms2 = tsutil.remove_hyd_from_adj_atms(
                        atms, adj_atms[site2[2]], othersite=site1)
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
                key1 = [site1[0], site1[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic = {0: (atms[atm][0], atm_vals[atm]-bnd_ord1, None)}
                bnd_dct = {frozenset({0, 1}): (bnd_ord1, None)}
                key1 = [site1[2], site1[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0]
                atm_dic[2] = (
                    atms[site1[2]][0], atm_vals[site1[2]]-bnd_ord2-1, None)
                atm_dic[1] = (
                    atms[site1[1]][0],
                    atm_vals[site1[1]]-bnd_ord1-bnd_ord2, None)
                bnd_dct[frozenset({1, 2})] = (bnd_ord2, None)
                key1 = [site2[0], site2[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic[3] = (
                    atms[site2[0]][0], atm_vals[site2[0]]-bnd_ord1-1, None)
                bnd_dct[frozenset({3, 4})] = (bnd_ord1, None)
                key1 = [site2[2], site2[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0]
                atm_dic[4] = (
                    atms[site2[1]][0],
                    atm_vals[site2[1]]-bnd_ord1-bnd_ord2, None)
                atm_dic[0] = (
                    atms[site2[2]][0], list(atm_dic[0])[1]-bnd_ord2, None)
                bnd_dct[frozenset({4, 0})] = (bnd_ord2, None)
                bnd_dct[frozenset({2, 3})] = (1, None)
            else:
                bnd_dct = {}
                atm_dic = {0: (atms[atm][0], int(atm_vals[atm]), None)}
            grai = (atm_dic, bnd_dct)
            try:
                grai = automol.graph.explicit(grai)
                key = 'exp_gra'
            except:
                key = 'ts_gra'
            newname = None
            repeat = False
            for name in frags:
                if key in frags[name]:
                    if automol.graph.full_isomorphism(frags[name][key], grai):
                        newname = name
                        repeat = True
            if not repeat:
                newname = len(frags.keys())
                frags[newname] = {}
                frags[newname][key] = grai
            util.add2dic(frags[newname], 'coeff', coeff)
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
        if (atms[atm][0] != 'H' or atm in site):
            if atm in (site[1], site[2]):
                continue
            coeff = 1.0
            if not bal:
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
            if atm == site[0]:
                key1 = [site[0], site[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic = {0: (atms[atm][0], atm_vals[atm]-bnd_ord1, None)}
                bnd_dct = {frozenset({0, 1}): (bnd_ord1, None)}
                key1 = [site[2], site[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0]
                atm_dic[2] = (
                    atms[site[2]][0], atm_vals[site[2]]-bnd_ord2, None)
                atm_dic[1] = (
                    atms[site[1]][0],
                    atm_vals[site[1]]-bnd_ord1-bnd_ord2, None)
                bnd_dct[frozenset({1, 2})] = (bnd_ord2, None)
            else:
                bnd_dct = {}
                atm_dic = {0: (atms[atm][0], int(atm_vals[atm]), None)}
            grai = (atm_dic, bnd_dct)
            try:
                grai = automol.graph.explicit(grai)
                key = 'exp_gra'
            except:
                key = 'ts_gra'
            newname = None
            repeat = False
            for name in frags:
                if key in frags[name]:
                    if automol.graph.full_isomorphism(frags[name][key], grai):
                        newname = name
                        repeat = True
            if not repeat:
                newname = len(frags.keys())
                frags[newname] = {}
                frags[newname][key] = grai
            util.add2dic(frags[newname], 'coeff', coeff)
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
    _, atms, bnd_ords, atm_vals, _ = tsutil.ts_graph(gra, site1, site2)
    # Determine CBHone fragments
    frags = {}
    if not site1[0] == site2[2]:
        site2, site1 = site1, site2
    for bnd in bnd_ords:
        atma, atmb = bnd
        unique_bond = False
        if atma not in site1 + site2 or atmb not in site1 + site2:
            coeff = 1.0
            if atmb in site1 + site2:
                atma, atmb = atmb, atma
            if (atma in site1 or atma in site2) and (atms[atmb][0] != 'H'):
                key1 = [site1[0], site1[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic = {
                    0: (atms[site1[0]][0], atm_vals[site1[0]]-bnd_ord1, None)}
                bnd_dct = {frozenset({0, 1}): (bnd_ord1, None)}
                key1 = [site1[2], site1[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0]
                atm_dic[2] = (
                    atms[site1[2]][0], atm_vals[site1[2]]-bnd_ord2-1, None)
                atm_dic[1] = (
                    atms[site1[1]][0],
                    atm_vals[site1[1]]-bnd_ord1-bnd_ord2, None)
                bnd_dct[frozenset({1, 2})] = (bnd_ord2, None)
                key1 = [site2[0], site2[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic[3] = (
                    atms[site2[0]][0], atm_vals[site2[0]]-bnd_ord1-1, None)
                bnd_dct[frozenset({3, 4})] = (bnd_ord1, None)
                key1 = [site2[2], site2[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0]
                atm_dic[4] = (
                    atms[site2[1]][0],
                    atm_vals[site2[1]]-bnd_ord1-bnd_ord2, None)
                atm_dic[0] = (
                    atms[site2[2]][0], list(atm_dic[0])[1]-bnd_ord2, None)
                bnd_dct[frozenset({4, 0})] = (bnd_ord2, None)
                bnd_dct[frozenset({2, 3})] = (1, None)
                # Add in the extra cbh1 bond
                key1 = [atma, atmb]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic[5] = (atms[atmb][0], atm_vals[atmb]-bnd_ord1, None)
                if atma == site1[0]:
                    atm_dic[0] = (
                        atm_dic[0][0], atm_dic[0][1] - bnd_ord1, None)
                    bnd_dct[frozenset({0, 5})] = (bnd_ord1, None)
                elif atma == site1[1]:
                    atm_dic[1] = (
                        atm_dic[1][0], atm_dic[1][1] - bnd_ord1, None)
                    bnd_dct[frozenset({1, 5})] = (bnd_ord1, None)
                elif atma == site1[2]:
                    atm_dic[2] = (
                        atm_dic[2][0], atm_dic[2][1] - bnd_ord1, None)
                    bnd_dct[frozenset({2, 5})] = (bnd_ord1, None)
                elif atma == site2[0]:
                    atm_dic[3] = (
                        atm_dic[3][0], atm_dic[3][1] - bnd_ord1, None)
                    bnd_dct[frozenset({3, 5})] = (bnd_ord1, None)
                elif atma == site2[1]:
                    atm_dic[4] = (
                        atm_dic[4][0], atm_dic[4][1] - bnd_ord1, None)
                    bnd_dct[frozenset({4, 5})] = (bnd_ord1, None)
                unique_bond = True
            elif atms[atma][0] != 'H' and atms[atmb][0] != 'H':
                bnd_ord = list(bnd_ords[bnd])[0]
                atm_dic = {
                    0: (atms[atma][0], int(atm_vals[atma])-bnd_ord, None)}
                atm_dic[1] = (
                    atms[atmb][0], int(atm_vals[atmb])-bnd_ord, None)
                bnd_dct = {frozenset({0, 1}): (bnd_ord, None)}
                unique_bond = True
            if unique_bond:
                grai = (atm_dic, bnd_dct)
                try:
                    grai = automol.graph.explicit(grai)
                    key = 'exp_gra'
                except:
                    key = 'ts_gra'
                newname = None
                repeat = False
                for name in frags:
                    if key in frags[name]:
                        if key == 'exp_gra':
                            if automol.graph.full_isomorphism(
                                   frags[name][key], grai):
                                newname = name
                                repeat = True
                        else:
                            if frags[name][key] == grai:
                                newname = name
                                repeat = True
                if not repeat:
                    newname = len(frags.keys())
                    frags[newname] = {}
                    frags[newname][key] = grai
                util.add2dic(frags[newname], 'coeff', coeff)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            zedfrags = cbhzed_elim(gra, site1, site2, bal=False)
            newfrags = frags.copy()
            for zedname in zedfrags:
                newname = None
                repeat = False
                if 'exp_gra' in zedfrags[zedname]:
                    key = 'exp_gra'
                    for onename in newfrags:
                        if 'exp_gra' in newfrags[onename]:
                            if automol.graph.full_isomorphism(
                                    newfrags[onename][key],
                                    zedfrags[zedname][key]):
                                newname = onename
                                repeat = True
                else:
                    key = 'ts_gra'
                if not repeat:
                    newname = len(newfrags.keys())
                    newfrags[newname] = {}
                    newfrags[newname][key] = zedfrags[zedname][key]
                util.add2dic(
                    newfrags[newname], 'coeff',  -zedfrags[zedname]['coeff'])
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
    _, atms, bnd_ords, atm_vals, _ = tsutil.ts_graph(gra, site)

    # Determine CBHone fragments
    frags = {}

    for bnd in bnd_ords:
        atma, atmb = bnd
        if ((atms[atma][0] != 'H' or atma in site) and
           (atms[atmb][0] != 'H' or atmb in site)):
            if atma in site and atmb in site:
                continue
            coeff = 1.0
            if atmb in site:
                atmb, atma = atma, atmb
            if atma in site:
                key1 = [site[0], site[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic = {
                    0: (atms[site[0]][0], atm_vals[site[0]]-bnd_ord1, None)}
                bnd_dct = {frozenset({0, 1}): (bnd_ord1, None)}
                key1 = [site[2], site[1]]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord2 = list(bnd_ords[key])[0]
                atm_dic[2] = (
                    atms[site[2]][0], atm_vals[site[2]]-bnd_ord2, None)
                atm_dic[1] = (
                    atms[site[1]][0],
                    atm_vals[site[1]]-bnd_ord1-bnd_ord2, None)
                bnd_dct[frozenset({1, 2})] = (bnd_ord2, None)
                key1 = [atma, atmb]
                key1.sort()
                key = frozenset({*key1})
                bnd_ord1 = list(bnd_ords[key])[0]
                atm_dic[3] = (atms[atmb][0], atm_vals[atmb]-bnd_ord1, None)
                if atma == site[0]:
                    atm_dic[0] = (
                        atm_dic[0][0], atm_dic[0][1] - bnd_ord1, None)
                    bnd_dct[frozenset({0, 3})] = (bnd_ord1, None)
                elif atma == site[2]:
                    atm_dic[2] = (
                        atm_dic[2][0], atm_dic[2][1] - bnd_ord1, None)
                    bnd_dct[frozenset({2, 3})] = (bnd_ord1, None)
                elif atma == site[1]:
                    atm_dic[1] = (
                        atm_dic[1][0], atm_dic[1][1] - bnd_ord1, None)
                    bnd_dct[frozenset({1, 3})] = (bnd_ord1, None)
            else:
                bnd_ord = list(bnd_ords[bnd])[0]
                atm_dic = {
                    0: (atms[atma][0], int(atm_vals[atma])-bnd_ord, None)}
                atm_dic[1] = (atms[atmb][0], int(atm_vals[atmb])-bnd_ord, None)
                bnd_dct = {frozenset({0, 1}): (bnd_ord, None)}
            grai = (atm_dic, bnd_dct)
            try:
                grai = automol.graph.explicit(grai)
                key = 'exp_gra'
            except:
                key = 'ts_gra'
            newname = None
            repeat = False
            for name in frags:
                if key in frags[name]:
                    if key == 'exp_gra':
                        if automol.graph.full_isomorphism(
                           frags[name][key], grai):
                            newname = name
                            repeat = True
                    else:
                        if frags[name][key] == grai:
                            newname = name
                            repeat = True
            if not repeat:
                newname = len(frags.keys())
                frags[newname] = {}
                frags[newname][key] = grai
            util.add2dic(frags[newname], 'coeff', coeff)
    frags = tsutil.simplify_gra_frags(frags)
    if bal:
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            zedfrags = cbhzed_habs(gra, site, bal=False)
            newfrags = frags.copy()
            for zedname in zedfrags:
                newname = None
                repeat = False
                if 'exp_gra' in zedfrags[zedname]:
                    key = 'exp_gra'
                    for onename in newfrags:
                        if 'exp_gra' in newfrags[onename]:
                            if automol.graph.full_isomorphism(
                               newfrags[onename][key],
                               zedfrags[zedname][key]):
                                newname = onename
                                repeat = True
                else:
                    key = 'ts_gra'
                if not repeat:
                    newname = len(newfrags.keys())
                    newfrags[newname] = {}
                    newfrags[newname][key] = zedfrags[zedname][key]
                util.add2dic(
                    newfrags[newname], 'coeff', -zedfrags[zedname]['coeff'])
            frags = newfrags
        balance_ = util.balance_ts(gra, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags_ts(gra, frags)
    frags = tsutil.simplify_gra_frags(frags)

    return frags
