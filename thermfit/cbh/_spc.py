""" Computes the Heat of Formation at 0 K for a given species

    Explanation of various basis determination schemes
"""

import automol.inchi
import automol.graph
import automol.formula
from thermfit.cbh import _util as util
from thermfit.cbh._basic import basic_basis


# Main callable function
def species_cbh_basis(ich, scheme, balance=True):
    """ Get the basis for the appropriate CBH scheme

        :param ich: InChI string for spc
        :type ich: str
        :param scheme: CBH Scheme used to generate basis
        :type scheme: str
        :param balance:
        :type balance: bool
        :rtype: (tuple(str), list(float))
    """

    if scheme == 'basic':
        frag_lst, coeff_lst = basic_basis(ich)
    else:
        frags = CBH_SCHEMES[scheme](ich, balance=balance)
        frag_lst, coeff_lst = (), ()
        for frag in frags:
            frag_lst += (frag,)
            coeff_lst += (frags[frag],)

    return frag_lst, coeff_lst


# Individual CBH-N calculators
def cbhzed(ich, balance=True):
    """
    Fragments molecule so that each heavy-atom is a seperate fragment
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.inchi.graph(ich)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    atms = automol.graph.atoms(gra)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)

    # Determine CBHzed fragments
    frags = {}
    for atm in atm_vals:
        coeff = 1
        if not balance:
            coeff = (
                util.branch_point(adj_atms[atm]) *
                util.terminal_moiety(adj_atms[atm])
            )
        if atm in rad_atms:
            atm_vals[atm] -= 1
        atm_dic = {0: (atms[atm][0], int(atm_vals[atm]), None)}
        gra = (atm_dic, {})
        frag = automol.graph.inchi(gra)
        util.add2dic(frags, frag, coeff)

    if balance:
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags(ich, frags)

    return frags


def cbhone(ich, balance=True):
    """
    Fragments molecule in a way that conserves each heavy-atom/heavy-atom bond
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments
    and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.inchi.graph(ich)
    atms = automol.graph.atoms(gra)
    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)

    # Determine CBHone fragments
    frags = {}
    for atm in atm_vals:
        for adj in list(adj_atms[atm]):
            if atm > adj:
                vali = atm_vals[atm]
                valj = atm_vals[adj]
                if atm in rad_atms:
                    vali -= 1
                if adj in rad_atms:
                    valj -= 1
                key = frozenset({atm, adj})
                bnd_ord = list(bnd_ords[key])[0]
                vali -= bnd_ord
                valj -= bnd_ord
                atm_dic = {0: (atms[atm][0], int(vali), None),
                           1: (atms[adj][0], int(valj), None)}
                bnd_dic = {frozenset({0, 1}): (1, None)}
                gra = (atm_dic, bnd_dic)
                frag = automol.graph.inchi(gra)
                util.add2dic(frags, frag)
    frags = {k: v for k, v in frags.items() if v}
    if not frags:
        frags = cbhzed(ich)
    # Balance
    if balance:
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            zedfrags = cbhzed(ich, balance=False)
            newfrags = frags.copy()
            for frag in zedfrags:
                util.add2dic(newfrags, frag, -zedfrags[frag])
            frags = {k: v for k, v in newfrags.items() if v}
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            frags = util.balance_frags(ich, frags)
    return frags


def cbhtwo(ich, balance=True):
    """
    Fragments molecule for each heavy-atom to stay bonded to its adjacent atoms
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name for fragments and
    value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.inchi.graph(ich)
    atms = automol.graph.atoms(gra)
    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)

    # Determine CBHtwo fragments
    frags = {}
    for atm in atms:
        vali = atm_vals[atm]
        if atm in rad_atms:
            vali -= 1
        # First loop over all atoms of this frag to get saturation of atomi
        for adj in list(adj_atms[atm]):
            key = frozenset({atm, adj})
            bnd_ord = list(bnd_ords[key])[0]
            vali -= bnd_ord
        atm_dic = {0: (atms[atm][0], int(vali), None)}
        bnd_dic = {}
        # Then start adding bonds to the bnddic and atomdic
        j = 0
        coeff = 1
        if not balance:
            coeff = (
                util.branch_point(adj_atms[atm]) *
                util.terminal_moiety(adj_atms[atm])
            )
        for adj in list(adj_atms[atm]):
            j += 1
            valj = atm_vals[adj]
            if adj in rad_atms:
                valj -= 1
            key = frozenset({atm, adj})
            bnd_ord = list(bnd_ords[key])[0]
            valj -= bnd_ord
            atm_dic[j] = (atms[adj][0], int(valj), None)
            bnd_dic[frozenset({0, j})] = (1, None)
        gra = (atm_dic, bnd_dic)
        frag = automol.graph.inchi(gra)
        util.add2dic(frags, frag, coeff)

    frags = {k: v for k, v in frags.items() if v}
    if not frags:
        frags = cbhone(frags)
    # Balance
    if balance:
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            newfrags = frags.copy()
            onefrags = cbhone(ich, balance=False)
            for frag in onefrags:
                util.add2dic(newfrags, frag, -onefrags[frag])
            frags = {k: v for k, v in newfrags.items() if v}
            balance_ = util.balance(ich, frags)
            balance_ = {k: v for k, v in balance_.items() if v}
            if balance_:
                newfrags = frags.copy()
                zedfrags = cbhzed(ich, balance=False)
                for frag in zedfrags:
                    util.add2dic(newfrags, frag, zedfrags[frag])
                frags = {k: v for k, v in newfrags.items() if v}
                balance_ = util.balance(ich, frags)
                balance_ = {k: v for k, v in balance_.items() if v}
                if balance_:
                    frags = util.balance_frags(ich, frags)

    return frags


def cbhthree(ich, balance=True):
    """
    Fragments molecule to retain each heavy-atom -- heavy-atom bond,
    and keep the bonds of each    atm1 b1 atm2 b2 atm3 b3 atm4 b4 atm5
    would give atm1 [b1] atm2 b2 atm3, atm1 b1 atm2 [b2] atm3 b3 atm4,
    atm2 b2 atm3 [b3] atm4 b4 atm5, and atm3 b3 atm4 [b4] atm5
    INPUT:
    ich --  STR inchi name for molecule
    OUTPUT
    frags -- DIC dictionary with keys as STR inchi name
    for fragments and value as INT their coefficient
    """

    # Graphical info about molecule
    gra = automol.inchi.graph(ich)
    atms = automol.graph.atoms(gra)
    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
    atm_vals = automol.graph.atom_element_valences(gra)
    adj_atms = automol.graph.atoms_neighbor_atom_keys(gra)

    # Determine CBHfour fragments
    frags = {}

    for bnd in list(bnd_ords):
        atm_dic = {}
        bnd_dic = {}
        bnd_dic[frozenset({0, 1})] = (list(bnd_ords[bnd])[0], None)
        for i, atm in enumerate(list(bnd)):
            vali = atm_vals[atm]
            if atm in rad_atms:
                vali -= 1
            for adj in list(adj_atms[atm]):
                key = frozenset({atm, adj})
                bnd_ord = list(bnd_ords[key])[0]
                vali -= bnd_ord
            atm_dic[i] = (atms[atm][0], int(vali), None)
            for j, adj in enumerate(list(adj_atms[atm]), start=1):
                if adj not in list(bnd):
                    valj = atm_vals[adj]
                    if adj in rad_atms:
                        valj -= 1
                    key = frozenset({atm, adj})
                    bnd_ord = list(bnd_ords[key])[0]
                    valj -= bnd_ord
                    atm_dic[i*4+j+1] = (atms[adj][0], int(valj), None)
                    bnd_dic[frozenset({i, i*4+j+1})] = (bnd_ord, None)
        gra = (atm_dic, bnd_dic)
        frag = automol.graph.inchi(gra)
        util.add2dic(frags, frag)

    if not frags:
        frags = cbhtwo(frags)

    if balance:
        balance_ = util.balance(ich, frags)
        balance_ = {k: v for k, v in balance_.items() if v}
        if balance_:
            newfrags = frags.copy()
            twofrags = cbhtwo(ich, balance=False)
            for frag in twofrags:
                util.add2dic(newfrags, frag, -twofrags[frag])
            frags = {k: v for k, v in newfrags.items() if v}
            # balance_ = util.balancec(ich, frags)
            # balance_ = {k: v for k, v in balance_.items() if v}
            # if balance_:
            #    newfrags = frags.copy()
            #    newerfrags = {}
            #    onefrags = cbhone(ich, balance=False)
            #    for frag in onefrags:
            #         util.add2dic(newfrags, frag, - onefrags[frag])
            #    frags = {k: v for k, v in newfrags.items() if v}
            #    balance_ = util.balancec(ich, frags)
            #    balance_ = {k: v for k, v in balance_.items() if v}
            #    if balance_:
            #        newfrags = frags.copy()
            #        zedfrags = cbhzed(ich, balance=False)
            #        for frag in zedfrags:
            #            util.add2dic(newfrags, frag, zedfrags[frag])
            #        frags = {k: v for k, v in newfrags.items() if v}
            #        balance_ = util.balancec(ich, frags)
            #        balance_ = {k: v for k, v in balance_.items() if v}
            #        if balance_:
            #            frags = util.balancec_frags(ich, frags)
    # balance = util.balancec(ich, frags)
    # balance_ = {k: v for k, v in balance_.items() if v}
    return frags

# def cbhfour(ich):
#    """
#    Fragments molecule for each heavy-atom to stay bonded to its
#    adjacent atoms, for those atoms to say bonded to their adjacent atoms
#    atm1 b1 atm2 b2 atm3 b3 atm4 b4 atm5
#    would give [atm1] b1 atm2 b2 atm3, atm1 b1 [atm2] b2 atm3 b3 atm4,
#    atm1 b1 atm2 b2  [atm3] b3 atm4 b4 amt5, atm2 b2 atm3 b3 [atm4] b4 atm5,
#    atm3 b3 atm4 b4 [atm5]
#    INPUT:
#    ich --  STR inchi name for molecule
#    OUTPUT
#    frags -- DIC dictionary with keys as STR inchi name for fragments and
#    value as INT their coefficient
#    """
#    #Graphical info about molecule
#    gra      = automol.inchi.graph(ich)
#    atms     = automol.graph.atoms(gra)
#    bnd_ords = automol.graph.one_resonance_dominant_bond_orders(gra)
#    rad_atms = list(automol.graph.sing_res_dom_radical_atom_keys(gra))
#    atm_vals = automol.graph.atom_element_valences(gra)
#    adj_atms = automol.graph.atom_neighbor_keys(gra)
#
#    #Determine CBHfour fragments
#    frags = {}
#    for atmi in atms:
#        vali = atm_vals[atmi]
#        if atmi in rad_atms:
#            vali -= 1
#        #Atoms j are adj to atoms i
#        #Loop over j one time to get the valence of i
#        for atmj in list(adj_atms[atmi]):
#            key     = frozenset({atmi, atmj})
#            bnd_ord = list(bnd_ords[key] )[0]
#            vali   -= bnd_ord
#        atm_dic = {0: (atms[atmi][0], int(vali), None)}
#        bnd_dic = {}
#        # Then Loop over j a second time to get bonds
#        # to i-j, valence of j, and bonds j-k
#        for j, atmj in enumerate(list(adj_atms[atmi]), start = 1):
#            valj = atm_vals[atmj]
#            if atmj in rad_atms:
#                valj -= 1
#            for atmk in list(adj_atms[atmj]):
#            #Loop over k first to get the valence of j
#                key     = frozenset({atmj, atmk})
#                bnd_ord = list(bnd_ords[key] )[0]
#                valj   -= bnd_ord
#            key     = frozenset({atmi, atmj})
#            bnd_ord = list(bnd_ords[key] )[0]
#
#            #Add i-j bonds and j atoms
#            atm_dic[j] = (atms[atmj][0], int(valj), None)
#            bnd_dic[frozenset({0, j})] =  (1, None)
#
#            #Loop over k to add atoms k and j-k bonds
#            for k, atmk in enumerate(list(adj_atms[atmj]), start = 1):
#               if atmk != atmi:
#                   valk = atm_vals[atmk]
#                   if atmk in rad_atms:
#                       valk -= 1
#                   key     = frozenset({atmj, atmk})
#                   bnd_ord = list(bnd_ords[key] )[0]
#                   valk   -= bnd_ord
#                   index   = k + len(list(adj_atms[atmi]))
#                   atm_dic[index] = (atms[atmk][0], int(valk), None)
#                   bnd_dic[frozenset({j, index})] =  (1, None)
#               else:
#                   k -= 1
#
#        gra     = (atm_dic, bnd_dic)
#        frag = automol.graph.inchi(gra)
#        util.add2dic(frags, frag)
#    frags =  {k: v for k, v in frags.items() if v}
#    ioprinter.info_message(frags)
#    return
#    #Balance
#    return frags


# Dictionary of schemes
CBH_SCHEMES = {
    'cbh0': cbhzed,
    'cbh1': cbhone,
    'cbh1_0': cbhone,
    'cbh1_1': cbhone,
    'cbh2': cbhtwo,
    'cbh2_1': cbhtwo,
    'cbh2_2': cbhtwo,
    'cbh3': cbhthree
}
