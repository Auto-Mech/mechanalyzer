import copy
from mechanalyzer.calculator import compare
from mechanalyzer.calculator import rates
from mechanalyzer.calculator.rates import check_p_t
from mechanalyzer.builder import _names as names
from ratefit.fit import _fit as fit
from automol import inchi
from automol import chi


def lump(rxns_by_idx):

    for rxn_by_idx, [rxns, params_lst] in rxns_by_idx.items():
        lump_dct = {}
        for idx, rxn in enumerate(rxns):  # rxns in terms of spc names
            rcts, prds, _ = rxn
            params = params_lst[idx]
            sort_rcts = tuple(sorted(rcts))  # alphabetize
            if sort_rcts in lump_dct:
                lump_dct[sort_rcts].append(params)
            else:
                lump_dct[sort_rcts] = [params]
         
        for lumped_rxn 
        print('lump_dct: ', lump_dct)                
       
    

def get_rxns_by_idx(rac_sets, rxn_param_dct):

    def get_spc_idxs(rcts_or_prds, rac_sets):
        spc_idxs = []
        for rct_or_prd in rcts_or_prds:
            for rac_idx, rac_set in enumerate(rac_sets):
                if rct_or_prd in rac_set:
                    spc_idxs.append(rac_idx)
                    break
        spc_idxs = sorted(spc_idxs)
        spc_idxs = tuple(spc_idxs)
    
        return spc_idxs
    
    rxns_by_idx = {}
    for rxn, params in rxn_param_dct.items():
        rcts, prds, thrd_bod = rxn
        rct_idxs = get_spc_idxs(rcts, rac_sets)
        prd_idxs = get_spc_idxs(prds, rac_sets)
        rxn_by_idx = (rct_idxs, prd_idxs, thrd_bod)
        if rxn_by_idx in rxns_by_idx:
            rxns_by_idx[rxn_by_idx][0].append(rxn)
            rxns_by_idx[rxn_by_idx][1].append(params)
        else:
            rxns_by_idx[rxn_by_idx] = [[rxn], [params]]

    return rxns_by_idx
            
   
def get_rac_sets(iso_sets, mech_spc_dct):

    rac_sets = []
    for iso_set in iso_sets:
        rac_dct = {}
        for iso in iso_set:
            orig_ich = mech_spc_dct[iso]['canon_enant_ich']
            rac_ich = chi.racemic(orig_ich)
            #if rac_ich == orig_ich and num_chiral_sites >= 1:  # need numchiral
            #    print('unchanged by racemization: ', iso)
            if rac_ich in rac_dct:
                rac_dct[rac_ich].append(iso)
            else:
                rac_dct[rac_ich] = [iso]
        for rac_set in rac_dct.values():
            rac_sets.append(rac_set)

    return rac_sets


def find_iso_sets(mech_spc_dct_strpd, canon_ent=False):
    """ Finds all sets of isomers in a stripped mech_spc_dct that are now the
        exact same species (since stereo has been stripped)

        :param mech_spc_dct_strpd: dct with only stereo specific spcs, but
            with stereo stripped from the inchis
        :type mech_spc_dct_strpd: dict
        :return iso_sets: list of all sets of stereoisomers
        :rtype: [[iso1, iso2, ...], [iso1, iso2, ...], ...]
    """

    # Get two things: (i) species and (ii) inchis
    spcs = tuple(mech_spc_dct_strpd.keys())
    ichs = ()
    for spc_dct in mech_spc_dct_strpd.values():
        if canon_ent:
            ichs += (spc_dct['canon_enant_ich'],)
        else:
            ichs += (spc_dct['inchi'],)

    # Get stereo sets, which are sets of species that are the same if one
    # ignores stereo (usually will be singles or doubles)
    iso_sets = []
    already_done = []
    for idx, ich in enumerate(ichs):
        # If this inchi has already been done, skip it
        if idx in already_done:
            continue
        # Otherwise, store ich and look through remaining ichs for match(es)
        iso_set = [spcs[idx]]
        for sub_idx in range(idx + 1, len(ichs)):
            curr_ich = ichs[sub_idx]
            if ich == curr_ich:
                iso_set.append(spcs[sub_idx])
                already_done.append(sub_idx)  # store that ich has been done
        iso_sets.append(iso_set)

    # Check that each iso_set has identical spc_dcts for all isomers
    for iso_set in iso_sets:
        for iso_idx, iso in enumerate(iso_set):
            if iso_idx == 0:  # if on first iso, store ref_spc_dct
                ref_spc_dct = mech_spc_dct_strpd[iso]
            else:  # if on later isos, check against reference
                spc_dct = mech_spc_dct_strpd[iso]
                same = _are_spc_dcts_same(ref_spc_dct, spc_dct)
                assert same, (f'In the set of stereoisomes {iso_set}, the '
                               'species dcts are not the same (even after '
                               'stripping stereo)')

    return iso_sets


def _are_spc_dcts_same(spc_dct1, spc_dct2):
    """ Checks if two spc dcts are the same
    """

    ich1 = spc_dct1['inchi']
    mlt1 = spc_dct1['mult']
    chg1 = spc_dct1['charge']
    exc1 = spc_dct1['exc_flag']
    fml1 = spc_dct1['fml']

    same = compare.are_spc_same(ich1, mlt1, chg1, exc1, fml1, spc_dct2)

    return same
