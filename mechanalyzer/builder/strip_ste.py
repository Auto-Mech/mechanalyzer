import copy
from mechanalyzer.calculator import compare
from mechanalyzer.calculator import combine
from mechanalyzer.calculator import rates
from mechanalyzer.builder import checker
from mechanalyzer.builder import _names as names
from ratefit.fit import _fit as fit
from automol import inchi


# Step 10
def comb_strpd_and_no_ste(re_mech_spc_dct, re_rxn_param_dct,
                             mech_spc_dct_no_ste, rxn_param_dct_no_ste):


#    comb_rxn_param_dct, comb_spc_nasa7_dct, _ = combine.comb_mechs(
#        re_rxn_param_dct, rxn_param_dct_no_ste,

    re_mech_spc_dct_comb = copy.deepcopy(re_mech_spc_dct)
    re_rxn_param_dct_comb = copy.deepcopy(re_rxn_param_dct)

    re_mech_spc_dct_comb.update(mech_spc_dct_no_ste)
    re_rxn_param_dct_comb.update(rxn_param_dct_no_ste)

    return re_mech_spc_dct_comb, re_rxn_param_dct_comb


# Step 9
def regenerate_names(mech_spc_dct_strpd_ich, rxn_param_dct_strpd):
    """ Takes the dcts in terms of inchis and regenerates names, then renames
    """

    map_dct = names.functional_group_name_dct(mech_spc_dct_strpd_ich)
    print('map_dct:\n', map_dct)
    re_mech_spc_dct, re_rxn_param_dct = names.remap_mechanism_names(
        mech_spc_dct_strpd_ich, rxn_param_dct_strpd, map_dct)

    return re_mech_spc_dct, re_rxn_param_dct


# Step 8
def join_rxns(iso_sets_par_avg):
    """ Takes iso_sets where each rxn has a single set of params and combines
        into a single rxn_param_dct
    """

    rxn_param_dct_strpd = {}
    for iso_set in iso_sets_par_avg:
        for rxn, params in iso_set.items():
            assert rxn not in rxn_param_dct_strpd, (
                f'The rxn {rxn} already exists!')
            rxn_param_dct_strpd[rxn] = params

    return rxn_param_dct_strpd


# Step 7
def get_avg_params(algn_iso_sets_rxns, temps_lst, pressures):

    iso_sets_par_avg = []
    for iso_set in algn_iso_sets_rxns:
        new_iso_set = {}
        for rxn, params_lst in iso_set.items():
            #new_iso_set[rxn] = params_lst[0]  # for now, just take first
            print(f'fitting rxn {rxn}')
            new_iso_set[rxn] = combine_params(params_lst, temps_lst,
                                              pressures)
        iso_sets_par_avg.append(new_iso_set)

    return iso_sets_par_avg


# Step 6
def align_rxns(iso_sets_rxns_ich):
    """ Align the dictionary

        NOTE: assumes that all reactions that should go together are written
        in the same direction! Prints a warning if exceptions are found
    """

    algn_iso_sets_rxns = []
    for iso_set_rxns in iso_sets_rxns_ich:
        algn_iso_set_rxns = compare.align_dcts(iso_set_rxns)
        algn_iso_sets_rxns.append(algn_iso_set_rxns)

    # Check to see if any reactions match within the current iso set
    for iso_set in algn_iso_sets_rxns:
        iso_set_copy = copy.deepcopy(iso_set)
        for rxn in iso_set.keys():
            iso_set_copy.pop(rxn)  # remove current rxn so no false matches
            matching_rxn, rev_rate = compare.assess_rxn_match(
                rxn, iso_set_copy)
            if matching_rxn:
                print(f'matching rxn {rxn} found where it should not be!')
                print('This occurs b/c the same rxn was written w/a diff. '
                      'permutation of rcts/prds or in reverse')
    
    # Check to see if any reactions match between different iso sets
    num_iso_sets = len(algn_iso_sets_rxns)
    # Loop over each iso set
    for iso_set_idx in range(num_iso_sets):
        # Copy the current algn_iso_sets_rxns (changes dynamically!)
        current_lst = copy.deepcopy(algn_iso_sets_rxns)
        # Get all iso_sets after the current one; these will be searched
        srch_lst = current_lst[(iso_set_idx + 1):]  
        # For each rxn in the current iso_set, look for matches
        for curr_rxn in current_lst[iso_set_idx]:
            # Look through each iso set in the search list
            for srch_idx, srch_iso_set in enumerate(srch_lst):
                match, rev = compare.assess_rxn_match(curr_rxn, srch_iso_set)
                # If match is found, add current params list to that in the 
                # original algn_iso_sets_rxns AND remove the rxn from the 
                # original algn_iso_sets_rxns
                if match and not rev:
                    params_lst = srch_iso_set[match]  # list of RxnParams 
                    algn_iso_sets_rxns[iso_set_idx][curr_rxn] += params_lst
                    algn_iso_sets_rxns[iso_set_idx + srch_idx + 1].pop(match)

    return algn_iso_sets_rxns


# Step 5
def rename_iso_sets_rxns(iso_sets_rxns, mech_spc_dct_strpd_ich,
                         mech_spc_dct_strpd):
    """ Renames all isomers according to their (stereo-stripped) inchis
    """

    # Get the rename instructions
    rename_instr = compare.get_rename_instr_v2(mech_spc_dct_strpd_ich,
                                               mech_spc_dct_strpd)
    # Rename each set of iso_rxns
    iso_sets_rxns_ich = []
    for iso_set_rxns in iso_sets_rxns:
        iso_set_rxns_ich = []
        for iso_rxns in iso_set_rxns:
            # Rename all reactions containing this species and store
            iso_rxns_ich, _ = compare.rename_species(iso_rxns, rename_instr)
            #print('iso_rxns_ich, in loop:\n', iso_rxns_ich)
            iso_set_rxns_ich.append(iso_rxns_ich)
        iso_sets_rxns_ich.append(iso_set_rxns_ich)

    return iso_sets_rxns_ich


# Step 4b
def get_no_ste_rxns(rxn_param_dct, iso_sets_rxns):
    """ Gets a rxn_param_dct with all reactions that have no stereoisomers
    """

    # Get list of all rxns involving stereoisomers
    all_iso_rxns = []
    for iso_set_rxns in iso_sets_rxns:
        for iso_rxns in iso_set_rxns:
            all_iso_rxns.extend(list(iso_rxns.keys()))

    rxn_param_dct_no_ste = copy.deepcopy(rxn_param_dct)
    for rxn in all_iso_rxns:
        if rxn in rxn_param_dct_no_ste:
            rxn_param_dct_no_ste.pop(rxn)

    return rxn_param_dct_no_ste


# Step 4
def get_ste_rxns(rxn_param_dct, iso_sets):
    """ For each isomer, gets a rxn param dct with all rxns containing that
        isomer
    """
    def search_rxns(rxn_param_dct, spc):
        """ Searches a rxn_param_dct for all rxns that contain a species
        """
        filt_rxn_param_dct = {}
        for rxn, params in rxn_param_dct.items():
            rcts, prds, _ = rxn
            if spc in rcts or spc in prds:
                filt_rxn_param_dct[rxn] = params

        return filt_rxn_param_dct

    iso_sets_rxns = []
    # Loop over each set of isomers
    for iso_set in iso_sets:
        iso_set_rxns = []
        for iso in iso_set:
            # Get all rxns that contain this isomer and store
            iso_rxn_params = search_rxns(rxn_param_dct, iso)
            iso_set_rxns.append(iso_rxn_params)
        iso_sets_rxns.append(iso_set_rxns)

    return iso_sets_rxns


# Step 3
def make_mech_spc_dct_ich(iso_sets, mech_spc_dct_strpd):
    """ Renames a mech_spc_dct_strpd (i.e., all stereo stripped) to have the
        species names be inchis
    """

    mech_spc_dct_strpd_ich = {}
    for iso_set in iso_sets:
        iso = iso_set[0]  # just take first iso b/c all spc_dcts are same
        spc_dct = mech_spc_dct_strpd[iso]
        ich = spc_dct['inchi']
        mech_spc_dct_strpd_ich[ich] = spc_dct

    return mech_spc_dct_strpd_ich


# Step 2
def find_iso_sets(mech_spc_dct_strpd):
    """ Finds all sets of isomers in a stripped mech_spc_dct that are now the
        exact same species (since stereo has been stripped)

        :return iso_sets: list of all sets of stereoisomers
        :type iso_sets: [[iso1, iso2, ...], [iso1, iso2, ...], ...]
    """

    # Get two things: (i) species and (ii) inchis
    spcs = tuple(mech_spc_dct_strpd.keys())
    ichs = ()
    for spc_dct in mech_spc_dct_strpd.values():
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
                same = are_spc_dcts_same(ref_spc_dct, spc_dct)
                assert same, (f'In the set of stereoisomes {iso_set}, the '
                               'species dcts are not the same (even after '
                               'stripping stereo)')

    return iso_sets


# Step 1
def strip_mech_spc_dct(mech_spc_dct):
    """ Removes stereochemistry from all species in a mech_spc_dct. Returns
        a new mech_spc_dct with all the stereo-specific species, but now
        stripped of the stereo. Also returns a separate mech_spc_dct of
        species that originally had no stereo in them (unchanged)

        :param mech_spc_dct: input mech_spc_dct
        :type mech_spc_dct: dict
        :return mech_spc_dct_strpd: dct with only stereo specific spcs, but
            with stereo stripped from the inchis
        :type mech_spc_dct_strpd: dict
        :return mech_spc_dct_no_ste: dct with all non-stereo spcs, unchanged
        :type mech_spc_dct_no_ste: dict
    """

    mech_spc_dct_strpd = {}
    mech_spc_dct_no_ste = {}
    for spc, spc_dct in mech_spc_dct.items():
        orig_ich = copy.copy(spc_dct['inchi'])
        strpd_ich = inchi.without_stereo(orig_ich)
        # If the species is stereo-free, save in the no_ste dct
        if orig_ich == strpd_ich:
            mech_spc_dct_no_ste[spc] = spc_dct
        # If the species had stereo, save stereo-stripped info in strpd dct
        else:
            # Get the smiles and inchikey without stereo
            strpd_smi = inchi.smiles(strpd_ich)
            strpd_ichkey = inchi.inchi_key(strpd_ich)
            # Store the stereo-stripped information
            spc_dct['smiles'] = strpd_smi
            spc_dct['inchikey'] = strpd_ichkey
            spc_dct['inchi'] = strpd_ich
            mech_spc_dct_strpd[spc] = spc_dct

    return mech_spc_dct_strpd, mech_spc_dct_no_ste


def are_spc_dcts_same(spc_dct1, spc_dct2):
    """ Checks if two spc dcts are the same
    """

    ich1 = spc_dct1['inchi']
    mlt1 = spc_dct1['mult']
    chg1 = spc_dct1['charge']
    exc1 = spc_dct1['exc_flag']
    fml1 = spc_dct1['fml']

    same = compare.are_spc_same(ich1, mlt1, chg1, exc1, fml1, spc_dct2)

    return same


def combine_params(params_lst, temps_lst, pressures, method='avg'):
    """ Takes a list of rate parameters, calculates the k(T,P) values, either
        averages or adds them, and then refits them
    """

    # Catch for when no params are given
    if len(params_lst) == 0:
        return None

    # Calculate rates for each params in the list
    ktp_dct_lst = []
    for params in params_lst:
        if params is not None:
            # Get rates
            ktp_dct = rates.eval_params(params, temps_lst, pressures)
            ktp_dct_lst.append(ktp_dct)
            # Get fit form (used later)
            fit_method = params.get_existing_forms()[0]  # first one

    # Average rate constants
    if method == 'avg':
        for idx, ktp_dct in enumerate(ktp_dct_lst):
            if idx == 0:
                summed_ktp_dct = copy.deepcopy(ktp_dct)
            else:
                summed_ktp_dct = rates.add_ktp_dcts(summed_ktp_dct, ktp_dct)
        factor = 1 / len(params_lst)
        comb_ktp_dct = rates.mult_by_factor(summed_ktp_dct, factor)

    # Fit the averaged rate constants
#    fit_methods = params_lst[0].get_existing_forms()
#    fit_method = fit_methods[0]
    comb_params, _ = fit.fit_ktp_dct(comb_ktp_dct, fit_method)

    return comb_params
