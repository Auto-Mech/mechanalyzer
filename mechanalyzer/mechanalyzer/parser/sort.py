"""
Sorter module - sorting of the mechanism according to various options
- reaction classes
- pes/subpes
- multiplicity
- species subsets
- submechanism (calling submech module)
"""

import sys
import copy
import pandas as pd
import numpy as np
import automol
from mechanalyzer.parser import ckin_ as ckin
from mechanalyzer.parser import submech
from mechanalyzer.parser import util
from mechanalyzer.parser import pes


class SortMech:
    '''
    class of functions to organize the mechanism according to given criteria
    from any step after initialization: call "return_mech_df"
    to get the current dataframe with mech info
    '''

    def __init__(self, mech_info, spc_dct):
        # extract data from mech info
        [_, formula_str_lst, rct_names_lst,
            prd_names_lst, rxn_name_lst] = mech_info
        # set dataframe: extract useful info
        molecularity = list(map(len, rct_names_lst))
        n_of_prods = list(map(len, prd_names_lst))
        rct_names_lst_ordered = util.order_rct_bystoich(
            rct_names_lst, spc_dct=spc_dct)  # put heavier reactant first
        prd_names_lst_ordered = util.order_rct_bystoich(
            prd_names_lst, spc_dct=spc_dct)  # put heavier product first
        rct_1, rct_2 = util.get_S1S2(rct_names_lst_ordered)
        data = np.array([rct_names_lst, prd_names_lst, rct_names_lst_ordered, prd_names_lst_ordered,
                         rct_1, rct_2, molecularity, n_of_prods, formula_str_lst], dtype=object).T
        self.mech_df = pd.DataFrame(data, index=rxn_name_lst, columns=[
                                    'rct_names_lst', 'prd_names_lst', 'rct_names_lst_ord',
                                    'prd_names_lst_ord', 'R1', 'R2', 'molecularity',
                                    'N_of_prods', 'pes'])

        self.spc_dct = spc_dct  # set for later use
        # empty list for initialization (otherwise pylint warning)
        self.species_list = []

    def sort(self, hierarchy, species_list):
        '''
        hierarchy = list of hierarchical criteria for the mech organization
        sorts the mechanism according to the given criteria
        species_list = list of species you want to isolate (empty to process the full mech).
        For extraction of a fuel submech:  ['speciesname','SUBMECH']
        '''

        # if species_list is not empty: pre-process the mechanism
        if len(species_list) > 0:
            if len(species_list) == 2 and species_list[-1] == 'submech':
                # select a subset of species appropriately according to the stoichiometries
                # specified in submech.py
                species_list, species_subset_df = submech.species_subset(
                    species_list[0], self.spc_dct)
                self.species_subset_df = species_subset_df
            elif len(species_list) > 2 and species_list[-1] == 'submech':
                print('Error: pyr/ox submech extraction available for only 1 species')
                sys.exit()

            self.mech_df, self.spc_dct = self.filter_byspecies(species_list)
            self.species_list = species_list

        elif len(species_list) == 0 and 'species' in hierarchy:
            print('*ERROR: sorting by species incompatible with empty isolate_species')
            sys.exit()

        # 0. look for keywords and sort accordingly
        # LIST OF AVAILABLE SORTING OPTIONS (BESIDES R1 AND PES, ALWAYS AVAILABLE)
        sort_optns_dct = {
            'species': self.group_species,
            'subpes': self.conn_chn,
            'mult': self.reac_mult,
            'rxn_class_broad': self.rxnclass_broad,
            'rxn_class_graph': self.rxnclass_graph,
            'submech': self.group_submech,
        }

        for optn, fun_name in sort_optns_dct.items():
            # call sorting function
            if any(optn == inp_crt for inp_crt in hierarchy):
                # generate corresponding dataframe
                df_optn = pd.DataFrame(
                    index=self.mech_df.index, columns=[optn])
                df_optn = fun_name(df_optn)
                # concaFtenate the new portion of dataframe
                self.mech_df = pd.concat([self.mech_df, df_optn], axis=1)

        # 1. sort
        try:
            self.mech_df = self.mech_df.sort_values(by=hierarchy[:-1])
        except KeyError as err:
            print(
                'WARNING: Reactions not sorted according to all criteria: missing {}'.format(err))
            sys.exit()

        # 2. assign class headers
        # set labels for all the possible criteria
        criteria_all = ['molecularity', 'N_of_prods', 'species', 'pes', 'subpes',
                        'submech', 'r1', 'mult', 'rxn_class_broad', 'rxn_class_graph']
        labels_all = ['NR', 'N_of_prods', 'SPECIES', 'PES', 'SUBPES', 'SUBMECH',
                      'Heavier rct', 'Total multiplicity', 'rxn type', 'rxn class']
        labels = pd.Series(labels_all, index=criteria_all)
        self.class_headers(hierarchy, labels)

    def filter_byspecies(self, species_list):
        """
        Find all reactions involving species of the species_list given as input
        make another spc_dct containing only the species appearing in the reactions of interest
        """
        mech_df = copy.deepcopy(self.mech_df)
        # check that all species selected are in the species dictionary
        if any(i not in self.spc_dct.keys() for i in species_list):
            print('Error in ISOLATE_SPECIES: not all species are in the species list ')
            sys.exit()

        spc_list = []
        # for all reactions in the dataframe: check if you have the species of the selected list.
        # Otherwise remove the reaction
        for rxn in mech_df.index:
            rcts = list(mech_df['rct_names_lst'][rxn])
            prds = list(mech_df['prd_names_lst'][rxn])
            # check if one of the species in the list is among reactants or products
            # of the reaction considered
            if (not any(rct == species for species in species_list for rct in rcts)
                    and not any(prd == species for species in species_list for prd in prds)):
                mech_df = mech_df.drop([rxn])
            else:
                # append all species to the list
                spc_list.extend(rcts)
                spc_list.extend(prds)

        # filter spc_list: unique elements
        spc_list = list(set(spc_list))
        # new spc_dct
        spc_dct_val = list(map(self.spc_dct.get, spc_list))
        spc_dct = dict(zip(spc_list, spc_dct_val))

        return mech_df, spc_dct

    def conn_chn(self, conn_chn_df):
        '''
        Identify connected channels
        Generate column 'subpes' in conn_chn_df
        '''

        for fml, peslist in self.mech_df.groupby('pes'):
            # print(peslist)
            # Set the names lists for the rxns and species needed below
            pes_rct_names_lst = peslist['rct_names_lst'].values
            pes_prd_names_lst = peslist['prd_names_lst'].values
            pes_rxn_name_lst = peslist.index
            connchnls = pes.find_conn_chnls(
                pes_rct_names_lst, pes_prd_names_lst, pes_rxn_name_lst)
            # write subpes in conn_chn_df
            for key, value in connchnls.items():
                rxns = peslist.iloc[value].index
                conn_chn_df['subpes'][rxns] = fml + '-' + str(key)

        return conn_chn_df

    def group_species(self, reac_sp_df):
        """
        Creates a new df column "species" - recognizes the species you set in the list
        Also in this case the species_list is hierarchical:
        if the first group contains also a species of the second group,
        the reaction remains in the first group
        """
        # if species list is not found: do nothing - species entry will remain empty
        if len(self.species_list) > 0:
            for rxn in reac_sp_df.index:
                rcts = list(self.mech_df['rct_names_lst'][rxn])
                prds = list(self.mech_df['prd_names_lst'][rxn])
                # check species hierarchically
                for sp_i in self.species_list:
                    if ((any(sp_i == rct for rct in rcts) or any(sp_i == prd for prd in prds))
                            and isinstance(reac_sp_df['species'][rxn], float)):
                        reac_sp_df['species'][rxn] = sp_i
        return reac_sp_df

    def group_submech(self, submech_df):
        """
        Creates a new df column "submech":
        Assigns a submechanism (fuel, fuel radical, fuel add ..) to the reaction
        according to the species taking part to it.
        The species classification is done according to self.species_subset_df
        THIS IS JUST A FIRST DRAFT
        """
        for rxn in submech_df.index:
            rcts = list(self.mech_df['rct_names_lst'][rxn])
            prds = list(self.mech_df['prd_names_lst'][rxn])
            # check species hierarchically (hierarchy fixed in species list)
            for sp_i in self.species_list:
                if ((any(sp_i == rct for rct in rcts) or any(sp_i == prd for prd in prds))
                        and isinstance(submech_df['submech'][rxn], float)):
                    submech_df['submech'][rxn] = self.species_subset_df[sp_i]
        return submech_df

    def reac_mult(self, reac_mult_df):
        '''
        Identify reaction multiplicity from spc_dct
        update column 'mult' in reac_mult_df
        '''
        # assign multiplicity values to each reactant
        for rxn in reac_mult_df.index:
            mult = 1
            rcts = self.mech_df['rct_names_lst'][rxn]
            mult = util.get_mult(rcts, self.spc_dct)
            reac_mult_df['mult'][rxn] = str(mult)

        return reac_mult_df

    def rxnclass_broad(self, rxncl_broad_df):
        '''
        Identify the reaction class by broad classification
        No use of the graph approach
        '''
        for rxn in rxncl_broad_df.index:
            rcts = self.mech_df['rct_names_lst_ord'][rxn]
            prds = self.mech_df['prd_names_lst_ord'][rxn]
            if len(rcts) > 1:
                # bimolecular reaction classification
                rxn_class_broad = submech.classify_bimol(
                    rcts, prds, self.spc_dct)
            else:
                # unimolecular reaction classification
                rxn_class_broad = submech.classify_unimol(
                    rcts, prds, self.spc_dct)
            rxncl_broad_df['rxn_class_broad'][rxn] = rxn_class_broad

        return rxncl_broad_df

    def rxnclass_graph(self, rxncl_graph_df):
        '''
        Identify the reaction class by the graph approach
        1. group by subpes
        2. identify rxn in each subpes
        3. classify well skipping channels
        '''
        # 1. group by subpes
        # check that SUBPES is present in indexes, otherwise generate corresponding dataframe
        if 'subpes' not in self.mech_df.columns:
            df_optn = pd.DataFrame(
                index=self.mech_df.index, columns=['subpes'])
            df_optn = self.conn_chn(df_optn)
            # concatenate the new portion of dataframe
            self.mech_df = pd.concat([self.mech_df, df_optn], axis=1)

        # 2. graph classification or each subpes
        for _, subpes_df in self.mech_df.groupby('subpes'):
            # sort by molecularity: analyze first unimolecular isomerizations,
            # unimolecular decompositions, and then bimolecular reactions
            subpes_df = subpes_df.sort_values(
                by=['molecularity', 'N_of_prods'])
            # REFER TO REORDERED SPECIES NAMES, OTHERWISE YOU MAY HAVE INCONSISTENT SPECIES NAMING
            # subpes species list
            species_subpes = list(set(list(
                subpes_df['rct_names_lst_ord'].values)+list(subpes_df['prd_names_lst_ord'].values)))

            # build dataframe: elementary reactivity matrix
            elem_reac_df = pd.DataFrame(np.zeros((len(species_subpes), len(
                species_subpes)), dtype='<U32'), index=species_subpes, columns=species_subpes)

            # graph classification
            for rxn in subpes_df.index:
                rct_names = subpes_df['rct_names_lst'][rxn]
                prd_names = subpes_df['prd_names_lst'][rxn]
                rct_names_ord = subpes_df['rct_names_lst_ord'][rxn]
                prd_names_ord = subpes_df['prd_names_lst_ord'][rxn]
                # exclude all reactions with more than 2 reactants or products (not elementary!)
                if len(rct_names) < 3 and len(prd_names) < 3:

                    rclass = classify_graph(self.spc_dct, rct_names, prd_names)

                    if rclass is None:
                        if (subpes_df['molecularity'][rxn] == 1
                                and subpes_df['N_of_prods'][rxn] == 1):
                            rclass = 'isomerization'
                            # TEMPORARY - BEFORE ALL ISOMERIZATION TYPES ARE ANALYZED
                        else:
                            rclass = 'unclassified'
                else:
                    rclass = 'unclassified - lumped'

                rxncl_graph_df['rxn_class_graph'][rxn] = rclass

                # store values in the elementary reactivity matrix
                # (for now contaminated with isomerizations)

                elem_reac_df[rct_names_ord][prd_names_ord] = rclass
                elem_reac_df[prd_names_ord][rct_names_ord] = rclass

            # 3. classify well skipping channels
            # reclassify the unclassified reactions A->B+C, B+C->D, B+C->E+F
            for rxn in subpes_df.index:
                if rxncl_graph_df['rxn_class_graph'][rxn] == 'unclassified':

                    # call external function for WS channel classification

                    rxn_type_ws = classify_ws(
                        subpes_df, elem_reac_df, species_subpes, rxn)
                    if rxn_type_ws is not None:
                        rxncl_graph_df['rxn_class_graph'][rxn] = rxn_type_ws

        return rxncl_graph_df

    ################## ASSIGN HEADERS ###################################################

    def class_headers(self, hierarchy, labels):
        """
        Read the hierarchy;
        assign classes based on the selected hierarchy
        adds columns in the df with appropriate comments/headers
        first N keywords --> main class --> comments_top (new col in df)
        other keywords --> other subclasses --> comments_inline
        """
        # df for comments_top and comments_inline
        ept_df = np.zeros((len(self.mech_df.index), 1), dtype=str)
        df_cmts_top = pd.DataFrame(
            ept_df, index=self.mech_df.index, columns=['cmts_top'])
        df_cmts_inline = pd.DataFrame(
            ept_df, index=self.mech_df.index, columns=['cmts_inline'])

        try:
            n_headers = int(hierarchy[-1])
        except ValueError:
            print(
                '*ERROR: Last line of sorting options is the N of criteria to use for class headers')
            sys.exit()
        ####### write topheader comments ######
        if n_headers > 0:
            for name, rxndf in self.mech_df.groupby(hierarchy[:n_headers]):
                # write rxn class as top header comments
                rxnclass = cmts_string(
                    name, labels[hierarchy[:n_headers]], 'class_head')
                idx0 = rxndf.index[0]
                df_cmts_top['cmts_top'][idx0] = rxnclass

                ##### write inline comments if necessary #########
                if n_headers < len(hierarchy)-1:
                    for name2, rxndf2 in rxndf.groupby(hierarchy[n_headers:-1]):
                        rxnclass = cmts_string(
                            name2, labels[hierarchy[n_headers:-1]], 'subclass')
                        df_cmts_inline['cmts_inline'][rxndf2.index] = rxnclass
        else:
            # write only inline comments
            for name, rxndf in self.mech_df.groupby(hierarchy[n_headers:-1]):
                rxnclass = cmts_string(
                    name, labels[hierarchy[n_headers:-1]], 'class')
                df_cmts_inline['cmts_inline'][rxndf.index] = rxnclass

        # concatenate DFs
        self.mech_df = pd.concat(
            [self.mech_df, df_cmts_top, df_cmts_inline], axis=1)
    ###################################### output dataframe ##############################

    def return_mech_df(self):
        '''
        Returns:
        - new_idx: reactants and products indices as tuples
        - cmts: dictionary containing comments of the corresponding reactions; indices are new_idx
        - self.spc_dct: species dictionary; may be different from the input
            if a subset of reactions is selected; useful for mech writing
        '''
        rct_names = self.mech_df['rct_names_lst'].values
        prd_names = self.mech_df['prd_names_lst'].values
        new_idx = list(zip(rct_names, prd_names))
        # store comments in dct
        cmts_df = pd.DataFrame(self.mech_df[['cmts_top', 'cmts_inline']].values, index=new_idx,
                               columns=['cmts_top', 'cmts_inline'])
        cmts = cmts_df.to_dict('index')

        return new_idx, cmts, self.spc_dct


######### functions specific for the sorter - non specific functions are in util ###########

########## functions for rxn graph classification ####################
def classify_graph(spc_dct, rct_names, prd_names):
    '''
    Classification of reactions from reactants and products names
    Requires species dictionary for inchis derivation
    '''
    if len(prd_names) >= 3:
        rclass = 'unclassified - lumped'
    else:
        # Get the inchis and graphs
        rct_ichs = list(spc_dct[rct]['inchi']
                        for rct in rct_names)
        rct_graph = list(map(automol.inchi.graph, rct_ichs))
        #rct_gras = list(
        #    map(automol.graph.without_stereo_parities, rct_graph))
        rct_gras = list(
            map(automol.graph._graph.explicit, rct_graph))
        # print(automol.graph.string(rct_gra))
        prd_ichs = list(spc_dct[prd]['inchi']
                        for prd in prd_names)
        prd_graph = list(map(automol.inchi.graph, prd_ichs))
        #prd_gras = list(
        #    map(automol.graph.without_stereo_parities, prd_graph))
        prd_gras = list(
            map(automol.graph._graph.explicit, prd_graph))
        # ID reaction
        rct_fmls = list(spc_dct[rct]['fml']
                        for rct in rct_names)
        prd_fmls = list(spc_dct[prd]['fml']
                        for prd in prd_names)
        print(rct_gras, prd_gras)
        print(rct_fmls,prd_fmls)
        if automol.formula.reac.is_valid_reaction(rct_fmls, prd_fmls):
            rclass = automol.reac._find.find(
                rct_gras, prd_gras)
            print(rclass)
        else:
            rclass = 'unclassified - Wrong Stoichiometry'
        # check stereo compatibility - I am not sure about this input
        # ret = automol.graph.trans.is_stereo_compatible(rclass, rct_graph, prd_graph)
    return rclass


def classify_ws(subpes_df, elem_reac_df, species_subpes, rxn):
    '''
    Classification of well skipping reaction channels of a given subpes
    Returns the wellskipping reaction type
    WARNING: STILL UNDER CONSTRUCTION - SOME TEMPORARY FEATURES
    '''
    # derive unimolecular species list
    mult_species_subpes = pd.Series(
        list(map(len, species_subpes)), index=species_subpes)
    unimol_species = mult_species_subpes[mult_species_subpes == 1].index

    rct_names = subpes_df['rct_names_lst_ord'][rxn]
    prd_names = subpes_df['prd_names_lst_ord'][rxn]
    # reactants: if bimolecular, find the label of the elementary reaction
    # going to unimolecular species; if unimol, label is 'isom'
    # isolate A+B->C and C->A+B connections
    rxn_types = elem_reac_df[rct_names][unimol_species]
    rxn_types = rxn_types[rxn_types != 'unclassified']
    rxn_types_1 = rxn_types[rxn_types != '']

    rxn_types = elem_reac_df[prd_names][unimol_species]
    rxn_types = rxn_types[rxn_types != 'unclassified']
    rxn_types_2 = rxn_types[rxn_types != '']

    try:
        # TEMPORARY: SHOULD RECONSTRUCT THE FULL PATH FROM REACTANTS TO PRODUCTS
        rxn_type_1 = rxn_types_1[0]
        # TEMPORARY: SHOULD RECONSTRUCT THE FULL PATH FROM REACTANTS TO PRODUCTS
        rxn_type_2 = rxn_types_2[0]

        # WRITE THE REACTION TYPE STRING
        rxn_type_ws = rxn_type_1 + '-' + rxn_type_2 + ' (WS)'
        return rxn_type_ws

    except IndexError:
        return None

######### functions for comments - called by the sorter ###############


def cmts_string(name, label, cltype):
    '''
    Return appropriate comment string depending on the type
    name: rxn class
    cltype: class types. options: class_head, class, subclass
    label: class label
    '''
    # assing top headers:
    tophead = '!!!!!!!!! class !!!!!!!!!\n'
    bottomhead = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
    if isinstance(name, str):
        name = [name]
    elif isinstance(name, int):
        name = [str(name)]
    else:
        name = np.array(name, dtype=str)

    if cltype == 'class_head':
        cmtlabel = '!       '+'_'.join(label)
        cmtlabel += '!       '+'_'.join(name)+'\n'
        rxnclass = tophead + cmtlabel + bottomhead
    else:
        cmtlabel = cltype + ': ' + ' _'.join(label) + '  '
        rxnclass = '! ' + cmtlabel + ' _'.join(name)

    return rxnclass
