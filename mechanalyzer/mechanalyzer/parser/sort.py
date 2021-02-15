"""
Sorter module - sorting of the mechanism according to various options
- reaction classes
- pes/subpes
- multiplicity
- species subsets
- submechanism
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
    """ class of methods to organize the mechanism according to given criteria
    """

    def __init__(self, mech_info, spc_dct):
        """ Initializes the mechanism dataframe and the species dictionary

        :param mech_info: list of mechanism info [formula dct, formulas,
                            rct_names, prd_names, rxn_names]
        :param spc_dct: species dictionary

        :returns: None, updates self.
                    self.mech_df: dataframe with mech info
                    self.spc_dct: species dictionary
        """
        # extract data from mech info
        [formula_dct_lst, formula_str_lst, rct_names_lst,
            prd_names_lst, rxn_name_lst] = mech_info
        # set dataframe: extract useful info
        pes_lst = util.count_atoms(formula_dct_lst)
        molecularity = list(map(len, rct_names_lst))
        n_of_prods = list(map(len, prd_names_lst))
        rct_names_lst_ordered = util.order_rct_bystoich(
            rct_names_lst, spc_dct=spc_dct)  # put heavier reactant first
        prd_names_lst_ordered = util.order_rct_bystoich(
            prd_names_lst, spc_dct=spc_dct)  # put heavier product first
        rct_1, rct_2 = util.get_S1S2(rct_names_lst_ordered)
        data = np.array([rct_names_lst, prd_names_lst, rct_names_lst_ordered, prd_names_lst_ordered,
                         rct_1, rct_2, molecularity, n_of_prods, pes_lst], dtype=object).T
        self.mech_df = pd.DataFrame(data, index=rxn_name_lst, columns=[
                                    'rct_names_lst', 'prd_names_lst', 'rct_names_lst_ord',
                                    'prd_names_lst_ord', 'R1', 'R2', 'molecularity',
                                    'N_of_prods', 'pes'])

        self.spc_dct = spc_dct  # set for later use
        # empty list for initialization (otherwise pylint warning)
        self.species_list = []

    def sort(self, hierarchy, species_list):
        """ Main flow of the sorter: takes a set of reactions and classifies
            them as indicated in hierarchy; possibly filters the species
            according to the species list

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param hierarchy: list of hierarchical criteria for mech organization
                            ['..','..','N']: N is the N of criteria to use for headers
                            all the other criteria are written as inline comments
        :param species_list: list of species subset; may be empty
                            if ['onename','submech']: extracts fuel submechanism
        :returns: None. updates self
        """

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
        labels_all = ['NR', 'N_of_prods', 'SPECIES', 'N_COH_PES', 'N_COH.subpes', 'SUBMECH',
                      'Heavier rct', 'Total multiplicity', 'rxn type', 'rxn class']
        labels = pd.Series(labels_all, index=criteria_all)
        self.class_headers(hierarchy, labels)

    def filter_byspecies(self, species_list):
        """ Find all reactions involving species of the species_list given as input
            Provides a new mechanism of the reactions of the selected subset
            and a reduced list of species

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param species_list: list of species subsets

        :returns: mech_df, spc_dct
        :rtype: dataframe, dict
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
        """ Identifies connected channels and assigns them to the same subpes

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param conn_chn_df: empty dataframe df[subpes][rxn]

        :returns: conn_chn_df dataframe[subpes][rxn]
        :rtype: dataframe[int][tuple]
        """

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
                fml_num = fml + key/100
                conn_chn_df['subpes'][rxns] = fml_num

        return conn_chn_df

    def group_species(self, reac_sp_df):
        """ Checks if the reactions in self.mech_df contain any species 
            of self.species_list and if so it marks the species.
            Assignment is hierarchical: the first species of species_list
            found determines the label of the reaction being checked
            WARNING FUNCTION STILL IN PROGRESS

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.species_list: list of subset of species considered
        :param reac_sp_df: empty dataframe index=rxns, column: 'species'

        :returns: reac_sp_df dataframe[species][rxn]
        :rtype: dataframe[str][tuple]
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
        """ Assigns a submechanism (fuel, fuel radical, fuel add ..) to the reaction
            according to the species taking part to it.
            WARNING FUNCTION STILL IN PROGRESS

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.species_list: list of subset of species considered
        :param self.species_subset_df: dataframe with species assigned to a 
                                        certain type (fuel, fuel radical..)
        :param submech_df: empty dataframe index=rxns, column: 'submech'

        :returns: submech_df dataframe with reactants mult dataframe[submech][rxn]
        :rtype: dataframe[str][tuple]
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
        """ determines the multiplicity of the reactants of the reactions
            considered

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param rxncl_mult_df: empty dataframe index=rxns, column: 'mult'

        :returns: reac_mult_df dataframe with reactants mult dataframe[mult][rxn]
        :rtype: dataframe[str][tuple]
        """
        # assign multiplicity values to each reactant
        for rxn in reac_mult_df.index:
            mult = 1
            rcts = self.mech_df['rct_names_lst'][rxn]
            mult = util.get_mult(rcts, self.spc_dct)
            reac_mult_df['mult'][rxn] = str(mult)

        return reac_mult_df

    def rxnclass_broad(self, rxncl_broad_df):
        """ assigns reaction classes based on stoichiometry and molecularity

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param rxncl_broad_df: empty dataframe index=rxns, column: 'rxn_class_broad'

        :returns: rxncl_broad_df dataframe with reaction classes dataframe[class][rxn]
        :rtype: dataframe[str][tuple]
        """
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
        """ assigns reaction classes using graph approach to all reactions 
            first subdivides the mech into subpeses; then classifies all rxn
            within the subpes, including wellskipping channels

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param rxncl_graph_df: empty dataframe index=rxns, column: 'rxn_class_graph'

        :returns: rxncl_graph_df dataframe with reaction classes dataframe[class][rxn]
        :rtype: dataframe[str][tuple]
        """
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

                    try:
                        rclass = classify_graph(
                            self.spc_dct, rct_names, prd_names)
                    except IndexError:
                        rclass = classify_graph_old(
                            self.spc_dct, rct_names, prd_names)

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
        """ assigns class headers based on the hierarchy provided in the input


        :param self.mech_df: dataframe with mech info
        :param hierarchy: list of strings with hierarchical order of sorting criteria
                            ['..','..','N']: N is the N of criteria to use for headers
                            all the other criteria are written as inline comments
        :param labels: labels corresponding to all possible sorting criteria
                        labels are then written as comments

        :returns: None. updates self.mech_df['cmts_top','cmts_inline']
                    'cmts_top': comments to write as header for a reaction
                    'cmts_inline': comments to write on the same line of the reaction
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
        """ provides sorted rxn indices and associated comments, sorted species dictionary

        :param self.mech_df: dataframe with mech info
        :param self.spc_dct: species dictionary for the reactions of the sorted mech

        :returns: sorted reaction names, comments, sorted species dictionary
        :rtype: list[tuple], dict{tuple(r1, r2, ): {cmts_top: str, cmts_inline: str}, ...}, dict
        """

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
    """ calls the graph classifier for a given reaction

    :param spc_dct: species dictionary
    :param rct_names: reactant names (r1, r2, )
    :param prd_names: product names (p1, p2, )

    :returns: reaction class (first of the possible identified classes)
    :rtype: str
    """

    if len(prd_names) >= 3:
        rclass = 'unclassified - lumped'
    else:

        # ID reaction
        rct_fmls = tuple(spc_dct[rct]['fml'] for rct in rct_names)
        prd_fmls = tuple(spc_dct[prd]['fml'] for prd in prd_names)

        rct_ichs = tuple(spc_dct[spc]['inchi'] for spc in rct_names)
        prd_ichs = tuple(spc_dct[spc]['inchi'] for spc in prd_names)

        if automol.formula.reac.is_valid_reaction(rct_fmls, prd_fmls):

            rxn_classes = automol.reac._find.find_from_inchis(
                rct_ichs, prd_ichs)

            if rxn_classes:
                # save only the first possible reaction type
                rclass = rxn_classes[0]
            else:
                rclass = 'unclassified'

        else:
            rclass = 'unclassified - Wrong Stoichiometry'
        # check stereo compatibility - I am not sure about this input
        # ret = automol.graph.trans.is_stereo_compatible(rclass, rct_graph, prd_graph)
    return rclass


def classify_ws(subpes_df, elem_reac_df, species_subpes, rxn):
    """ classifies well skipping channels of a given subpes
        WARNING: STILL UNDER CONSTRUCTION - SOME TEMPORARY FEATURES

    :param subpes_df: dataframe with subpes info
    :param elem_reac_df: dataframe with elementary reaction channels of the subpes
    :param species_subpes: list of subpes species
    :param rxn: string with rxn belonging to the subpes

    :returns: reaction class of the WS channel considered
    :rtype: str
    """
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
    """ assign comment strings based on reaction class

    :param name: reaction class (list, string, int, float)
    :param label: labels corresponding to reaction class (list)
    :param cltype: format type. options: class_head, class, subclass (string)

    :returns: comments string for writing in the mechanism
    :rtype: str
    """
    # assing top headers:
    tophead = '!!!!!!!!! class !!!!!!!!!\n'
    bottomhead = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
    if isinstance(name, str):
        name = [name]
    elif isinstance(name, int):
        name = [str(name)]
    elif isinstance(name, float):
        name = ['{:.2f}'.format(name)]
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

######################## Archived Functions ######################


def classify_graph_old(spc_dct, rct_names, prd_names):
    """ graph classification based on the first possible reaction class identified
    :param spc_dct: species dictionary
    :param rct_names: names of reactants (r1, r2, )
    :param prd_names: names of products (p1, p2, )

    :returns: reaction class
    :rtype: str
    """
    if len(prd_names) >= 3:
        rclass = 'unclassified - lumped'
    else:
        rct_graph = get_graph(spc_dct, rct_names)
        rct_gras = format_graph(rct_graph)

        prd_graph = get_graph(spc_dct, prd_names)
        prd_gras = format_graph(prd_graph)
        print('rct:', rct_gras, '\nprd:', prd_gras, '\n')
        # ID reaction
        rct_fmls = list(spc_dct[rct]['fml']
                        for rct in rct_names)
        prd_fmls = list(spc_dct[prd]['fml']
                        for prd in prd_names)

        if automol.formula.reac.is_valid_reaction(rct_fmls, prd_fmls):

            rxn_info = automol.reac._find.find(
                rct_gras, prd_gras)

            if rxn_info:
                # save only the first possible reaction type
                rclass = rxn_info[0].class_
            else:
                rclass = 'unclassified'

        else:
            rclass = 'unclassified - Wrong Stoichiometry'
        # check stereo compatibility - I am not sure about this input
        # ret = automol.graph.trans.is_stereo_compatible(rclass, rct_graph, prd_graph)
    return rclass


def get_graph(spc_dct, names):
    """ converst a list of species to a list of graphs
    :param spc_dct: species dictionary
    :param names: species names list [s1, s2, ]

    :returns: list of graphs
    :rtype: list[graph]
    """
    # Get the inchis and graphs
    spc_ichs = list(spc_dct[spc]['inchi']
                    for spc in names)
    graph_list = list(map(automol.inchi.graph, spc_ichs))

    return graph_list


def format_graph(gra):
    """ formats a list of graphs appropriately for rxn classification processing
        : no explicit hydrogens, no stereo, unique keys
    :param gra: list of graphs

    :returns: list of graphs suitable for classifier processing
    :rtype: list[graph]
    """
    '''
    gra: sequence of graphs
    returns: gra = sequence of graphs with appropriate format for the classifier
    - explicit hydrogens
    - no stereo parities
    - unique keys    
    '''
    gra_explicit = list(
        map(automol.graph._graph.explicit, gra))
    gra_nostereo = list(
        map(automol.graph._graph.without_stereo_parities, gra_explicit))
    gra_std_keys, _ = automol.graph._graph.standard_keys_for_sequence(
        gra_nostereo)

    return gra_std_keys
