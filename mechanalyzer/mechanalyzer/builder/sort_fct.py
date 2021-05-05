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
from mechanalyzer.builder import submech
from mechanalyzer.parser import pes
from mechanalyzer.parser._util import count_atoms
from mechanalyzer.parser._util import order_rct_bystoich
from mechanalyzer.parser._util import extract_spc
from mechanalyzer.parser._util import get_mult


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

        # Extract data from mech info
        [formula_dct_lst, formulas, rct_names_lst,
            prd_names_lst, thrdbdy_lst, rxn_name_lst, param_vals] = mech_info

        rxn_index = list(zip(rxn_name_lst, thrdbdy_lst))

        # Set dataframe: extract useful info
        pes_lst = count_atoms(formula_dct_lst)
        isthrdbdy = np.array(
            [(t[0] is not None and '(' not in t[0]) for t in thrdbdy_lst],
            dtype=int)
        molecularity = np.array(
            list(map(len, rct_names_lst)), dtype=int) + isthrdbdy
        n_of_prods = list(map(len, prd_names_lst))
        rct_names_lst_ordered = order_rct_bystoich(
            rct_names_lst, spc_dct=spc_dct)  # put heavier reactant first
        prd_names_lst_ordered = order_rct_bystoich(
            prd_names_lst, spc_dct=spc_dct)  # put heavier product first
        rct_1, rct_2 = extract_spc(rct_names_lst_ordered)
        data = np.array(
            [rct_names_lst, prd_names_lst,
             rct_names_lst_ordered, prd_names_lst_ordered,
             rct_1, rct_2, molecularity, n_of_prods, pes_lst, formulas,
             thrdbdy_lst, param_vals],
            dtype=object).T
        self.mech_df = pd.DataFrame(
            data, index=rxn_index,
            columns=['rct_names_lst', 'prd_names_lst', 'rct_names_lst_ord',
                     'prd_names_lst_ord', 'r1', 'r2', 'molecularity',
                     'N_of_prods', 'pes', 'formulas',
                     'thrdbdy', 'param_vals'])

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
        # set labels for all the possible criteria
        criteria_all = ['molecularity', 'N_of_prods', 'species',
                        'pes', 'subpes', 'submech',
                        'r1', 'mult',
                        'rxn_class_broad', 'rxn_class_graph',
                        'rxn_max_vals', 'rxn_max_ratio']
        labels_all = ['NR', 'N_of_prods', 'SPECIES',
                      'N_COH_PES', 'N_COH.subpes', 'SUBMECH',
                      'Heavier rct', 'Total multiplicity',
                      'rxn type', 'rxn class', 'max val', 'max ratio']
        # series: ascending/descending values
        asc_val = [True]*len(criteria_all)
        asc_val[-2:] = [False, False]
        asc_series = pd.Series(asc_val, index=criteria_all)

        # if species_list is not empty: pre-process the mechanism
        if len(species_list) > 0:
            if len(species_list) == 1 and 'submech' in hierarchy:
                # Select subset of species according to stoichiometries
                # specified in submech.py
                species_list, species_subset_df = submech.species_subset(
                    species_list[0], self.spc_dct)
                self.species_subset_df = species_subset_df
            elif len(species_list) > 1 and 'submech' in hierarchy:
                print('Error: pyr/ox submech extraction available ',
                      'for only 1 species')
                sys.exit()

            self.mech_df, self.spc_dct = self.filter_byspecies(species_list)
            self.species_list = species_list

        # 0. Look for keywords and sort accordingly
        # List of available sorting options with specific sorting functions
        sort_optns_dct = {
            'species': self.group_species,
            'subpes': self.conn_chn,
            'mult': self.reac_mult,
            'rxn_class_broad': self.rxnclass_broad,
            'rxn_class_graph': self.rxnclass_graph,
            'submech': self.group_submech,
            'rxn_max_vals': self.rxn_max_vals,
            'rxn_max_ratio': self.rxn_max_ratio
        }

        for optn, fun_name in sort_optns_dct.items():
            # Call sorting function
            if any(optn == inp_crt for inp_crt in hierarchy):
                # Generate corresponding dataframe
                df_optn = pd.DataFrame(
                    index=self.mech_df.index, columns=[optn])
                df_optn = fun_name(df_optn)
                # Concatenate the new portion of dataframe
                self.mech_df = pd.concat([self.mech_df, df_optn], axis=1)

        # 1. Sort
        try:
            self.mech_df = self.mech_df.sort_values(
                by=hierarchy[:-1], ascending=list(asc_series[hierarchy[:-1]].values))
        except KeyError as err:
            print(
                'WARNING: Reactions not sorted according ',
                'to all criteria: missing {}'.format(err))
            sys.exit()

        # 2. assign class headers
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
            print('Error in ISOLATE_SPECIES: ',
                  'not all species are in the species list ')
            sys.exit()

        spc_list = []
        # For all reactions in dataframe: check if species in selected list.
        # Otherwise remove the reaction
        for rxn in mech_df.index:
            rcts = list(mech_df['rct_names_lst'][rxn])
            prds = list(mech_df['prd_names_lst'][rxn])
            # Check if one species in the list is among reactants or products
            # of the reaction considered
            _rchk = any(rct == _spc for _spc in species_list for rct in rcts)
            _pchk = any(prd == _spc for _spc in species_list for prd in prds)
            if not _rchk and not _pchk:
                mech_df = mech_df.drop(index=[rxn])
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
            Generate pes dictionary for each reaction and save for later use

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param conn_chn_df: empty dataframe df[subpes][rxn]

        :returns: conn_chn_df dataframe[subpes][rxn]
        :rtype: dataframe[int][tuple]
        """
        pes_index = 0
        pes_dct_df = pd.DataFrame(
            index=self.mech_df.index,
            columns=['pes_dct'],
            dtype=object)
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
                # store in pes_dct_df
                fml_str = self.mech_df['formulas'][rxns].values[0]
                for rxn in rxns:
                    pes_dct_df['pes_dct'][rxn] = (fml_str, pes_index, key)

            pes_index += 1

        conn_chn_df = pd.concat([conn_chn_df, pes_dct_df], axis=1)

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

        # if species list is not found: do nothing, species entry remain empty
        if len(self.species_list) > 0:
            for rxn in reac_sp_df.index:
                rcts = list(self.mech_df['rct_names_lst'][rxn])
                prds = list(self.mech_df['prd_names_lst'][rxn])
                # check species hierarchically
                for sp_i in self.species_list:
                    if ((any(sp_i == rct for rct in rcts) or
                         any(sp_i == prd for prd in prds))
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

        :returns: submech_df dataframe with reactants
            mult dataframe[submech][rxn]
        :rtype: dataframe[str][tuple]
        """

        for rxn in submech_df.index:
            rcts = list(self.mech_df['rct_names_lst'][rxn])
            prds = list(self.mech_df['prd_names_lst'][rxn])
            # check species hierarchically (hierarchy fixed in species list)
            for sp_i in self.species_list:
                _rchk = any(sp_i == rct for rct in rcts)
                _pchk = any(sp_i == prd for prd in prds)
                _isfloat = isinstance(submech_df['submech'][rxn], float)
                if (_rchk or _pchk) and _isfloat:
                    submech_df['submech'][rxn] = self.species_subset_df[sp_i]
        return submech_df

    def reac_mult(self, reac_mult_df):
        """ determines the multiplicity of the reactants of the reactions
            considered

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param rxncl_mult_df: empty dataframe index=rxns, column: 'mult'

        :returns: reac_mult_df dataframe with reactants mult
            dataframe[mult][rxn]
        :rtype: dataframe[str][tuple]
        """
        # assign multiplicity values to each reactant
        for rxn in reac_mult_df.index:
            mult = 1
            rcts = self.mech_df['rct_names_lst'][rxn]
            mult = get_mult(rcts, self.spc_dct)
            reac_mult_df['mult'][rxn] = str(mult)

        return reac_mult_df

    def rxnclass_broad(self, rxncl_broad_df):
        """ assigns reaction classes based on stoichiometry and molecularity

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param rxncl_broad_df: empty dataframe
            index=rxns, column: 'rxn_class_broad'

        :returns: rxncl_broad_df dataframe with reaction classes
            dataframe[class][rxn]
        :rtype: dataframe[str][tuple]
        """
        for rxn in rxncl_broad_df.index:
            rcts = self.mech_df['rct_names_lst_ord'][rxn]
            prds = self.mech_df['prd_names_lst_ord'][rxn]
            _smol = (self.mech_df['molecularity'][rxn] == 1)
            _tbody = (self.mech_df['molecularity'][rxn] == 2
                      and self.mech_df['thrdbdy'][rxn][0] is not None)
            if _smol or _tbody:
                # unimolecular reaction classification
                rxn_class_broad = submech.classify_unimol(
                    rcts, prds, self.spc_dct)
            else:
                # bimolecular reaction classification
                rxn_class_broad = submech.classify_bimol(
                    rcts, prds, self.spc_dct)
            rxncl_broad_df['rxn_class_broad'][rxn] = rxn_class_broad

        return rxncl_broad_df

    def rxnclass_graph(self, rxncl_graph_df):
        """ assigns reaction classes using graph approach to all reactions
            first subdivides the mech into subpeses; then classifies all rxn
            within the subpes, including wellskipping channels

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param rxncl_graph_df: empty dataframe
            index=rxns, column: 'rxn_class_graph'

        :returns: rxncl_graph_df dataframe with reaction classes
            dataframe[class][rxn]
        :rtype: dataframe[str][tuple]
        """

        # 1. Group by subpes
        # Check SUBPES present in idxs, otherwise make corresponding dataframe
        if 'subpes' not in self.mech_df.columns:
            df_optn = pd.DataFrame(
                index=self.mech_df.index, columns=['subpes'])
            df_optn = self.conn_chn(df_optn)
            # concatenate the new portion of dataframe
            self.mech_df = pd.concat([self.mech_df, df_optn], axis=1)

        # 2. Graph classification or each subpes
        for _, subpes_df in self.mech_df.groupby('subpes'):
            # sort by molecularity: analyze first unimolecular isomerizations,
            # unimolecular decompositions, and then bimolecular reactions
            subpes_df = subpes_df.sort_values(
                by=['molecularity', 'N_of_prods'])
            # REFER TO REORDERED SPECIES NAMES,
            # OTHERWISE YOU MAY HAVE INCONSISTENT SPECIES NAMING
            # subpes species list
            species_subpes = list(set(
                list(subpes_df['rct_names_lst_ord'].values) +
                list(subpes_df['prd_names_lst_ord'].values)))

            # build dataframe: elementary reactivity matrix
            _nsubpes = len(species_subpes)
            elem_reac_df = pd.DataFrame(
                np.zeros((_nsubpes, _nsubpes), dtype='<U32'),
                index=species_subpes,
                columns=species_subpes)

            # graph classification
            for rxn in subpes_df.index:
                rct_names = subpes_df['rct_names_lst'][rxn]
                prd_names = subpes_df['prd_names_lst'][rxn]
                rct_names_ord = subpes_df['rct_names_lst_ord'][rxn]
                prd_names_ord = subpes_df['prd_names_lst_ord'][rxn]
                # Exclude rxns with more than 2 rcts or prds (not elementary!)
                if len(rct_names) < 3 and len(prd_names) < 3:

                    rclass = classify_graph(
                        self.spc_dct, rct_names, prd_names)

                    if rclass == 'unclassified':
                        if (subpes_df['molecularity'][rxn] == 1
                                and subpes_df['N_of_prods'][rxn] == 1):
                            rclass = 'isomerization'
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

    def rxn_max_vals(self, rxn_maxvals_df):
        """ Determines the maximum value of the rates in ktp dictionary.

        :param rxn_maxvals_df:
            empty dataframe index=rxns, column: 'rxn_max_vals'

        :returns: rxn_maxvals_df dataframe with overall max value of the rate
        :rtype: dataframe[float][tuple]
        """
        # extract maximum value for each ktp dictionary
        for rxn in rxn_maxvals_df.index:
            param_vals_dct = self.mech_df['param_vals'][rxn]
            max_val = get_max_aligned_values(param_vals_dct)
            rxn_maxvals_df['rxn_max_vals'][rxn] = max_val

        return rxn_maxvals_df

    def rxn_max_ratio(self, rxn_maxratio_df):
        """ Determines the maximum value of the ratios between
            different rates of ktp dct.

        :param rxn_maxratio_df:
            empty dataframe index=rxns, column: 'rxn_max_ratio'

        :returns: rxn_max_ratio dataframe with overall max value of the ratio
        :rtype: dataframe[float][tuple]
        """
        # extract maximum ratio for each set ktp dictionary
        for rxn in rxn_maxratio_df.index:
            param_vals_dct = self.mech_df['param_vals'][rxn]
            # get the ratio:
            param_ratio_dct = get_aligned_rxn_ratio_dct(param_vals_dct)
            max_val = get_max_aligned_values(param_ratio_dct)
            rxn_maxratio_df['rxn_max_ratio'][rxn] = max_val

        return rxn_maxratio_df

    # ASSIGN HEADERS #

    def class_headers(self, hierarchy, labels):
        """ assigns class headers based on the hierarchy provided in the input


        :param self.mech_df: dataframe with mech info
        :param hierarchy: list of strings with hierarchical order
            of sorting criteria ['..','..','N']: N is of criteria for headers
            all the other criteria are written as inline comments
        :param labels: labels corresponding to all possible sorting criteria
                        labels are then written as comments

        :returns: None. updates self.mech_df['cmts_top','cmts_inline']
                    'cmts_top': comments to write as header for a reaction
                    'cmts_inline': comments to write on same line of reaction
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
                '*ERROR: Last line of sorting options is the N ',
                'of criteria to use for class headers')
            sys.exit()

        # Write topheader comments
        if n_headers > 0:
            for name, rdf in self.mech_df.groupby(hierarchy[:n_headers]):
                # Write rxn class as top header comments
                rxnclass = cmts_string(
                    name, labels[hierarchy[:n_headers]], 'class_head')
                idx0 = rdf.index[0]
                df_cmts_top['cmts_top'][idx0] = rxnclass

                # Write inline comments if necessary
                if n_headers < len(hierarchy)-1:
                    for name2, rdf2 in rdf.groupby(hierarchy[n_headers:-1]):
                        rxnclass = cmts_string(
                            name2, labels[hierarchy[n_headers:-1]], 'subclass')
                        df_cmts_inline['cmts_inline'][rdf2.index] = rxnclass
        else:
            # Write only inline comments
            for name, rdf in self.mech_df.groupby(hierarchy[n_headers:-1]):
                rxnclass = cmts_string(
                    name, labels[hierarchy[n_headers:-1]], 'class')
                df_cmts_inline['cmts_inline'][rdf.index] = rxnclass

        # Concatenate DFs
        self.mech_df = pd.concat(
            [self.mech_df, df_cmts_top, df_cmts_inline], axis=1)

    # OUTPUT DATAFRAMES and DICTIONARIES #
    def return_mech_df(self):
        """ Provides sorted rxn indices and associated comments, as well as
            sorted species dictionary.

        :param self.mech_df: dataframe with mech info
        :param self.spc_dct: species dictionary for reactions of sorted mech

        :returns: sorted reaction names, comments, sorted species dictionary
        :rtype: list[tuple],
                dict{tuple(r1, r2,): {cmts_top: str, cmts_inline: str}, ...},
                dict
        """

        rct_names = self.mech_df['rct_names_lst'].values
        prd_names = self.mech_df['prd_names_lst'].values
        thrdbdy = self.mech_df['thrdbdy'].values
        new_idx = list(zip(rct_names, prd_names, thrdbdy))
        # store comments in dct
        cmts_df = pd.DataFrame(
            self.mech_df[['cmts_top', 'cmts_inline']].values,
            index=new_idx,
            columns=['cmts_top', 'cmts_inline'])
        cmts = cmts_df.to_dict('index')

        return new_idx, cmts, self.spc_dct

    def return_pes_dct(self):
        """ returns a PES dictionary

            chn idx is set to to OVERALL PES, not the SUB_PES

        :returns: pes_dct: {('fml', n_pes, n_subpes): (chn_idx, (rcts, prds)), ...}
        :rtype: dct{tuple: tuple}
        """

        # Check SUBPES present in idxs, otherwise make corresponding dataframe
        if 'subpes' not in self.mech_df.columns:
            df_optn = pd.DataFrame(
                index=self.mech_df.index, columns=['subpes'])
            df_optn = self.conn_chn(df_optn)
            # concatenate the new portion of dataframe
            self.mech_df = pd.concat([self.mech_df, df_optn], axis=1)

        # get the pes dictionary
        pes_dct = {}
        prev_idx = -1
        for _, pes_dct_df in self.mech_df.groupby('subpes'):

            # Get the ('fml', n_pes, n_subpes) for dict key
            pes_dct_key = pes_dct_df['pes_dct'].values[0]

            # Get the names of the reaction reagents
            rct_names = pes_dct_df['rct_names_lst'].values
            prd_names = pes_dct_df['prd_names_lst'].values
            # thrdbdy = pes_dct_df['thrdbdy'].values
            # new_idx = list(zip(rct_names, prd_names, thrdbdy))
            rxn_names = tuple(zip(rct_names, prd_names))

            # Count number of channels on the previous SUB-PES (if there is one)
            # Used to increment the channel idxs on the current SUB-PES
            # Checks the idx to see if it has changed from previous
            _, pidx, _ = pes_dct_key
            if pidx == prev_idx:
                nchnls += len(rxn_names)
            else:
                nchnls = 0
            prev_idx = pidx

            # Build the full channel list for the SUB-PES
            chnl_lst = tuple((idx+nchnls, rxn)
                             for idx, rxn in enumerate(rxn_names))

            pes_dct[pes_dct_key] = chnl_lst

        return pes_dct


# FUNCTIONS FOR RXN GRAPH CLASSIFICATION #


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
            try:
                rxn_objs = automol.reac.rxn_objs_from_inchi(
                    rct_ichs, prd_ichs)
                rxn_classes = tuple(obj[0].class_ for obj in rxn_objs)
            except AssertionError:
                rxn_classes = ('AssertionError', )
            except TypeError:
                print('geoms of rxn classifier fail for rxn: {} = {}'.format(
                    rct_ichs, prd_ichs))
                rxn_classes = ('TypeError', )

            if rxn_classes:
                # save only the first possible reaction type
                # rclass = rxn_classes[0]
                rclass = '/'.join(set(rxn_classes))
            else:
                rclass = 'unclassified'

        else:
            rclass = 'unclassified - Wrong Stoichiometry'

    return rclass


def classify_ws(subpes_df, elem_reac_df, species_subpes, rxn):
    """ classifies well skipping channels of a given subpes
        WARNING: STILL UNDER CONSTRUCTION - SOME TEMPORARY FEATURES

    :param subpes_df: dataframe with subpes info
    :param elem_reac_df: dataframe with elementary reaction channels of subpes
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
        # TEMPORARY: SHOULD RECONSTRUCT FULL PATH FROM REACTANTS TO PRODUCTS
        rxn_type_1 = rxn_types_1[0]
        # TEMPORARY: SHOULD RECONSTRUCT FULL PATH FROM REACTANTS TO PRODUCTS
        rxn_type_2 = rxn_types_2[0]

        # WRITE THE REACTION TYPE STRING
        rxn_type_ws = rxn_type_1 + '-' + rxn_type_2 + ' (WS)'
        return rxn_type_ws

    except IndexError:
        return None


# FUNCTIONS FOR COMMENTS - CALLED BY THE SORTER #

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


# FUNCTIONS WORKING WITH KTP DICTIONARIES #

def get_aligned_rxn_ratio_dct(aligned_rxn_dct_entry):
    """ converts the entry of the aligned_rxn_ktp_dictionary to the ratios
        between the reference rate and the given rate

    :param aligned_rxn_dct_entry: entry of aligned_rxn_ktp/ratio_dct
    :type aligned_rxn_dct_entry:
        list[dct{pressure: np.array(temps), np.array(values)}]
    :return aligned_ratio_dct_entry: aligned dictionary entry
    :rtype: list(dct)
    """

    ref_ktp_dct = aligned_rxn_dct_entry[0]
    ratio_dct_entry = []
    for mech_idx, ktp_dct in enumerate(aligned_rxn_dct_entry):
        # If (1) on the first ktp_dct, (2) the ref_ktp_dct is None,
        # or (3) the current_ktp_dct is None, set the ratio_dct to None
        if mech_idx == 0 or ref_ktp_dct is None or ktp_dct is None:
            ratio_dct = None
        # Otherwise, calculate the ratio_dct
        else:
            ratio_dct = {}
            for pressure, (temps, kts) in ktp_dct.items():
                # If pressure defined in ref ktp_dct: calculate and store ratio
                if pressure in ref_ktp_dct.keys():
                    _, ref_kts = ref_ktp_dct[pressure]
                    ratios = kts / ref_kts
                    ratio_dct[pressure] = (temps, ratios)
            if ratio_dct == {}:  # account for when no pressures contain ratios
                ratio_dct = None

        # Append the current ratio_dct
        ratio_dct_entry.append(ratio_dct)

    return ratio_dct_entry


def get_max_aligned_values(aligned_rxn_dct_entry):
    """ Gets the maximum values for each reaction from an entry (value) of
        either an aligned_rxn_ktp_dct or an aligned_rxn_ratio_dct

    :param aligned_rxn_dct_entry: entry of aligned_rxn_ktp/ratio_dct
    :type aligned_rxn_dct_entry:
        list[dct{pressure: np.array(temps), np.array(values)}]
    :return max_val: max value
    :rtype: float
    """

    max_val = 0
    for single_dct in aligned_rxn_dct_entry:
        if single_dct is not None:
            for _, (_, values) in single_dct.items():
                if max(values) > max_val:
                    max_val = max(values)

    return max_val
