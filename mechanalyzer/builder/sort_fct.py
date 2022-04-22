"""
Sorter module - sorting of the mechanism according to various options
- reaction classes
- pes/subpes
- multiplicity
- species subsets
- submechanism
"""
import time
import sys
import copy
import pandas as pd
import numpy
import automol
from mechanalyzer.builder import submech
from mechanalyzer.parser import pes
from mechanalyzer.parser._util import count_atoms
from mechanalyzer.parser._util import order_rct_bystoich
from mechanalyzer.parser._util import extract_spc
from mechanalyzer.parser._util import get_mult
from mechanalyzer.parser.spc import name_inchi_dct
from mechanalyzer.parser._util import get_fml


class SortMech:
    """ class of methods to organize the mechanism according to given criteria
    """

    def __init__(self, rxn_param_dct, spc_dct):
        """ Initializes the mechanism dataframe and the species dictionary

        :param rxn_param_dct: rxn param dct for info extraction
        :param spc_dct: species dictionary

        :returns: None, updates self.
                    self.mech_df: dataframe with mech info
                    self.spc_dct: species dictionary
        """

        # Extract data from mech info
        [formula_dct_lst, formulas, rct_names_lst,
            prd_names_lst, thrdbdy_lst, rxn_name_lst, param_vals] = mech_info(rxn_param_dct, spc_dct)

        rxn_index = list(zip(rxn_name_lst, thrdbdy_lst))

        # Set dataframe: extract useful info
        pes_lst = count_atoms(formula_dct_lst)
        isthrdbdy = numpy.array(
            [(t[0] is not None and '(' not in t[0]) for t in thrdbdy_lst],
            dtype=int)
        molecularity = numpy.array(
            list(map(len, rct_names_lst)), dtype=int) + isthrdbdy
        n_of_prods = list(map(len, prd_names_lst))
        rct_names_lst_ordered = order_rct_bystoich(
            rct_names_lst, spc_dct=spc_dct)  # put heavier reactant first
        prd_names_lst_ordered = order_rct_bystoich(
            prd_names_lst, spc_dct=spc_dct)  # put heavier product first
        rct_1, rct_2 = extract_spc(rct_names_lst_ordered)
        data = numpy.array(
            [rct_names_lst, prd_names_lst,
             rct_names_lst_ordered, prd_names_lst_ordered,
             rct_1, rct_2, molecularity, n_of_prods, pes_lst, formulas,
             isthrdbdy, thrdbdy_lst, param_vals, rxn_name_lst],
            dtype=object).T
        self.mech_df = pd.DataFrame(
            data, index=rxn_index,
            columns=['rct_names_lst', 'prd_names_lst', 'rct_names_lst_ord',
                     'prd_names_lst_ord', 'r1', 'r2', 'molecularity',
                     'N_of_prods', 'pes_fml', 'formulas', 'isthrdbdy',
                     'thrdbdy', 'param_vals', 'rxn_names'])

        # reindex pes
        df_pes = pd.DataFrame(
            index=self.mech_df.index, columns=['pes'])
        pes_index = 1
        for _, pes_fml in self.mech_df.groupby('pes_fml'):
            rxns = pes_fml.index
            df_pes['pes'][rxns] = pes_index
            pes_index += 1
        # Concatenate the new portion of dataframe
        self.mech_df = pd.concat([self.mech_df, df_pes], axis=1)

        self.spc_dct = spc_dct  # set for later use
        # empty list for initialization (otherwise pylint warning)
        self.species_subset_df = ()
        self.species_list = ()

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
        # maybe turn into dataframe with criteria, labels, fct names, and ascending opts
        criteria_all = ['molecularity', 'N_of_prods', 'species',
                        'pes', 'subpes', 'chnl',
                        'r1', 'mult',
                        'rxn_class_broad', 'rxn_class_graph',
                        'submech', 'submech_prompt', 'submech_ext',
                        'rxn_max_vals', 'rxn_max_ratio',
                        ]
        labels_all = ['NR', 'N_of_prods', 'SPECIES',
                      'pes', 'subpes', 'channel',
                      'Heavier_rct', 'Totalmultiplicity',
                      'rxntype', 'rxnclass',
                      'SUBMECH', 'submech_prompt', 'submech_ext',
                      'maxval', 'maxratio',
                      ]

        # check on pes/subpes criteria:
        # if subpes alone, also add pes
        if (('subpes' in hierarchy and 'pes' not in hierarchy) or
                ('chnl' in hierarchy and ['pes', 'subpes'] not in hierarchy)):
            try:
                idx_insert_pes = hierarchy.index('subpes')
            except ValueError:
                idx_insert_pes = hierarchy.index('chnl')
            hierarchy.insert(idx_insert_pes, 'pes')
            # chenge N of headers if necessary
            hierarchy[-1] += 1*(hierarchy[-1] > idx_insert_pes)
            print('pes criterion added: necessary for subpes/chnls')

        # if species_list is not empty: pre-process the mechanism
        if len(species_list) > 0:
            if len(species_list) >= 1 and 'submech_prompt' in hierarchy:
                filtertype = 'submech_prompt'
                # add pes, subpes, chnl:
                hierarchy_new = ['pes', 'subpes', 'chnl']
                [hierarchy_new.append(h)
                 for h in hierarchy if h not in hierarchy_new]
                hierarchy = hierarchy_new
            elif len(species_list) == 1 and 'submech' in hierarchy:
                # Select subset of species according to stoichiometries
                # specified in submech.py
                species_list, species_subset_df = submech.species_subset(
                    species_list[0], self.spc_dct)
                self.species_subset_df = species_subset_df
                filtertype = 'submech'

            elif len(species_list) == 1 and 'submech_ext' in hierarchy:
                species_list, species_subset_df = submech.species_subset_ext(
                    species_list[0], self.spc_dct)
                self.species_subset_df = species_subset_df
                filtertype = 'submech_ext'

            elif len(species_list) > 1 and 'submech' in hierarchy:
                print('Error: pyr/ox submech extraction available ',
                      'for only 1 species')
                sys.exit()

            else:
                filtertype = 'submech'

            self.mech_df, self.spc_dct = self.filter_byspecies(
                species_list, filtertype)
            self.species_list = species_list

        # 0. Look for keywords and sort accordingly
        # List of available sorting options with specific sorting functions
        sort_optns_dct = {
            'species': self.group_species,
            'subpes': self.chnl,
            'chnl': self.chnl,
            'mult': self.reac_mult,
            'rxn_class_broad': self.rxnclass_broad,
            'rxn_class_graph': self.rxnclass_graph,
            'submech': self.group_submech,
            'submech_prompt': self.group_prompt,
            'submech_ext': self.group_submech,
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
        # series: ascending/descending values
        asc_val = [True]*len(criteria_all + ['rxn_names'])
        # rxn vals, ratio and names are descending
        asc_val[-3:] = [False, False, False]
        asc_series = pd.Series(asc_val, index=criteria_all + ['rxn_names'])

        try:
            # last "standard" criterion: rxn name
            asclst = list(asc_series[hierarchy[:-1] + ['rxn_names']].values)
            self.mech_df = self.mech_df.sort_values(
                by=hierarchy[:-1] + ['rxn_names'],
                ascending=asclst)
        except KeyError as err:
            print(
                'Error: Reactions not sorted according ',
                f'to all criteria: missing {err}, exiting')
            sys.exit()

        # 2. assign class headers
        labels = pd.Series(labels_all, index=criteria_all)
        self.class_headers(hierarchy, labels)

    def filter_byspecies(self, species_list, filtertype):
        """ Find all reactions involving species of the species_list given as input
            Provides a new mechanism of the reactions of the selected subset
            and a reduced list of species

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param species_list: list of species subsets
        :param filtertype: can be 'submech', 'submech_ext', 'submech_prompt'

        :returns: mech_df, spc_dct
        :rtype: dataframe, dict
        """
        # add classification by subpes first
        if filtertype == 'submech_prompt':
            self.mech_df = pd.concat(
                [self.mech_df, self.chnl('')], axis=1)  # add subpes
            self.mech_df['prompt'] = ''
        # deepcopy
        mech_df = copy.deepcopy(self.mech_df)
        # check that all species selected are in the species dictionary
        if any(i not in self.spc_dct.keys() for i in species_list):
            print('Error in ISOLATE_SPECIES: ',
                  'not all species are in the species list. Exiting')
            sys.exit()

        # For all reactions in dataframe: check if species in selected list.
        # Otherwise remove the reaction
        spc_list = []

        for rxn in mech_df.index:
            rcts = list(mech_df['rct_names_lst'][rxn])
            prds = list(mech_df['prd_names_lst'][rxn])
            # Check if one species in the list is among reactants or products
            # of the reaction considered
            if filtertype in ['submech', 'submech_prompt']:
                _rchk = any(
                    rct == _spc for _spc in species_list for rct in rcts)
                _pchk = any(
                    prd == _spc for _spc in species_list for prd in prds)
                chk = int(_rchk or _pchk)

            elif filtertype == 'submech_ext':
                # filter out bimol/bimol reactions if some bimol rcts/prds are not in the species list
                # equivalent to saying: at least all reactants (also single react works) or all products
                # must be in the species list
                chk = 0
                chk += int(all(any(rct == _spc for _spc in species_list)
                               for rct in rcts))
                chk += int(all(any(prd == _spc for _spc in species_list)
                               for prd in prds))*int(len(prds) <= 2)

            if filtertype == 'submech_prompt' and chk >= 1:
                # subpes list to keep
                # rxn is bimol on both sides
                if len(rcts) == 2 and len(prds) == 2:
                    mech_df['prompt'][rxn] = 'RAD_GEN'
                elif len(prds) > 2 and any(pr in species_list for pr in prds):
                    mech_df['prompt'][rxn] = 'PROMPT_LUMPED'
                # rxn is unimol deco/formation of the radical
                elif ((len(rcts) == 1 and rcts[0] in species_list and len(prds) == 2) or
                      (len(prds) == 1 and prds[0] in species_list and len(rcts) == 2)):
                    mech_df['prompt'][rxn] = 'RAD_DECO'

                continue  # next for loop

            if chk >= 1:  # keep reaction and add to species list
                spc_list.extend(rcts)
                spc_list.extend(prds)
            elif chk == 0:  # reaction filtered out
                mech_df = mech_df.drop(index=[rxn])

        if filtertype == 'submech_prompt':
            # submech_prompt: if RAD_GEN/RAD_DECO in list, keep the rxns
            for _, subpes_df in mech_df.groupby(['pes', 'subpes']):
                if any(pr in ['RAD_GEN', 'RAD_DECO', 'PROMPT_LUMPED'] for pr in subpes_df['prompt'].values):
                    [spc_list.extend(list(rcts_tup))
                     for rcts_tup in mech_df['rct_names_lst'].values]
                    [spc_list.extend(list(prds_tup))
                     for prds_tup in mech_df['prd_names_lst'].values]
                else:
                    # drop the corrisponding pes/subpes
                    mech_df = mech_df.drop(index=subpes_df.index)

        # filter spc_list: unique elements
        spc_list = list(set(spc_list))
        # new spc_dct
        spc_dct_val = list(map(self.spc_dct.get, spc_list))
        spc_dct = dict(zip(spc_list, spc_dct_val))
        # print('line 298 in sort_fct', mech_df)
        return mech_df, spc_dct

    def conn_chn(self, conn_chn_df):
        """ Identifies connected channels and assigns them to the same subpes
            Generate pes dictionary for each reaction and save for later use

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param conn_chn_df: empty dataframe df[subpes][rxn]

        :returns: conn_chn_df dataframe[subpes][rxn]
        :rtype: dataframe[int][tuple]
        """

        pes_dct_df = pd.DataFrame(
            index=self.mech_df.index,
            columns=['chnl', 'pes_chnl_tuple'],
            # columns=['pes_dct', 'pes_chnl_tuple', 'pes_chnl'],
            dtype=object)

        for _, peslist in self.mech_df.groupby('pes'):
            idx_start = 0
            # Set the names lists for the rxns and species needed below
            peslist = peslist.sort_values(by=['rxn_names'])
            pes_rct_names_lst = peslist['rct_names_lst'].values
            pes_prd_names_lst = peslist['prd_names_lst'].values
            pes_rxn_name_lst = peslist.index
            connchnls = pes.find_conn_chnls(
                pes_rct_names_lst, pes_prd_names_lst, pes_rxn_name_lst)

            # Write subpes in conn_chn_df
            for key, value in connchnls.items():
                rxns = peslist.iloc[value].sort_values(
                    by=['rxn_names'], ascending=False).index
                conn_chn_df['subpes'][rxns] = key+1
                # reorder by rxn name before assigning the channel index

                for chnl_idx, rxn in enumerate(rxns):
                    rxn_name = tuple(
                        (peslist['rct_names_lst'][rxn],
                         peslist['prd_names_lst'][rxn]))
                    pes_dct_df['chnl'][rxn] = chnl_idx+idx_start+1
                    # pes_dct_df['pes_dct'][rxn] = (fml_str, pes_index, key)
                    pes_dct_df['pes_chnl_tuple'][rxn] = (
                        chnl_idx+idx_start, rxn_name)
                    # pes_dct_df['pes_chnl'][rxn] = ','.join(
                    #    [str(pes_index+1), str(key+1),
                    #     str(chnl_idx+idx_start+1)])

                idx_start += len(rxns)

        conn_chn_df = pd.concat([conn_chn_df, pes_dct_df], axis=1)

        return conn_chn_df

    def chnl(self, _):
        """ Calls subpes if not done
        """
        if 'subpes' not in self.mech_df.columns:
            df_optn = pd.DataFrame(
                index=self.mech_df.index, columns=['subpes'])
            ret_df = self.conn_chn(df_optn)
        else:  # no need to do anything; put empty df to avoid duplicate idxing
            ret_df = pd.DataFrame(
                index=self.mech_df.index, columns=[''])
        return ret_df

    def group_species(self, reac_sp_df):
        """ Checks if the reactions in self.mech_df contain any species
            of self.species_list and if so it marks the species.
            Assignment is hierarchical: the first species of species_list
            found determines the label of the reaction being checked

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

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.species_list: list of subset of species considered
        :param self.species_subset_df: dataframe with species assigned to a
                                        certain type (fuel, fuel radical..)
        :param submech_df: empty dataframe index=rxns, column: 'submech'
                                                        or 'submech_ext'
        :returns: submech_df dataframe with reactants
            mult dataframe[submech][rxn]
        :rtype: dataframe[str][tuple]
        """
        # lbl as submech_df columns -> so it works with both submech and submech_ext
        lbl_col = submech_df.columns[0]

        for rxn in submech_df.index:
            rcts = list(self.mech_df['rct_names_lst'][rxn])
            prds = list(self.mech_df['prd_names_lst'][rxn])
            # check species hierarchically (hierarchy fixed in species list)
            for sp_i in self.species_list:
                _rchk = any(sp_i == rct for rct in rcts)
                _pchk = any(sp_i == prd for prd in prds)
                _isfloat = isinstance(submech_df[lbl_col][rxn], float)
                if (_rchk or _pchk) and _isfloat:
                    submech_df[lbl_col][rxn] = self.species_subset_df[sp_i]

        return submech_df

    def group_prompt(self, submech_df):
        """ group reactions according to prompt dissociation channels
            builds group dictionary list:
            :param grps: [{'grp': 0, 'idxs': [], 'peds': [], 'hot': []}, ...]
            :type grps: list(dct)
        """
        # assign group values:
        # 1. groups by PES
        # 2. if PES is rad_deco: renames as rad_deco_radname
        # 3. if PES produces hot radical: adds it to a group ['grp N']
        #       then specifies the name in ['groupname'] nb add to same list if same subpes
        submech_df['submech_prompt'] = 'unclassified'
        grps = []
        species_deco_dct = dict.fromkeys(self.species_list)
        # filter by RAD_DECO: save the pes/subpes value
        hot_df = self.mech_df[self.mech_df['prompt'] == 'RAD_DECO']
        for rxn in hot_df.index:
            for sp in self.species_list:
                if sp in (hot_df['rct_names_lst'][rxn][0], hot_df['prd_names_lst'][rxn][0]):
                    species_deco_dct[sp] = '{}:{}'.format(
                        hot_df['pes'][rxn], hot_df['subpes'][rxn])
                    submech_df['submech_prompt'][rxn] = 'RAD_DECO_{}'.format(
                        sp)

        prompt_df = self.mech_df[self.mech_df['prompt'] == 'PROMPT_LUMPED']
        for rxn in prompt_df.index:
            for sp in self.species_list:
                if sp in prompt_df['prd_names_lst'][rxn]:
                    submech_df['submech_prompt'][rxn] = 'PROMPT_LUMPED_{}'.format(
                        sp)

        grp_dct_template = {'grp': 0, 'idxs': [], 'peds': [], 'hot': []}
        grpN = 0
        ped_df = self.mech_df[self.mech_df['prompt'] == 'RAD_GEN']
        for pes, pesdf in ped_df.groupby('pes'):
            grpN += 1
            hotspc_dct = {}  # 'pes:subpes':[sp1, sp2]
            grp_dct = copy.deepcopy(grp_dct_template)
            grp_dct['grp'] = grpN

            for subpes, subpesdf in pesdf.groupby('subpes'):
                peds = []
                grp_dct['hot'].append([])
                grp_dct['idxs'].append('{}:{}'.format(pes, subpes))

                for rxn in subpesdf.index:
                    # sp in products
                    if any(sp in subpesdf['prd_names_lst'][rxn] for sp in self.species_list):
                        sp = self.species_list[[
                            i in subpesdf['prd_names_lst'][rxn] for i in self.species_list].index(True)]
                        rxn_ped = rxn[0]
                    # species in reactants
                    elif any(sp in subpesdf['rct_names_lst'][rxn] for sp in self.species_list):
                        sp = self.species_list[[
                            i in subpesdf['rct_names_lst'][rxn] for i in self.species_list].index(True)]
                        r1, r2 = subpesdf['rct_names_lst'][rxn]
                        p1, p2 = subpesdf['prd_names_lst'][rxn]
                        # switch rcts and prds
                        rxn_ped = '{}+{}={}+{}'.format(p1, p2, r1, r2)
                    submech_df['submech_prompt'][rxn] = 'RAD_GEN_{}'.format(sp)
                    peds.append(rxn_ped)
                    # deal with hotspecies
                    hotpes = species_deco_dct[sp]
                    if hotpes in hotspc_dct.keys():
                        hotspc_dct[hotpes].append(sp)
                    else:
                        hotspc_dct[hotpes] = [sp]
                grp_dct['peds'].append(peds)
            # assign hot species
            for hotpes in hotspc_dct.keys():
                grp_dct['peds'].append([])
                grp_dct['idxs'].append(hotpes)
                grp_dct['hot'].append(list(set(hotspc_dct[hotpes])))

            grps.append(grp_dct)

        self.grps = grps

        return submech_df

    def filter_groups_prompt(self, therm_dct, threshold=30.):
        """ filter prompt groups according to
            A+B->B+C->B+(D+E) must be < 30 kcal/mol endothermic
            otherwise remove from dct
        """
        # 1 find min endothermicity of deco reactions of hotspecies
        dh_min_hot = dict.fromkeys(self.species_list)
        for hot_sp in self.species_list:
            hot_sp_df = self.mech_df[self.mech_df['submech_prompt']
                                     == 'RAD_DECO_{}'.format(hot_sp)]
            dh = []
            for rxn in hot_sp_df.index:
                rcts = hot_sp_df['rct_names_lst'][rxn]
                prds = hot_sp_df['prd_names_lst'][rxn]
                dh_rcts = [therm_dct[rct][1][0]/1000 for rct in rcts]
                dh_prds = [therm_dct[prd][1][0]/1000 for prd in prds]
                dh_rxn = sum(dh_prds) - sum(dh_rcts)
                if hot_sp in prds:  # revert sign
                    dh_rxn = -dh_rxn
                dh.append(dh_rxn)
                # print(rcts, prds, dh_rcts, dh_prds, dh_rxn)
            dh_min_hot[hot_sp] = min(dh)

        # loop over groups and delete too endothermic reactions
        filtered_grps = []
        grp_idx = 0
        for grp in self.grps:
            check = 0
            # lists, potentially more than 1
            for n, ped in enumerate(grp['peds']):
                for ped_i in ped:
                    # extract rxn exo/endo thermicity
                    rcts = ped_i.split('=')[0].split('+')
                    prds = ped_i.split('=')[1].split('+')
                    dh_rcts = [therm_dct[rct][1][0]/1000 for rct in rcts]
                    dh_prds = [therm_dct[prd][1][0]/1000 for prd in prds]
                    dh = sum(dh_prds) - sum(dh_rcts)
                    # check with threshold
                    hot_sp = list(set(self.species_list).intersection(prds))[0]
                    dh_tot = dh + dh_min_hot[hot_sp]
                    # print(ped_i, dh, dh_min_hot[hot_sp], dh_tot)
                    if dh_tot > threshold:
                        ped.remove(ped_i)
                        print('*Warning: rxn {} removed from pes_groups \
                              min overall DH of {:.2f} kcal/mol'.format(ped_i, dh_tot))
                    else:
                        check += 1

                if not ped and not grp['hot'][n]:  # remove all indices from dct
                    del grp['idxs'][n]
                    del grp['peds'][n]
                    del grp['hot'][n]

            if check > 0:
                grp_idx += 1
                grp['grp'] = grp_idx
                filtered_grps.append(grp)

        self.grps = filtered_grps

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
                      and self.mech_df['isthrdbdy'][rxn] == 1)

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

        # 1. Group by subpes if absent
        self.mech_df = pd.concat([self.mech_df, self.chnl('')], axis=1)

        # 2. Graph classification or each subpes
        for _, subpes_df in self.mech_df.groupby(['pes', 'subpes']):
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
                numpy.zeros((_nsubpes, _nsubpes), dtype='<U32'),
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
        ept_df = numpy.zeros((len(self.mech_df.index), 1), dtype=str)
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

            chn idx is set to the OVERALL PES, not the SUB_PES

        :returns: pes_dct: {(fml, n_pes, n_subpes): ((chn_idx, (rcts, prds)),)}
        :rtype: dct{tuple: tuple}
        """

        # Check by subpes if absent
        self.mech_df = pd.concat([self.mech_df, self.chnl('')], axis=1)

        # get the pes dictionary
        pes_dct = {}

        for pes_idx, pes_all_df in self.mech_df.groupby('pes'):
            fml_str = pes_all_df['formulas'].values[0]

            for subpes_idx, pes_dct_df in pes_all_df.groupby('subpes'):
                # Get the ('fml', n_pes, n_subpes) for dict key
                pes_dct_df = pes_dct_df.sort_values(by='chnl')
                pes_dct_key = (fml_str, pes_idx-1, subpes_idx-1)
                pes_chnl_tuples = pes_dct_df['pes_chnl_tuple'].values
                chnl_lst = tuple(pes_chnl_tuples)
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
                print('geoms of rxn classifier fail for rxn: '
                      f'{rct_ichs} = {prd_ichs}')
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
    head = '!!!!!!!!!!!!!!!!!!!!!!!!!\n'

    try:
        labeldct = dict(zip(name, label.values))
    except TypeError:
        labeldct = dict(zip([name], [label.values]))

    def formatname(name):
        if isinstance(name, int):
            name = str(name)
        elif isinstance(name, numpy.int64):
            name = str(name)
        elif isinstance(name, float):
            if labeldct[name] not in ['maxval', 'maxratio']:
                name = f'{name:.2f}'
            else:
                name = f'{name:.2e}'

        return name

    if isinstance(name, tuple):
        namenew = []
        for _nm in name:
            namenew.append(formatname(_nm))
        name = namenew
    else:
        name = [formatname(name)]

    if cltype == 'class_head':
        cmtlabel = '!       '+'.'.join(label)
        cmtlabel += '!       '+'.'.join(name)+'\n'
        rxnclass = head + cmtlabel + head
    else:
        cmtlabel = '.'.join(label) + '  '
        rxnclass = '! ' + cmtlabel + '.'.join(name)

    return rxnclass


# FUNCTIONS WORKING WITH KTP DICTIONARIES #

def get_aligned_rxn_ratio_dct(aligned_rxn_dct_entry):
    """ converts the entry of the aligned_rxn_ktp_dictionary to the ratios
        between the reference rate and the given rate

    :param aligned_rxn_dct_entry: entry of aligned_rxn_ktp/ratio_dct
    :type aligned_rxn_dct_entry:
        list[dct{pressure: numpy.array(temps), numpy.array(values)}]
    :return aligned_ratio_dct_entry: aligned dictionary entry
    :rtype: list(dct)
    """

    ref_ktp_dct = aligned_rxn_dct_entry[0]
    ratio_dct_entry = []
    for mech_idx, ktp_dct in enumerate(aligned_rxn_dct_entry):
        # If (1) you are on the first ktp_dct, (2) the ref_ktp_dct is None,
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
        list[dct{pressure: numpy.array(temps), numpy.array(values)}]
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


# EXTRACT MECH INFO - PREVIOUSLY IN MECHANALYZER PARSER

def mech_info(rxn_param_dct, spc_dct):
    """ Build mech_info object for mech sorting

        :param spc_dct: species dictionary
        :type spc_dct: dict[?:?]
        :param rxn_dct: parameter dictionary
        :type rxn_dct: dict[?:?]
        :return mech_info: objects with mech info
        :rtype: list
    """

    def _check_names(rct_names, prd_names, all_spc_names):
        """ Assess if the reactant and product names provided in the
            rxn_param_dct exist in the spc_dct
        """
        all_mech_names = ()
        for _rct_names, _prd_names in zip(rct_names, prd_names):
            all_mech_names += _rct_names
            all_mech_names += _prd_names
        all_mech_names = set(all_mech_names)

        missing_names = all_mech_names - all_spc_names
        if missing_names:
            print('Names in provided in mechanism, '
                  'but not provided in species list (likely from .csv file):')
            for name in missing_names:
                print('  ', name)
            print('Unable to finish parsing mechanism. Exiting...')
            sys.exit()

    def _inf(rct_names, prd_names, ich_dct):
        """ Sort reactant and product name lists by formula to facilitate
            multichannel, multiwell rate evaluations
        """
        rxn_name_lst, formula_str_lst, formula_dct_lst = [], [], []
        for _rct_names, _prd_names in zip(rct_names, prd_names):
            rxn_name = '='.join(['+'.join(_rct_names), '+'.join(_prd_names)])
            rxn_name_lst.append(rxn_name)
            rct_ichs = list(map(ich_dct.__getitem__, _rct_names))
            formula_dct, formula_str = get_fml(rct_ichs)
            formula_dct_lst.append(formula_dct)
            formula_str_lst.append(formula_str)

        return formula_dct_lst, formula_str_lst, rxn_name_lst

    # Extract info from dictionary
    rcts, prds, thrdbdy = zip(*rxn_param_dct.keys())
    rct_names, prd_names, thrdbdy_lst = list(rcts), list(prds), list(thrdbdy)

    # Check if the rxn_param dct may be fully parsed
    _check_names(rct_names, prd_names, set(spc_dct.keys()))

    # formulas and reaction names (repplace with the mech info from ckin
    ich_dct = name_inchi_dct(spc_dct)
    formula_dct, formula_str, rxn_name = _inf(rct_names, prd_names, ich_dct)

    return [formula_dct, formula_str,
            rct_names, prd_names, thrdbdy_lst,
            rxn_name, list(rxn_param_dct.values())]
