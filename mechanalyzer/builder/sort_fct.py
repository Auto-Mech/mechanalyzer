"""
Sorter module - sorting of the mechanism according to various options
- reaction classes
- pes/subpes
- multiplicity
- species subsets
- submechanism
"""
import enum
from multiprocessing.sharedctypes import Value
from re import T
import time
import sys
import copy
import pandas as pd
import numpy
import automol
from mechanalyzer.builder import submech
from mechanalyzer.builder import rxnclass
from mechanalyzer.calculator import rates as calc_rates
from mechanalyzer.calculator import thermo
from mechanalyzer.calculator import ene_partition
from mechanalyzer.calculator import nonboltz
from mechanalyzer.calculator import ktp_util
from mechanalyzer.parser import pes
from mechanalyzer.parser.spc import name_inchi_dct
from mechanalyzer.parser._util import count_atoms
from mechanalyzer.parser._util import order_rct_bystoich
from mechanalyzer.parser._util import extract_spc
from mechanalyzer.parser._util import get_mult
from mechanalyzer.parser._util import get_fml


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

    if not all([isinstance(val, dict) or isinstance(val, list) for val in rxn_param_dct.values()]):
        print(
            '*Warning: ktp dct vals for sorting purposes derived at [300:10:2500] K at 1 atm')
        rxn_ktp_dct = calc_rates.eval_rxn_param_dct(rxn_param_dct, [numpy.arange(300, 2010, 10)],
                                                    [1])
        for key, val in rxn_ktp_dct.items():
            if isinstance(val, dict):
                rxn_ktp_dct[key] = [val]
    else:
        rxn_ktp_dct = rxn_param_dct  # it means you already provided a ktp dct as input

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
            rxn_name, list(rxn_ktp_dct.values())]


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
                        'submech_del',  'rxn_max_vals', 'rxn_max_ratio',
                        ]
        labels_all = ['NR', 'N_of_prods', 'SPECIES',
                      'pes', 'subpes', 'channel',
                      'Heavier_rct', 'Totalmultiplicity',
                      'rxntype', 'rxnclass',
                      'SUBMECH', 'submech_prompt', 'submech_ext',
                      'submech_del', 'maxval', 'maxratio',
                      ]

        # check on pes/subpes criteria:
        # if subpes alone, also add pes
        if (('subpes' in hierarchy and 'pes' not in hierarchy) or
                ('chnl' in hierarchy and 'pes' not in hierarchy and
                 'subpes' not in hierarchy)):
            try:
                idx_insert_pes = hierarchy.index('subpes')
            except ValueError:
                idx_insert_pes = hierarchy.index('chnl')
            hierarchy.insert(idx_insert_pes, 'pes')
            # chenge N of headers if necessary
            hierarchy[-1] += 1*(hierarchy[-1] > idx_insert_pes)
            print('pes criterion added: necessary for subpes/chnls')

        # if species_list is not empty: pre-process the mechanism
        self.preproc_specieslist(species_list, hierarchy)

        # sort and label
        self.sort_and_label(criteria_all, labels_all)
        self.criteria_all = criteria_all
        self.labels_all = labels_all

    def sort_and_label(self, criteria_all, labels_all):
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
            'submech_del': self.group_submech,
            'rxn_max_vals': self.rxn_max_vals,
            'rxn_max_ratio': self.rxn_max_ratio
        }

        for optn, fun_name in sort_optns_dct.items():
            # Call sorting function
            if any(optn == inp_crt for inp_crt in self.hierarchy):
                # Generate corresponding dataframe
                df_optn = pd.DataFrame(
                    index=self.mech_df.index, columns=[optn])
                df_optn = fun_name(df_optn)
                # Concatenate the new portion of dataframe
                self.mech_df = pd.concat([self.mech_df, df_optn], axis=1)

        # drop ''
        if '' in self.mech_df.columns:
            self.mech_df = self.mech_df.drop([''], axis=1)
        # 1. Sort
        # series: ascending/descending values
        asc_val = [True]*len(criteria_all + ['rxn_names'])
        # rxn vals, ratio and names are descending
        asc_val[-3:] = [False, False, False]
        asc_series = pd.Series(asc_val, index=criteria_all + ['rxn_names'])

        try:
            # last "standard" criterion: rxn name
            asclst = list(
                asc_series[self.hierarchy[:-1] + ['rxn_names']].values)
            self.mech_df = self.mech_df.sort_values(
                by=self.hierarchy[:-1] + ['rxn_names'],
                ascending=asclst)
        except KeyError as err:
            print(
                'Error: Reactions not sorted according ',
                f'to all criteria: missing {err}, exiting')
            sys.exit()

        # 2. assign class headers
        labels = pd.Series(labels_all, index=criteria_all)
        self.class_headers(self.hierarchy, labels)

    def preproc_specieslist(self, species_list, hierarchy):

        if len(species_list) == 0 and 'submech_prompt' in hierarchy:
            # species list includes all radicals in the mech
            print('Prompt selected w/o species specification: \
                all radicals analyzed ...')
            for sp_i in self.spc_dct.keys():
                mult = get_mult(sp_i, self.spc_dct)
                atoms = sum(automol.chi.formula(
                    self.spc_dct[sp_i]['inchi']).values())
                if mult > 1 and atoms > 2:
                    species_list.append(sp_i)

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

            elif len(species_list) == 1 and 'submech_del' in hierarchy:
                species_list, species_subset_df = submech.species_subset_del(
                    species_list[0], self.spc_dct)
                self.species_subset_df = species_subset_df
                filtertype = 'submech_ext'
                
            elif len(species_list) > 1 and 'submech' in hierarchy:
                print('Error: pyr/ox submech extraction available ',
                      'for only 1 species')
                sys.exit()

            else:
                filtertype = 'submech'

            self.mech_df_full = copy.deepcopy(self.mech_df)
            self.spc_dct_full = copy.deepcopy(self.spc_dct)

            self.mech_df, self.spc_dct = self.filter_byspecies(
                species_list, filtertype)
            self.species_list = species_list

        self.hierarchy = hierarchy

    def filter_byspecies(self, species_list, filtertype):
        """ Find all reactions involving species of the species_list given as input
            Provides a new mechanism of the reactions of the selected subset
            and a reduced list of species

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param species_list: list of species subsets
        :param filtertype: can be 'submech', 'submech_ext/del', 'submech_prompt'

        :returns: mech_df, spc_dct
        :rtype: dataframe, dict
        """

        # add classification by subpes first
        if filtertype == 'submech_prompt':
            if 'submech_prompt' not in self.mech_df_full.columns:
                self.mech_df_full = pd.concat(
                    [self.mech_df_full, self.chnl('')], axis=1)  # add subpes
                self.mech_df_full[['submech_prompt', 'rxn_ped']] = ''
        # deepcopy
        mech_df = copy.deepcopy(self.mech_df_full)
        # check that all species selected are in the species dictionary
        if any(i not in self.spc_dct_full.keys() for i in species_list):
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

            elif filtertype =='submech_ext':
                # filter out bimol/bimol reactions if some bimol rcts/prds are not in the species list
                # equivalent to saying: at least all reactants (also single react works) or all products
                # must be in the species list
                chk = 0
                chk += int(all(any(rct == _spc for _spc in species_list)
                               for rct in rcts))
                chk += int(all(any(prd == _spc for _spc in species_list)
                               for prd in prds))*int(len(prds) <= 2)
                """
                chk = int(all(any(rct == _spc for _spc in species_list)
                               for rct in rcts))*int(all(any(prd == _spc for _spc in species_list)
                        for prd in prds))
                """
            if filtertype == 'submech_prompt' and chk >= 1:
                for sp in species_list:
                    if len(rcts) == 2 and any(rct == sp for rct in rcts) and len(prds) <= 2:
                        mech_df['submech_prompt'][rxn] = 'RAD_GEN_{}'.format(
                            sp)
                        mech_df['rxn_ped'][rxn] = '{}={}'.format(
                            rxn[0].split('=')[-1], rxn[0].split('=')[0])
                    elif len(prds) == 2 and any(pr == sp for pr in prds):
                        mech_df['submech_prompt'][rxn] = 'RAD_GEN_{}'.format(
                            sp)
                        mech_df['rxn_ped'][rxn] = rxn[0]

                    elif len(prds) > 2 and any(pr == sp for pr in prds):
                        mech_df['submech_prompt'][rxn] = 'PROMPT_LUMPED_{}'.format(
                            sp)
                    # rxn is unimol deco/formation of the radical
                    elif ((len(rcts) == 1 and rcts[0] == sp and len(prds) == 2) or
                          (len(prds) == 1 and prds[0] == sp and len(rcts) == 2)):
                        mech_df['submech_prompt'][rxn] = 'RAD_DECO_{}'.format(
                            sp)
                    elif (len(rcts) == 1 and rcts[0] == sp and len(prds) > 2):
                        mech_df['submech_prompt'][rxn] = 'RAD_DECO_LUMPED_{}'.format(
                            sp)

                continue  # next for loop

            if chk >= 1:  # keep reaction and add to species list
                spc_list.extend(rcts)
                spc_list.extend(prds)
            elif chk == 0:  # reaction filtered out
                mech_df = mech_df.drop(index=[rxn])

        if filtertype == 'submech_prompt':
            # submech_prompt: if RAD_GEN/RAD_DECO in list, keep the rxns
            for _, subpes_df in mech_df.groupby(['pes', 'subpes']):
                if any(len(CHECK.split('_')) > 1   # if len > 1, value was assigned
                       for CHECK in subpes_df['submech_prompt'].values):
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
        spc_dct_val = list(map(self.spc_dct_full.get, spc_list))
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
        # return nothing because already done in the species filtering,
        # nonsense to do it again
        submech_df = pd.DataFrame(
            index=self.mech_df.index, columns=[''])

        # assign group values if not present:
        # 1. groups by PES
        # 2. if PES is rad_deco: renames as rad_deco_radname
        # 3. if PES produces hot radical: adds it to a group ['grp N']
        #       then specifies the name in ['groupname'] nb add to same list if same subpes
        try:
            self.grps
        except AttributeError:
            grps = []
            self.species_deco_dct = dict.fromkeys(self.species_list)
            # filter by RAD_DECO: save the pes/subpes value

            for hot_label in ['RAD_DECO_{}'.format(sp) for sp in self.species_list]:
                hot_df = self.mech_df[self.mech_df['submech_prompt'] == hot_label]
                sp = hot_label.split('_')[-1]
                self.species_deco_dct[sp] = '{}:{}'.format(
                    hot_df['pes'].iloc[0], hot_df['subpes'].iloc[0])

            grp_dct_template = {'grp': 0, 'idxs': [], 'peds': [], 'hot': []}
            grpN = 0
            rad_gen_list = ['RAD_GEN_{}'.format(sp) for sp in self.species_list]
            ped_df = self.mech_df[self.mech_df['submech_prompt'].isin(
                rad_gen_list)]

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
                        peds.append(subpesdf['rxn_ped'][rxn])
                        # deal with hotspecies
                        sp = subpesdf['submech_prompt'][rxn].split('_')[-1]
                        hotpes = self.species_deco_dct[sp]
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
                    grp_dct['modeltype'] = 'rovib_dos'

                grps.append(grp_dct)

            self.grps = grps
        return submech_df

    def filter_groups_prompt(self, therm_dct, DFG, T0=300.):
        """ filter prompt groups according to DFG dct criteria:
            A+B->B+C->B+(D+E) must be < X kcal/mol endothermic
            If DH1 is added to dissociation reaction: k(T*)/k(T) > RATIO
            ratio: (DH1+DH3)/DH3 < threshold
            Absolute threshold for k(T*): should be > k*abs (only for exothermic rxns not fulfilling other criteria)
            NB all criteria should be strict --> by default, no reaction is kept.
            the criterion introduced should be the strict one
            T0 is the T at which DH are considered
            Tref is the T at which the decomposition rate of the radical is extracted
            otherwise remove from dct
        """

        # define criteria
        self.DHmax = DFG['DH']
        self.H5H3ratio = DFG['H5H3ratio']
        self.kratio = DFG['kratio']
        self.kabs = DFG['kabs']
        self.Tref = DFG['Tref']
        self.keepfiltered = int(DFG['keepfiltered'])
        # convert therm dct to dataframe
        self.therm_df = thermo.spc_therm_dct_df(therm_dct)
        # 1 find min endothermicity of deco reactions of hotspecies
        self.dh_min_hot = dict.fromkeys(self.species_list)
        self.k_max_hot = dict.fromkeys(self.species_list)
        self.labels_hot = dict.fromkeys(self.species_list)
        hot_sp_df_dct = dict.fromkeys(self.species_list)
        for hot_sp in self.species_list:
            hot_sp_df_dct[hot_sp] = self.mech_df[self.mech_df['submech_prompt']
                                                 == 'RAD_DECO_{}'.format(hot_sp)]
            self.k_max_hot[hot_sp], self.dh_min_hot[hot_sp], self.labels_hot[hot_sp] = nonboltz.get_max_reactivity(
                hot_sp, hot_sp_df_dct[hot_sp], self.therm_df, T0, self.Tref)

        # loop over groups and delete too endothermic reactions
        filtered_grps = []
        grp_idx = 0
        fmt_lbls = numpy.array(['%40s', '%15s', '%.1f', '%.1f',
                                '%.1f', '%1.1e', '%1.1e', '%.0f', '%s'], dtype=object)
        lbls = numpy.array(['prompt rxn{}\t'.format(' '*40), 'hot species', 'dh1({:.0f}K)*phi'.format(T0), 'dh2({:.0f}K)'.format(T0),
                            'dhtot({:.0f}K)'.format(T0), 'k({:.0f}K)'.format(self.Tref), 'k*(T*({:.0f}K))'.format(self.Tref), 'T*({:.0f}K)'.format(self.Tref), 'keep?'], dtype=object)

        self.rxns_dh = numpy.vstack((fmt_lbls, lbls))
        for grp in self.grps:
            grp_new = {'grp': 0, 'idxs': [],
                       'peds': [], 'hot': [], 'modeltype': ''}
            self.grp_add = {'idxs': [], 'peds': [], 'hot': []}
            check = 0
            check_filtered = 0
            active_hotsp = []
            # lists, potentially more than 1

            for n, ped in enumerate(grp['peds']):
                for ped_i in ped:
                    # extract rxn exo/endo thermicity
                    print('anaylze ped .. {}'.format(ped_i))
                    rcts = ped_i.split('=')[0].split('+')
                    prds = ped_i.split('=')[1].split('+')
                    # check with threshold
                    hot_sp = list(set(self.species_list) & set(prds))[0]
                    nonhot = list(set([hot_sp]) ^ set(prds))[0]

                    phi = ene_partition.phi_equip_fromdct(
                        hot_sp, nonhot, self.spc_dct)

                    dh = thermo.extract_deltaX_therm(
                        self.therm_df, rcts, prds, 'H')/1000

                    # calculate equivalent T* at Tref and corresponding k* rate
                    # phi: fraction of en transferred to prods - very approixmate
                    try:
                        T_star, k_star, dh_tot = nonboltz.estimate_hot_hk(
                            dh*phi*1000, self.Tref, self.therm_df[hot_sp]['Cp'], self.k_max_hot[hot_sp], self.dh_min_hot[hot_sp]*1000)
                    except TypeError:
                        print('dh failed for:', rcts, prds, hot_sp)
                        continue

                    if dh_tot[T0] < self.DHmax or dh_tot[T0]/self.dh_min_hot[hot_sp][T0] < self.H5H3ratio \
                        or (k_star/self.k_max_hot[hot_sp][self.Tref] > self.kratio) \
                            or (k_star > self.kabs and dh[T0] < 0):
                        check += 1
                        keep = 'YES'
                        active_hotsp.append(hot_sp)
                        # BUILD REACTION CHAINS STARTING FROM HERE
                        # look for chains
                        self.rxn_chain_prompt2(
                            T0, dh_tot, hot_sp, hot_sp_df_dct)
                    else:
                        if self.keepfiltered == 0:
                            ped.remove(ped_i)
                        else:
                            active_hotsp.append(hot_sp)
                            check += 1
                            check_filtered += 1  # count filtered rxns
                        keep = 'NO'

                    array_info = numpy.array(
                        [ped_i, hot_sp, dh[T0]*phi, self.dh_min_hot[hot_sp][T0], dh_tot[T0], self.k_max_hot[hot_sp][self.Tref], k_star, T_star, keep], dtype=object)
                    self.rxns_dh = numpy.vstack((self.rxns_dh, array_info))

                if ped:
                    grp_new['idxs'].append(grp['idxs'][n])
                    grp_new['peds'].append(ped)
                    grp_new['hot'].append([])
                    # CHECK FOR NEW CHAINS AND UPDATE GRP NEW
                    # NB REMOVE hotsp FROM ACTIVE_HOTSP TO AVOID DOUBLE COUNTING

            for n, hot in enumerate(grp['hot']):
                grp['hot'][n] = [hot_i for hot_i in hot if hot_i in active_hotsp]
                grp['hot'][n].sort()  # keep same order
                if grp['hot'][n]:
                    grp_new['idxs'].append(grp['idxs'][n])
                    # APPEND ALSO HOT SPECIES IF NEW CHAINS CONSIDERED
                    grp_new['peds'].append([])
                    grp_new['hot'].append(grp['hot'][n])

            # these are single indices
            for n, idx in enumerate(self.grp_add['idxs']):
                # merge other groups -
                # if you have same idxs: it is ped + hot, so take the
                # hot of the grp_new (more complete) and the ped label of the latter
                try:
                    i_idx = grp_new['idxs'].index(idx)
                    grp_new['peds'][i_idx] = self.grp_add['peds'][n]
                except ValueError:
                    grp_new['idxs'].append(idx)
                    # APPEND ALSO HOT SPECIES IF NEW CHAINS CONSIDERED
                    grp_new['peds'].append(self.grp_add['peds'][n])
                    grp_new['hot'].append(self.grp_add['hot'][n])

            if check > 0:
                grp_idx += 1
                grp_new['grp'] = grp_idx
                if check_filtered == check:
                    grp_new['modeltype'] = 'thermal'
                else:
                    grp_new['modeltype'] = 'rovib_dos'
                filtered_grps.append(grp_new)

        self.grps = filtered_grps
        # resort because you added reactions
        self.sort_and_label(self.criteria_all, self.labels_all)

    def rxn_chain_prompt(self, T0, dh_tot, rad, sp_df_dct):
        """ from a given hot product (rad), derive prompt reaction chain complying with thresholds
        """
        # append to self groups only if pes not present as ped
        if self.species_deco_dct[rad] not in self.grp_add['idxs']:
            self.grp_add['idxs'].append(self.species_deco_dct[rad])
            self.grp_add['peds'].append([])
            self.grp_add['hot'].append([rad])
            
        check_break = 0
        while check_break == 0:
            # radical generating reaction (and hot sp)

            dh_start = copy.deepcopy(dh_tot)
            # rad is the potentially prompt species:
            # 1. CHECK decomposition channels:
            # 1.1 DECO CHANNEL: IF IT PRODUCES RADICALS, CONTINUE
            prds = self.labels_hot[rad].split('=')[1].split('+')
            # potentially change the algorithm:
            # if you have 2 products, just analyze both; select the one
            # most 'compliant' with the conditions (e.g., smallest dhtot or ratio, largest kstar or k ratio)
            # and go on with that.

            if len(prds) == 2 and get_mult(prds, self.spc_dct) in [2, 3]:
                # 1 radical
                hot_sp = [
                    prd for prd in prds if self.spc_dct[prd]['mult'] > 1][0]
                phi = ene_partition.phi_equip_fromdct(
                    hot_sp, list(set([hot_sp]) ^ set(prds))[0], self.spc_dct)
            elif len(prds) == 2 and get_mult(prds, self.spc_dct) > 3:
                # 1.2 IF 2 RADICALS: SELECT THE ONE WITH LARGER PHI
                phis = [ene_partition.phi_equip_fromdct(
                        prd, list(set([prd]) ^ set(prds))[0], self.spc_dct)
                        for prd in prds]
                phi = max(phis)
                hot_sp = prds[phis.index(phi)]
            else:
                break
            # formula check
            if sum(automol.chi.formula(
                    self.spc_dct[hot_sp]['inchi']).values()) < 3:
                break
            
            print('prompt chain for {}'.format(self.labels_hot[rad]), \
                'hot species is {}'.format(hot_sp),
                'fraction of en transferred: {:.1f}'.format(phi))
            # 2. NEW HOT SP: NEW PROMPT
            # 2.1 SEARCH FOR DECO OF THE NEW HOT SP
            if hot_sp not in self.k_max_hot.keys():
                hot_mech_df, hot_spc_dct = self.filter_byspecies(
                    [hot_sp], 'submech_prompt')
                sp_df_dct[hot_sp] = hot_mech_df[hot_mech_df['submech_prompt']
                                                == 'RAD_DECO_{}'.format(hot_sp)]
                self.species_deco_dct[hot_sp] = '{}:{}'.format(
                    sp_df_dct[hot_sp]['pes'].iloc[0], sp_df_dct[hot_sp]['subpes'].iloc[0])
                self.k_max_hot[hot_sp], self.dh_min_hot[hot_sp], self.labels_hot[hot_sp] = nonboltz.get_max_reactivity(
                    hot_sp, sp_df_dct[hot_sp], self.therm_df, T0, self.Tref)  # high T to get bimol faster
                # update mechanism
                self.mech_df = pd.concat(
                    [self.mech_df, sp_df_dct[hot_sp]], axis=0)
                self.spc_dct.update(hot_spc_dct)
            # 2.2 IF COMPLIANT WITH CONDITIONS:
            try:
                T_star, k_star, dh_tot = nonboltz.estimate_hot_hk(
                    dh_start*phi*1000, self.Tref, self.therm_df[hot_sp]['Cp'], self.k_max_hot[hot_sp], self.dh_min_hot[hot_sp]*1000)
            except TypeError:
                print('dh failed for: {}'.format(self.labels_hot[rad]))
                break

            if dh_tot[T0] < self.DHmax or dh_tot[T0]/self.dh_min_hot[hot_sp][T0] < self.H5H3ratio \
                or (k_star/self.k_max_hot[hot_sp][self.Tref] > self.kratio) \
                    or (k_star > self.kabs and dh_start[T0] < 0):
                keep = 'YES'
                array_info = numpy.array(
                    [self.labels_hot[rad], hot_sp, dh_start[T0]*phi, self.dh_min_hot[hot_sp][T0], dh_tot[T0], self.k_max_hot[hot_sp][self.Tref], k_star, T_star, keep], dtype=object)
                self.rxns_dh = numpy.vstack((self.rxns_dh, array_info))
                peds = [self.labels_hot[rad]]
                
            else:
                # keep species only as hot
                check_break = 1
                peds = []

            if peds != []:
                # find index and replace ped
                idx = self.grp_add['idxs'].index(self.species_deco_dct[rad])
                self.grp_add['peds'][idx] = peds
                
                if  self.species_deco_dct[hot_sp] not in self.grp_add['idxs']:
                    self.grp_add['idxs'].append(self.species_deco_dct[hot_sp])
                    self.grp_add['peds'].append([])
                    self.grp_add['hot'].append([hot_sp])
           
            # updates for next cycle
            #     - new cycle: rad = hot_sp
            # UPDATE SELF.MECH_DF AND SPECIES TO ADD NEW PESs
                
            rad = copy.deepcopy(hot_sp)
            
    def rxn_chain_prompt2(self, T0, dh_tot, rad, sp_df_dct):
        """ from a given hot product (rad), derive prompt reaction chain complying with thresholds
        """
        # append to self groups only if pes not present as ped
        if self.species_deco_dct[rad] not in self.grp_add['idxs']:
            self.grp_add['idxs'].append(self.species_deco_dct[rad])
            self.grp_add['peds'].append([])
            self.grp_add['hot'].append([rad])
            
        check_break = 0
        while check_break == 0:
            # radical generating reaction (and hot sp)

            dh_start = copy.deepcopy(dh_tot)
            # rad is the potentially prompt species:
            # 1. CHECK decomposition channels:
            # 1.1 DECO CHANNEL: IF IT PRODUCES RADICALS, CONTINUE
            prds = self.labels_hot[rad].split('=')[1].split('+')
            # potentially change the algorithm:
            # if you have 2 products, just analyze both; select the one
            # most 'compliant' with the conditions (e.g., smallest dhtot ratio)
            # and go on with that.

            if len(prds) != 2:
                continue
            
            phis = dict.fromkeys(prds)
            T_stars = dict.fromkeys(prds)
            k_stars = dict.fromkeys(prds)
            dh_tots = dict.fromkeys(prds)
            add_check = dict.fromkeys(prds)
            for prd in prds:
                if sum(automol.chi.formula(
                        self.spc_dct[prd]['inchi']).values()) < 3:
                    continue
                phis[prd] = ene_partition.phi_equip_fromdct(
                    prd, list(set([prd]) ^ set(prds))[0], self.spc_dct)
                
                if prd not in self.k_max_hot.keys():
                    add_check[prd] = 1
                    hot_mech_df, hot_spc_dct = self.filter_byspecies(
                        [prd], 'submech_prompt')
                    sp_df_dct[prd] = hot_mech_df[hot_mech_df['submech_prompt']
                                                    == 'RAD_DECO_{}'.format(prd)]
                    self.species_deco_dct[prd] = '{}:{}'.format(
                        sp_df_dct[prd]['pes'].iloc[0], sp_df_dct[prd]['subpes'].iloc[0])
                    self.k_max_hot[prd], self.dh_min_hot[prd], self.labels_hot[prd] = nonboltz.get_max_reactivity(
                        prd, sp_df_dct[prd], self.therm_df, T0, self.Tref)  # high T to get bimol faster
                else:
                    add_check[prd] = 0
                try:
                    T_stars[prd], k_stars[prd], dh_tots[prd] = nonboltz.estimate_hot_hk(
                        dh_start*phis[prd]*1000, self.Tref, self.therm_df[prd]['Cp'], self.k_max_hot[prd], self.dh_min_hot[prd]*1000)
                except TypeError:
                    print('dh failed for: {}'.format(self.labels_hot[rad]))
                    break
                
            # identify hot species by min dh_tot
            if not any(phis.values()):
                break
            elif not all(phis.values()):
                hot_sp = [key for key in phis.keys() if phis[key]][0]
            else:
                dh_tots_T0 = [dh_tots[prds[0]][T0], dh_tots[prds[1]][T0]]
                min_dhtot = min(dh_tots_T0)
                hot_sp = prds[dh_tots_T0.index(min_dhtot)]
            dh_tot = dh_tots[hot_sp]
            print('prompt chain for {}'.format(self.labels_hot[rad]), \
                'hot species is {}'.format(hot_sp),
                'fraction of en transferred: {:.1f}'.format(phis[hot_sp]))
            
            ##################################################
            if add_check[hot_sp] == 1:
            # 2.1 SEARCH FOR DECO OF THE NEW HOT SP
                # update mechanism
                self.mech_df = pd.concat(
                    [self.mech_df, sp_df_dct[hot_sp]], axis=0)
                self.spc_dct.update(hot_spc_dct)

            # 2.2 IF COMPLIANT WITH CONDITIONS:

            if dh_tot[T0] < self.DHmax or dh_tot[T0]/self.dh_min_hot[hot_sp][T0] < self.H5H3ratio \
                or (k_stars[hot_sp]/self.k_max_hot[hot_sp][self.Tref] > self.kratio) \
                    or (k_stars[hot_sp] > self.kabs and dh_start[T0] < 0):
                keep = 'YES'
                array_info = numpy.array(
                    [self.labels_hot[rad], hot_sp, dh_start[T0]*phis[hot_sp], self.dh_min_hot[hot_sp][T0], dh_tot[T0], self.k_max_hot[hot_sp][self.Tref], k_stars[hot_sp], T_stars[hot_sp], keep], dtype=object)
                self.rxns_dh = numpy.vstack((self.rxns_dh, array_info))
                peds = [self.labels_hot[rad]]
                
            else:
                # keep species only as hot
                check_break = 1
                peds = []

            if peds != []:
                # find index and replace ped
                idx = self.grp_add['idxs'].index(self.species_deco_dct[rad])
                self.grp_add['peds'][idx] = peds
                
                if  self.species_deco_dct[hot_sp] not in self.grp_add['idxs']:
                    self.grp_add['idxs'].append(self.species_deco_dct[hot_sp])
                    self.grp_add['peds'].append([])
                    self.grp_add['hot'].append([hot_sp])
           
            # updates for next cycle
            #     - new cycle: rad = hot_sp
            # UPDATE SELF.MECH_DF AND SPECIES TO ADD NEW PESs
                
            rad = copy.deepcopy(hot_sp)
            
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

                    rclass = rxnclass.classify_graph(
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

                    rxn_type_ws = rxnclass.classify_ws(
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
            max_val = ktp_util.get_max_aligned_values(param_vals_dct)
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
            param_ratio_dct = ktp_util.get_aligned_rxn_ratio_dct(
                param_vals_dct)
            max_val = ktp_util.get_max_aligned_values(param_ratio_dct)
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
        # if comments already in mech_df: remove them:
        if 'cmts_inline' in self.mech_df.columns and 'cmts_top' in self.mech_df.columns:
            self.mech_df = self.mech_df.drop(
                ['cmts_inline', 'cmts_top'], axis=1)
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
