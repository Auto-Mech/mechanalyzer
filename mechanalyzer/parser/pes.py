"""
Read the mechanism file
"""

import automol
from automol.graph._graph import explicit
from mechanalyzer.parser import ckin_ as ckin
import pandas as pd
import numpy as np
import copy
import os
import ioformat

# PARSE THE MECHANISM FIle 
def read_mechanism_file(mech_str, mech_type, spc_dct, sort_rxns=False):
    """ Get the reactions and species from the mechanism input
    """

    # Parse the info from the chemkin file
    if mech_type == 'chemkin':
        formulas_dct, formulas, rct_names, prd_names, rxn_names = ckin.parse(
            mech_str, spc_dct,sort_rxns)
    else:
        raise NotImplementedError

    return [formulas_dct, formulas, rct_names, prd_names, rxn_names]
    # list, list of tuples, list of tuples, list


class SORT_MECH:
    '''
    class of functions to organize the mechanism according to given criteria
    from any step after initialization: call "return_mech_df" to get the current dataframe with mech info
    '''
    def __init__(self,mech_info,spc_dct):
        # extract data from mech info
        [formulas_dct,formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst] = mech_info         
        # set dataframe
        # index = rxn name
        # cols = reac1, reac2, prod1, prod2, formula
        print(rct_names_lst)
        rct_names_lst_ordered = order_rct_bystoich(rct_names_lst,spc_dct) # put heavier reactant first
        R1,R2 = get_S1S2(rct_names_lst_ordered) 
        numC,numN = count_C_N(formulas_dct)
        print(rct_names_lst)
        data = np.array([rct_names_lst,prd_names_lst,R1,R2,formula_str_lst,numC,numN],dtype=object).T
        self.mech_df = pd.DataFrame(data,index=rxn_name_lst,columns=['rct_names_lst','prd_names_lst','R1','R2','PES','numC','numN'])
        self.spc_dct = spc_dct # set for later use
        

    def sort(self,hierarchy,species_list):
        '''
        hierarchy = list of hierarchical criteria for the mech organization
        sorts the mechanism according to the given criteria
        species_list = list of species you want to isolate (empty if you want to process the full mech)
        '''
        
        # if species_list is not empty: pre-process the mechanism
        if len(species_list) > 0:
            self.filter_byspecies(species_list)
            self.species_list = species_list

        # 0. look for keywords and sort accordingly
        # LIST OF AVAILABLE SORTING OPTIONS (BESIDES R1 AND PES, ALWAYS AVAILABLE)
        sort_optns_dct = {
            'SPECIES':self.group_species,
            'SUBPES':self.conn_chn,
            'MULT_R1':self.reac_mult,
            'RXN_CLASS_BROAD':self.rxn_class_broad,
            'RXN_CLASS_GRAPH':self.rxn_class_graph}


        for optn,fun_name in sort_optns_dct.items():
            # call sorting function 
            if any(optn == inp_crt for inp_crt in hierarchy):
                # generate corresponding dataframe
                df_optn = pd.DataFrame(index=self.mech_df.index,columns=[optn])
                df_optn = fun_name(df_optn)
                # concatenate the new portion of dataframe
                self.mech_df = pd.concat([self.mech_df,df_optn],axis=1)

        # 1. sort
        try:
            self.mech_df = self.mech_df.sort_values(by=hierarchy[:-1])
        except KeyError as e:
            print('WARNING: The reactions were not sorted according to all desired criteria: missing {}'.format(e))
            return None

        # 2. assign class headers
        # set labels for all the possible criteria
        criteria_all = ['SPECIES','PES','SUBPES','numC','R1','MULT_R1','RXN_CLASS_BROAD','RXN_CLASS_GRAPH']
        labels_all = ['SPECIES','PES','SUBPES','N of C atoms','Heavier rct','Multiplicity of rct1','rxn type broad','rxn type']
        labels = pd.Series(labels_all,index=criteria_all)
        self.class_headers(hierarchy,labels)

    def filter_byspecies(self,species_list):
        """
        Find all reactions involving species of the species_list given as input
        """
        # check that all species selected are in the species dictionary
        if any(i not in self.spc_dct.keys() for i in species_list):
            print('Error in ISOLATE_SPECIES: not all species are in the species list ')
            exit()

        # for all reactions in the dataframe: check if you have the species of the selected list. otherwise remove the reaction
        for ii in self.mech_df.index:
            rcts = list(self.mech_df['rct_names_lst'][ii])
            prds = list(self.mech_df['prd_names_lst'][ii])
            # check if one of the species in the list is among reactants or products of the reaction considered
            if (any(rct == species for species in species_list for rct in rcts) == False
                and any(prd == species for species in species_list for prd in prds) == False):
                self.mech_df = self.mech_df.drop([ii])


    def class_headers(self,hierarchy,labels):
        """
        Read the hierarchy;
        assign classes based on the selected hierarchy
        adds columns in the df with appropriate comments/headers
        first N keywords --> main class --> comments_top (new col in df)
        other keywords --> other subclasses --> comments_inline
        """
        # df for comments_top and comments_inline
        ept_df = np.zeros((len(self.mech_df.index),1),dtype=str)
        df_cmts_top = pd.DataFrame(ept_df,index=self.mech_df.index,columns=['cmts_top'])
        df_cmts_inline = pd.DataFrame(ept_df,index=self.mech_df.index,columns=['cmts_inline'])

        try:
            N = int(hierarchy[-1])
        except ValueError:
            print('Last line of sorting options must be the N of criteria to be used for class headers')
            exit()
        ####### write topheader comments ######
        if N > 0:
            for name,rxndf in self.mech_df.groupby(hierarchy[:N]):
                # write rxn class as top header comments
                rxnclass = cmts_string(name,labels[hierarchy[:N]],'class_head')
                idx0 = rxndf.index[0]
                df_cmts_top['cmts_top'][idx0] = rxnclass

                ##### write inline comments if necessary #########
                if N < len(hierarchy)-1:
                    for name2,rxndf2 in rxndf.groupby(hierarchy[N:-1]):
                        rxnclass = cmts_string(name2,labels[hierarchy[N:-1]],'subclass')
                        df_cmts_inline['cmts_inline'][rxndf2.index] = rxnclass
        else:
            # write only inline comments
            for name,rxndf in self.mech_df.groupby(hierarchy[N:-1]):            
                rxnclass = cmts_string(name,labels[hierarchy[N:-1]],'class')
                df_cmts_inline['cmts_inline'][rxndf.index] = rxnclass

        # concatenate DFs
        self.mech_df = pd.concat([self.mech_df,df_cmts_top,df_cmts_inline],axis=1)

    def conn_chn(self,conn_chn_df):
        '''
        Identify connected channels
        Generate column 'SUBPES' in conn_chn_df
        '''

        for fml,peslist in self.mech_df.groupby('PES'):
            #print(peslist)
            # Set the names lists for the rxns and species needed below
            pes_rct_names_lst = peslist['rct_names_lst'].values
            pes_prd_names_lst = peslist['prd_names_lst'].values
            pes_rxn_name_lst = peslist.index
            connchnls = find_conn_chnls(pes_rct_names_lst,pes_prd_names_lst,pes_rxn_name_lst)
            # write subpes in conn_chn_df
            for key,value in connchnls.items():
                rxns = peslist.iloc[value].index
                conn_chn_df['SUBPES'][rxns] = fml + '-' + str(key)

        return conn_chn_df

    def group_species(self,reac_sp_df):
        """
        Creates a new df column "SPECIES" - recognizes the species you set in the list
        Also in this case the species_list is hierarchical: if the first group contains also a species of the second group,
        the reaction remains in the first group
        """
        # if species list is not found: do nothing - species entry will remain empty
        if len(self.species_list) > 0:
            for rxn in reac_sp_df.index:
                rcts = list(self.mech_df['rct_names_lst'][rxn])
                prds = list(self.mech_df['prd_names_lst'][rxn])
                # check species hierarchically             
                for sp in self.species_list:
                    if (any(sp==rct for rct in rcts) or any(sp==prd for prd in prds)) and type(reac_sp_df['SPECIES'][rxn])==float:
                        reac_sp_df['SPECIES'][rxn] = sp 
        return reac_sp_df

    def reac_mult(self,reac_mult_df):
        '''
        Identify reactant multiplicity from spc_dct
        update column 'MULT_R1' in reac_mult_df
        '''
        # assign multiplicity values to each reactant
        for rxn in reac_mult_df.index:
            R1 = self.mech_df['R1'][rxn]
            reac_mult_df['MULT_R1'][rxn] = str(self.spc_dct[R1]['mult'])

        return reac_mult_df

    def rxn_class_broad(self,rxn_clB_df):
        '''
        Identify the reaction class by the graph approach
        '''
        print('classification by broad reaction classes is not available yet')
        exit()

        return rxn_clB_df

    def rxn_class_graph(self,rxn_clG_df):
        '''
        Identify the reaction class by the graph approach
        '''
        for rxn in rxn_clG_df.index:
            rct_names = self.mech_df['rct_names_lst'][rxn]
            prd_names = self.mech_df['prd_names_lst'][rxn]

            # exclude all reactions with more than 2 reactants or products (not elementary!)
            if len(rct_names) < 3 and len(prd_names) < 3:
                # Get the inchis and graphs
                rct_ichs = list(self.spc_dct[rct]['inchi'] for rct in rct_names)
                rct_graph = graph_from_ichs(rct_ichs)
                # delete stereo information before the classification
                rct_gras = list(map(automol.graph.without_stereo_parities, rct_graph))
                # print(automol.graph.string(rct_gra)) - 
                prd_ichs = list(self.spc_dct[prd]['inchi'] for prd in prd_names)
                prd_graph = graph_from_ichs(prd_ichs)
                prd_gras = list(map(automol.graph.without_stereo_parities, prd_graph))

                # ID reaction
                try:
                    _, _, _, rclass = automol.graph.reac.classify(rct_gras, prd_gras)
                # check stereo compatibility - I am not sure about this input
                # ret = automol.graph.trans.is_stereo_compatible(rclass, rct_graph, prd_graph)
                except AssertionError:
                    rclass = 'unclassified - Assert Error'

                if rclass == None:
                    rclass = 'unclassified'
            else:
                rclass = 'unclassified - lumped'

            rxn_clG_df['RXN_CLASS_GRAPH'][rxn] = rclass


        return rxn_clG_df

        
    ###################################### output dataframe ##############################
    def return_mech_df(self):
        '''
        Returns the dataframe in the current status
        Indexes become the reactant and product tuples --> adapt to rxn_param_dct
        '''
        # print(self.mech_df)
        rct_names = self.mech_df['rct_names_lst'].values
        prd_names = self.mech_df['prd_names_lst'].values
        new_idx = list(zip(rct_names, prd_names))
        # store comments in dct
        cmts_df = pd.DataFrame(self.mech_df[['cmts_top','cmts_inline']].values,index=new_idx,columns=['cmts_top','cmts_inline'])
        cmts = cmts_df.to_dict('index')
        
        return new_idx,cmts

########################## useful functions run in the class #######################
def order_rct_bystoich(rct_names_lst,spc_dct):
    '''
    reorder reactants and products based on the higher number of atoms
    '''
    rct_names_lst_ordered = copy.deepcopy(rct_names_lst)
    ich_dct = {}
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich_dct[key] = spc_dct[key]['inchi']

    for key,val in enumerate(rct_names_lst_ordered):
        rct_names= val
        rct_ichs = list(map(ich_dct.__getitem__, rct_names))
        fml_rct = list(map(automol.inchi.formula,rct_ichs))
        atoms_rct = list(map(automol.formula.atom_count,fml_rct))
        if len(rct_names)>1:
            if atoms_rct[1] > atoms_rct[0]:
                # swap places of reactants 1 and 2
                rct_names_lst_ordered[key] = (rct_names[1],rct_names[0])

    return rct_names_lst_ordered
    
def cmts_string(name,label,cltype):
    '''
    Return appropriate comment string depending on the type
    name: rxn class
    cltype: class types. options: class_head, class, subclass
    label: class label
    '''
    # assing top headers:
    tophead = '!!!!!!!!! class !!!!!!!!!\n'
    bottomhead = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
    if isinstance(name,str):
        name = [name]
    elif isinstance(name,int):
        name = [str(name)]
    else:
        name = np.array(name,dtype=str)

    if cltype == 'class_head':
        cmtlabel = '!       '+'_'.join(label)+'\n'
        rxnlabel = '!       '+'_'.join(name)+'\n'
        rxnclass = tophead + cmtlabel + rxnlabel + bottomhead
    else:
        cmtlabel = cltype + ': ' + ' _ '.join(label) + '  '
        rxnclass = '! '+ cmtlabel +' _ '.join(name)

    return rxnclass

def count_C_N(fml_list):
    '''
    count C and N atoms in formula list
    should go somewhere else?
    '''
    count_C_lst = []
    count_N_lst = []
    for fml in fml_list:
        count_C = automol.formula.element_count(fml,'C')
        count_N = automol.formula.element_count(fml,'N')
        count_C_lst.append(count_C)
        count_N_lst.append(count_N)

    return count_C_lst, count_N_lst


def graph_from_ichs(ichs):
    '''
    Generate graphs to be processed by the rxn classifier without stereo parities involved
    Calls sets of functions of automol
    '''
    rxn_graph = list(map(automol.inchi.graph, ichs))
    rxn_graph = list(map(explicit, rxn_graph))
    # reorder the keys for multiple reactants
    if len(rxn_graph) == 2:
        rxn_graph, _ = automol.graph.standard_keys_for_sequence(rxn_graph)

    return rxn_graph


def get_S1S2(SPECIES):
    '''
    extract species 1 from tuple
    '''
    S1 = []
    S2 = []
    for S in SPECIES:
        if len(S) > 1:
            # bimol species
            S1.append(S[0])
            S2.append(S[1])
        else:
            # unimol species
            S1.append(S[0])
            S2.append('')

    return S1,S2

def find_conn_chnls(pes_rct_names_lst,pes_prd_names_lst,pes_rxn_name_lst):
    '''
    Given rxn names, reactants, products belonging to 1 PES:
    generate SUB PESs dictionaries
    {0: [0,1,2] , 1:[3,4,5]}...
    corresponding to each subpes
    '''
    # Split up channels into a connected sub-pes within a formula
    subpes_idx = 0
    conndct = {}
    connchnls = {}
    for chnl_idx, _ in enumerate(pes_rxn_name_lst):
        connected_to = []
        chnl_species = [list(pes_rct_names_lst[chnl_idx]),
                        list(pes_prd_names_lst[chnl_idx])]
        for conn_chnls_idx in conndct:
            for spc_pair in chnl_species:
                if len(spc_pair) == 1:
                    if spc_pair in conndct[conn_chnls_idx]:
                        if conn_chnls_idx not in connected_to:
                            connected_to.append(conn_chnls_idx)
                    elif spc_pair[::-1] in conndct[conn_chnls_idx]:
                        if conn_chnls_idx not in connected_to:
                            connected_to.append(conn_chnls_idx)
        if not connected_to:
            conndct[subpes_idx] = chnl_species
            connchnls[subpes_idx] = [chnl_idx]
            subpes_idx += 1
        else:
            conndct[connected_to[0]].extend(chnl_species)
            connchnls[connected_to[0]].append(chnl_idx)
            if len(connected_to) > 1:
                for cidx, cval in enumerate(connected_to):
                    if cidx > 0:
                        conn_specs = conndct.pop(cval, None)
                        conn_chnls = connchnls.pop(cval, None)
                        conndct[connected_to[0]].extend(conn_specs)
                        connchnls[connected_to[0]].extend(conn_chnls)
            for cidx in conndct:
                conndct[cidx].sort()
                conndct[cidx] = [
                    conndct[cidx][i] for i in
                    range(len(conndct[cidx])) if i == 0 or
                    conndct[cidx][i] != conndct[cidx][i-1]]

    return connchnls   

#################### functions working with dictionaries   
# FUNTIONS FOR THE PES DICT OBJECTS CONTAINING INFO FOR THE REACTIONS ON PES
def build_pes_dct(formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst):
    """ Build a dictionary of the PESs
    """

    pes_dct = {}
    current_formula = ''
    for fidx, formula in enumerate(formula_str_lst):
        if current_formula == formula:
            pes_dct[formula]['rct_names_lst'].append(rct_names_lst[fidx])
            pes_dct[formula]['prd_names_lst'].append(prd_names_lst[fidx])
            pes_dct[formula]['rxn_name_lst'].append(rxn_name_lst[fidx])
        else:
            current_formula = formula
            pes_dct[formula] = {}
            pes_dct[formula]['rct_names_lst'] = [rct_names_lst[fidx]]
            pes_dct[formula]['prd_names_lst'] = [prd_names_lst[fidx]]
            pes_dct[formula]['rxn_name_lst'] = [rxn_name_lst[fidx]]

    return pes_dct


# FUNCTIONS FOR THE CHANNELS DICT OBJECTS
def connected_channels_dct(pes_dct):
    """ Determine all the connected reaction channels for each PES
        Build a dictionary for each PES with lists of connected channels:
            dct[PES_FORMULA] = [ [SUB_PES_1], [SUB_PES_2], ... , [SUB_PES_N] ]
            where each SUB_PES = [n1, n2, ... , nN],
            where n1 to nN correspond to ixds for channels that are
            connected to each other
        For efficiency we only determine channels for PESs we wish to run.
    """
    conn_chn_dct = {}
    for _, formula in enumerate(pes_dct):
        # Set the names lists for the rxns and species needed below
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']

        connchnls = find_conn_chnls(pes_rct_names_lst,pes_prd_names_lst,pes_rxn_name_lst)

        # Add connected channels list to the dictionary
        conn_chn_dct[formula] = connchnls

    return conn_chn_dct

