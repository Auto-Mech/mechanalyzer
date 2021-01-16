"""
Read the mechanism file
"""

import automol

from automol.graph._graph import explicit
from mechanalyzer.parser import ckin_ as ckin
from mechanalyzer.parser import submech as submech
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
        # set dataframe: extract useful info
        molecularity = list(map(len,rct_names_lst))
        N_of_prods = list(map(len,prd_names_lst))
        rct_names_lst_ordered = order_rct_bystoich(rct_names_lst,spc_dct=spc_dct) # put heavier reactant first
        prd_names_lst_ordered = order_rct_bystoich(prd_names_lst,spc_dct=spc_dct) # put heavier product first
        R1,R2 = get_S1S2(rct_names_lst_ordered) 
        numC,numN = count_atoms(formulas_dct)
        data = np.array([rct_names_lst,prd_names_lst,rct_names_lst_ordered,prd_names_lst_ordered,R1,R2,molecularity,N_of_prods,formula_str_lst],dtype=object).T
        self.mech_df = pd.DataFrame(data,index=rxn_name_lst,columns=['rct_names_lst','prd_names_lst','rct_names_lst_ord','prd_names_lst_ord','R1','R2','molecularity','N_of_prods','PES'])
        self.spc_dct = spc_dct # set for later use
        

    def sort(self,hierarchy,species_list):
        '''
        hierarchy = list of hierarchical criteria for the mech organization
        sorts the mechanism according to the given criteria
        species_list = list of species you want to isolate (empty if you want to process the full mech). if you want to extract a fuel submech:  ['speciesname','SUBMECH']
        '''
        
        # if species_list is not empty: pre-process the mechanism
        if len(species_list) > 0:
            if species_list[-1] == 'SUBMECH':
                #select a subset of species appropriately according to the stoichiometries specified in submech.py
                species_list, species_subset_df = submech.species_subset(species_list[0],self.spc_dct)
                # check compatibility with selected hierarchy: cannot have headers
                if hierarchy[-1] != 0:
                    print('error: submech extraction already provides headers. incompatible with other headers. please select 0')
                
            self.mech_df,self.spc_dct = self.filter_byspecies(species_list)
            self.species_list = species_list

        # 0. look for keywords and sort accordingly
        # LIST OF AVAILABLE SORTING OPTIONS (BESIDES R1 AND PES, ALWAYS AVAILABLE)
        sort_optns_dct = {
            'SPECIES':self.group_species,
            'SUBPES':self.conn_chn,
            'MULT':self.reac_mult,
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
        criteria_all = ['molecularity','N_of_prods','SPECIES','PES','SUBPES','R1','MULT','RXN_CLASS_BROAD','RXN_CLASS_GRAPH']
        labels_all = ['molecularity','N_of_prods','SPECIES','PES','SUBPES','Heavier rct','Total multiplicity','rxn type broad','rxn type']
        labels = pd.Series(labels_all,index=criteria_all)
        self.class_headers(hierarchy,labels)

    def filter_byspecies(self,species_list):
        """
        Find all reactions involving species of the species_list given as input
        make another spc_dct containing only the species appearing in the reactions of interest
        """
        mech_df = copy.deepcopy(self.mech_df)
        # check that all species selected are in the species dictionary
        if any(i not in self.spc_dct.keys() for i in species_list):
            print('Error in ISOLATE_SPECIES: not all species are in the species list ')
            exit()

        spc_list = []
        # for all reactions in the dataframe: check if you have the species of the selected list. otherwise remove the reaction
        for ii in mech_df.index:
            rcts = list(mech_df['rct_names_lst'][ii])
            prds = list(mech_df['prd_names_lst'][ii])
            # check if one of the species in the list is among reactants or products of the reaction considered
            if (any(rct == species for species in species_list for rct in rcts) == False
                and any(prd == species for species in species_list for prd in prds) == False):
                mech_df = mech_df.drop([ii])
            else:
                # append all species to the list
                spc_list.extend(rcts)
                spc_list.extend(prds)
        
        # filter spc_list: unique elements
        spc_list = list(set(spc_list))
        # new spc_dct
        spc_dct_val = list(map(self.spc_dct.get,spc_list))
        spc_dct = dict(zip(spc_list,spc_dct_val))

        return mech_df,spc_dct


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
        Identify reaction multiplicity from spc_dct
        update column 'MULT' in reac_mult_df
        '''
        # assign multiplicity values to each reactant
        for rxn in reac_mult_df.index:
            mult = 1
            for Ri in self.mech_df['rct_names_lst'][rxn]:
                mult *= self.spc_dct[Ri]['mult']
            reac_mult_df['MULT'][rxn] = str(mult)

        return reac_mult_df

    def rxn_class_broad(self,rxn_clB_df):
        '''
        Identify the reaction class by the graph approach
        '''
        for rxn in rxn_clB_df.index:
            rcts = self.mech_df['rct_names_lst_ord'][rxn]
            prds = self.mech_df['prd_names_lst_ord'][rxn]
            if len(rcts) > 1:
                # bimolecular reaction classification
                rxn_class_broad = submech.classify_bimol(rcts,prds,self.spc_dct)
            else:
                # unimolecular reaction classification
                rxn_class_broad = submech.classify_unimol(rcts,prds,self.spc_dct)
            rxn_clB_df['RXN_CLASS_BROAD'][rxn] = rxn_class_broad

        return rxn_clB_df

    def rxn_class_graph(self,rxn_clG_df):
        '''
        Identify the reaction class by the graph approach
        1. group by subpes
        2. identify rxn in each subpes
        3. classify well skipping channels
        '''
        # 1. group by subpes
        # check that SUBPES is present in indexes, otherwise generate corresponding dataframe
        if 'SUBPES' not in self.mech_df.columns:
                df_optn = pd.DataFrame(index=self.mech_df.index,columns=['SUBPES'])
                df_optn = self.conn_chn(df_optn)
                # concatenate the new portion of dataframe
                self.mech_df = pd.concat([self.mech_df,df_optn],axis=1)

        # 2. graph classification or each subpes
        for subpes_name,subpes_df in self.mech_df.groupby('SUBPES'):
            # sort by molecularity: analyze first unimolecular isomerizations, unimo decompositions, and then bimolecular reactions
            subpes_df = subpes_df.sort_values(by=['molecularity','N_of_prods'])
            # REFER TO REORDERED SPECIES NAMES, OTHERWISE YOU MAY HAVE INCONSISTENT SPECIES NAMING
            # subpes species list
            rcts = list(subpes_df['rct_names_lst_ord'].values)
            prds = list(subpes_df['prd_names_lst_ord'].values)
            species_subpes = list(set(rcts+prds))

            #build dataframe: elementary reactivity matrix
            elem_reac_df = pd.DataFrame(np.zeros((len(species_subpes),len(species_subpes)),dtype='<U32'),index=species_subpes,columns=species_subpes)
            mult_species_subpes = pd.Series(list(map(len,species_subpes)),index=species_subpes)
            unimol_species = mult_species_subpes[mult_species_subpes==1].index
            # graph classification
            for rxn in subpes_df.index:
                rct_names = subpes_df['rct_names_lst'][rxn]
                prd_names = subpes_df['prd_names_lst'][rxn]
                rct_names_ord = subpes_df['rct_names_lst_ord'][rxn]
                prd_names_ord = subpes_df['prd_names_lst_ord'][rxn]
                # exclude all reactions with more than 2 reactants or products (not elementary!)
                if len(rct_names) < 3 and len(prd_names) < 3:
                    # Get the inchis and graphs
                    rct_ichs = list(self.spc_dct[rct]['inchi'] for rct in rct_names)
                    rct_graph = list(map(automol.inchi.graph, rct_ichs))
                    rct_gras = list(map(automol.graph.without_stereo_parities, rct_graph))
                    # print(automol.graph.string(rct_gra))  
                    prd_ichs = list(self.spc_dct[prd]['inchi'] for prd in prd_names)
                    prd_graph = list(map(automol.inchi.graph, prd_ichs))
                    prd_gras = list(map(automol.graph.without_stereo_parities, prd_graph))
                    # ID reaction
                    if automol.graph.reac.is_valid_reaction(rct_gras, prd_gras):
                        rclass = automol.graph.reac.classify_simple(rct_gras, prd_gras)
                    else:
                        rclass = 'unclassified - Wrong Stoichiometry'
                    # check stereo compatibility - I am not sure about this input
                    # ret = automol.graph.trans.is_stereo_compatible(rclass, rct_graph, prd_graph)
                    if rclass == None:
                        if subpes_df['molecularity'][rxn] == 1 and subpes_df['N_of_prods'][rxn] == 1:
                            rclass = 'isomerization' # TEMPORARY - BEFORE ALL ISOMERIZATION TYPES ARE ANALYZED
                        else:
                            rclass = 'unclassified'
                else:
                    rclass = 'unclassified - lumped'

                rxn_clG_df['RXN_CLASS_GRAPH'][rxn] = rclass

                # store values in the elementary reactivity matrix (for now contaminated with isomerizations)

                elem_reac_df[rct_names_ord][prd_names_ord] = rclass
                elem_reac_df[prd_names_ord][rct_names_ord] = rclass

            # 3. classify well skipping channels
            # reclassify the unclassified reactions A->B+C, B+C->D, B+C->E+F
            for rxn in subpes_df.index:
                if rxn_clG_df['RXN_CLASS_GRAPH'][rxn] == 'unclassified':

                    rct_names = subpes_df['rct_names_lst_ord'][rxn]
                    prd_names = subpes_df['prd_names_lst_ord'][rxn]   
                    # reactants: if bimolecular, find the label of the elementary reaction going to unimolecular species; if unimol, label is 'isom'
                    # isolate A+B->C and C->A+B connections
                    rxn_types = elem_reac_df[rct_names][unimol_species]
                    rxn_types = rxn_types[rxn_types != 'unclassified']
                    rxn_types_1 = rxn_types[rxn_types != '']

                    rxn_types = elem_reac_df[prd_names][unimol_species]
                    rxn_types = rxn_types[rxn_types != 'unclassified']
                    rxn_types_2 = rxn_types[rxn_types != '']

                    try:
                        rxn_type_1 = rxn_types_1[0] # TEMPORARY: SHOULD RECONSTRUCT THE FULL PATH FROM REACTANTS TO PRODUCTS
                        rxn_type_2 = rxn_types_2[0] # TEMPORARY: SHOULD RECONSTRUCT THE FULL PATH FROM REACTANTS TO PRODUCTS
                        rxn_clG_df['RXN_CLASS_GRAPH'][rxn] = rxn_type_1 + '-' + rxn_type_2 + ' (WS)'
                    except IndexError:
                        continue # do not change the reaction type
                        

        return rxn_clG_df

        
    ###################################### output dataframe ##############################
    def return_mech_df(self):
        '''
        Returns:
        - new_idx: reactants and products indices as tuples
        - cmts: dictionary containing comments of the corresponding reactions; indices are new_idx
        - self.spc_dct: species dictionary; may be different from the input if a subset of reactions is selected; useful for mech writing
        '''
        rct_names = self.mech_df['rct_names_lst'].values
        prd_names = self.mech_df['prd_names_lst'].values
        new_idx = list(zip(rct_names, prd_names))
        # store comments in dct
        cmts_df = pd.DataFrame(self.mech_df[['cmts_top','cmts_inline']].values,index=new_idx,columns=['cmts_top','cmts_inline'])
        cmts = cmts_df.to_dict('index')
        
        return new_idx,cmts,self.spc_dct

########################## useful functions run in the class #######################
def order_rct_bystoich(rct_names_lst,spc_dct=None):
    '''
    reorder reactants and products based on the higher number of atoms
    If no species dictionary is given as input: reorder just according to name length, 
    if length or stoichiometry is the same, by alphabetical order
    '''
    rct_names_lst_ordered = copy.deepcopy(rct_names_lst)
    ich_dct = {}
    if spc_dct:
        for key in spc_dct.keys():
            if 'ts' not in key and 'global' not in key:
                ich_dct[key] = spc_dct[key]['inchi']

        for key,val in enumerate(rct_names_lst_ordered):
            rct_names= val
            rct_ichs = list(map(ich_dct.__getitem__, rct_names))
            fml_rct = list(map(automol.inchi.formula,rct_ichs))
            atoms_rct = list(map(automol.formula.atom_count,fml_rct))
            if len(rct_names)==2:
                if atoms_rct[1] > atoms_rct[0]:
                    # swap places of reactants 1 and 2
                    rct_names_lst_ordered[key] = (rct_names[1],rct_names[0])
                elif atoms_rct[1] == atoms_rct[0]:
                    rct_names = list(rct_names)
                    rct_names.sort()
                    rct_names_lst_ordered[key] = tuple(rct_names)
                
    else:
        for key,val in enumerate(rct_names_lst_ordered):
            rct_names= val
            if len(rct_names)==2:
                if len(rct_names[1]) > len(rct_names[0]):
                    # swap places of reactants 1 and 2
                    rct_names_lst_ordered[key] = (rct_names[1],rct_names[0])  
                elif len(rct_names[1]) == len(rct_names[0]):
                    rct_names = list(rct_names)
                    rct_names.sort()
                    rct_names_lst_ordered[key] = tuple(rct_names)
    

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
        cmtlabel = cltype + ': ' + '_'.join(label) + '  '
        rxnclass = '! '+ cmtlabel +'_'.join(name)

    return rxnclass

def count_atoms(fml_list):
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
    conndct = {0: [['S1','S2'],['S3','S4'],[S5']], 1:[['S6','S7'],['S8','S9']]}
    connchnls = {0: [0,1,2] , 1:[3,4,5]}...
    corresponding to each subpes
    '''
    # preprocessing:
    # order (bimol) reactants and products in the same fashion
    # example if you have A+B and B+A they will be ordered in the same way
    pes_rct_names_lst = order_rct_bystoich(pes_rct_names_lst)
    pes_prd_names_lst = order_rct_bystoich(pes_prd_names_lst)
    # put everything in a dataframe. indices are the numbers given by enumerate
    len_rct_prd = np.array(list(map(len,pes_rct_names_lst)))+np.array(list(map(len,pes_prd_names_lst)))
    pes_df = pd.DataFrame(np.array([pes_rct_names_lst,pes_prd_names_lst,len_rct_prd],dtype=object).T,index=np.arange(0,len(pes_rxn_name_lst)),columns=['rcts','prds','N_rcts_prds'])
    # order according to the total number of species (N of reactants + N of products)
    pes_df = pes_df.sort_values(by='N_rcts_prds')
    # Split up channels into a connected sub-pes within a formula
    subpes_idx = 0
    conndct = {}
    connchnls = {}

    for chnl_idx in pes_df.index:
        connected_to = []
        chnl_species = [list(pes_df['rcts'][chnl_idx]),
                        list(pes_df['prds'][chnl_idx])]

        for conn_chnls_idx in conndct:
            for spc_pair in chnl_species:
                if len(spc_pair) == 1:
                    # this works for unimol species; need also to verify bimol wellskipping channels
                    if spc_pair in conndct[conn_chnls_idx]:
                        if conn_chnls_idx not in connected_to:
                            connected_to.append(conn_chnls_idx)
                    elif spc_pair[::-1] in conndct[conn_chnls_idx]:
                        if conn_chnls_idx not in connected_to:
                            connected_to.append(conn_chnls_idx)

            if len(chnl_species[0]) == 2 and len(chnl_species[1]) == 2:
                # bimol bimol reactions
                if (chnl_species[0] in conndct[conn_chnls_idx]) and (chnl_species[1] in conndct[conn_chnls_idx]):
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

