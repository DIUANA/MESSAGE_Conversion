# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 15:08:53 2019

@author: Fabio Diuana
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:21:26 2019

@author: Fabio Diuana

"""
# =============================================================================
#packages required
from datetime import datetime
import os
import pandas as pd
import numpy as np
import re
from itertools import cycle
import itertools
import pickle
import sys

idx = pd.IndexSlice #pandas option to slice list
 
#functions that I am using in the main function
#sys.path.append('D:FabioDiuana_Pdrive_IIASAModelsNEST-BLUESInput_data_scripts')
#sys.path.append('D:/Fabio/Diuana_Pdrive_IIASA/Models/NEST-BLUES/Input_data_scripts')
#from xfuncs import *
#%%
def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]
#%%
from math import ceil, floor
def float_round(num, places = 0, direction = floor):
    return direction(num * (10**places)) / float(10**places)
#%%
def hasNumbers(inputString, any_element=False):
    if any_element:
        a = any(char.isdigit() for char in inputString)
    else:
        a = all(char.isdigit() for char in inputString)
    return a
#%%
def map_level(df, dct, level=0):
    index = df.index
    index.set_levels([[dct.get(item, item) for item in names] if i==level else names
                      for i, names in enumerate(index.levels)], inplace=True)
    
#%%
def find_between(s, first, last ): #function to find everything between two strings
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""    
    
#%%
def sort_list(list1, list2): 
  
    zipped_pairs = zip(list2, list1) 
  
    z = [x for _, x in sorted(zipped_pairs)] 
      
    return z   
#%%
def r_t_z(x):
    """
    removing trailing zeros
    """
    return str(x).rstrip('0').rstrip('.')

#remove_trailing_zeros(1.23)  # '1.23'
#%%
# =============================================================================
def expand_grid_name(name,df,l,*args):
    """
        Expand a dataframe according to the list content
        replicating every row times the the number
        of elements in the list
    """

    ncols = sum([c.shape[1] if type(c)==np.ndarray else 1 for c in [df,l,args]])
    if [isinstance(k, pd.DataFrame) for k in [df,l,args]].count(True) !=0:
        raise ValueError('pass the dataframe variable as value = df.values')
        if len(name) != ncols:
            raise ValueError('variable name must have the same number of parameters as the rest of the arguments in the function.')
    rows = itertools.product(df,l)        
    a = pd.DataFrame.from_records(rows, columns=['A_'+str(i) for i in range(2)])
    #b=pd.DataFrame(a['df'].tolist(), columns=df.columns.tolist())
    c = a
    if len(args) > 0:
        for args,nargs in zip(args,range(len(args))):
            rows = itertools.product(c.values,args)
            a = pd.DataFrame.from_records(rows, columns=['C_0','C_'+str(nargs+1)])
            b = pd.DataFrame(a['C_0'].tolist(), columns=c.columns.tolist())
            c = pd.concat([b,a.iloc[:,-1]], axis=1)
            #return c
    #else:
    #    b = pd.DataFrame(a['A_0'].tolist())
    #    c = pd.concat([b,a.iloc[:,-1]], axis=1)
            #res = pd.DataFrame(columns=name)
    ka=0
    for nc in range(c.shape[1]):
        if type(c.iloc[:,nc].values[0]) == np.ndarray:
            aux = pd.DataFrame(c[c.columns[nc]].tolist())
            aux_df = pd.concat([c.iloc[:,:nc], aux, c.iloc[:,nc+1:]], axis=1)
            ka+=1
    #print('forDONE')
    if ka>0:
        aux_df.columns = name
    else:
        aux_df = c
        aux_df.columns = name
    return aux_df
#%%
def df_generator():

    """
        This function generates a multindex dataframe with all parameters observed in BLUES model
        It takes the parameters as columns and technologies, nodes and modes as indexes
    """
# =============================================================================
# Generating multindex dataframe
# =============================================================================
    
    my_index = pd.MultiIndex(levels=[['0'],['0'], ['0']],
                                 codes=[['0'],['0'],['0']],
                                 names=['tec',  'node', 'mode']) #defining tec, nod and mode as index    
    
    t_col1 = ['letter_id', 'type', #'node', 
              'node_in', 'node_out', 'years', 'time', 'vintages'] #first tire of columns
    
    #creating multindex dataframe for cols 1
    t_col_1 = pd.MultiIndex.from_product([t_col1, ['1']], names = ['parameter', 'num'])
    #, sort_levels=False) #THIS IS MANNUALY MODIFIED BASED ON https://github.com/pandas-dev/pandas/pull/14062/files#diff-16bad7e1686739db18ced89420aaa349
    
    dft = pd.DataFrame(columns=t_col_1, index=my_index)

    #second tire of columns
    t_col2 = [                    'capacity_factor', 'optm_t',
                                   'minimum_utilization_factor', 
                                   'construction_time',
                                   'lifetime',
                                   #'investment_cost_years', 
                                   'investment_cost', #'fixed_cost_years', 
                                   'fixed_cost', 
                                   #'variable_cost_mode', 'variable_cost_years', 
                                   'variable_cost', 
                                   'historical_new_capacity_years', 'historical_new_capacity', 
                                   #'growth_new_capacity_up_years', 'growth_new_capacity_up', #NO VALUES FOR THESE PARAMETERS
                                   #'bound_new_capacity_up_years', 
                                   'bound_new_capacity_lo', 'bound_new_capacity_up',
                                   'bound_total_capacity_lo', 'bound_total_capacity_up',#'bound_activity_lo_years', 
                                   'bound_activity_lo',
                                   #'bound_activity_up_years', 
                                   'bound_activity_up', 
                                   'emissions',  'emission_values', #'P_L_F',                                   
                                   'act_constraints', 'act_constraint_values', 
                                   'cap_constraints', 'cap_constraint_values',
                                    'first_year', 'last_year','land_constraint',
                                    'historical_land_constraint',
                                   'comments'
                                   ]
    
    #creating multindex dataframe for cols 2
    t_col_2 = pd.MultiIndex.from_product([t_col2, ['1']], names = ['parameter', 'num'])
    #,sort_levels=False) #THIS IS MANNUALY MODIFIED BASED ON https://github.com/pandas-dev/pandas/pull/14062/files#diff-16bad7e1686739db18ced89420aaa349
    
    dft2 = pd.DataFrame(columns=t_col_2, index=my_index)
    
    #adding columns multindex for input
    inpc = ['input_commodities', 'input_levels', 'input_values'] 
    inp_v = [str(a) for a in list(range(1,20))]
    
    #adding columns multindex for output    
    outc = ['output_commodities', 'output_levels', 'output_values']
    out_v = [str(a) for a in list(range(1,20))]
    
    inp_tpa = [] #aux list
    out_tpa = [] #aux list
    
    #getting the multindex of commodities, levels and values for input and output
    for pair in itertools.product(inp_v, inpc): #input (commodity,1; commoodity, 2 ....)
        inp_tpa.append(pair)
    
    inp_tp = [(t[-1], t[0]) for t in inp_tpa]   #adjusting order
     
    for pair in itertools.product(out_v, outc, ): #output (commodity,1; commoodity, 2 ....)
        out_tpa.append(pair)
        
    out_tp = [(t[-1], t[0]) for t in out_tpa]   #adjusting order
    
    id_col = inp_tp+out_tp #put them together
    
    #generating df from tuples    
    s = pd.MultiIndex.from_tuples(id_col, names = ['parameter', 'num'])
    df0 = pd.DataFrame(columns=s, index=my_index) #create multindex df
    
    #DOING THAT BC OTHERWISE IT MISS THE CORRECT ORDER
    df1 = dft.join(df0, how='outer', sort=False) #MERGING
    dfg = df1.join(dft2, how='outer', sort=False) #MERGING
    
    dfg = dfg.astype('object')
    return dfg
#%%
def conv_data(tpp, adb_files, ldb, vtgs,year_act, BLUES, sc2_e, df_generator):#adb_files,
    """
        This function gets the info from the message adb/ldb file
        and organize it in tables and dictionaries
        a dictionary with all technologies constraints variables of each ldb file (each region) is exported
        as well as a csv file with all information of all technologies for each region
        The function return a dictionary with all technologies and their parameters in a list for each region
        It is necessary to define the path in which your projects are stored and the name of the project
    """
    startTime = datetime.now() #starting time
    print('\nworking on conversion')
    task, path, project = tpp #getting data
    idx = pd.IndexSlice    #pandas slice
    
    #MAIN DF
    df_f = df_generator()
    
    #Getting all the Energy Levels and Energy forms into a dictionary to relate the name and the letters
    #aa = np.load(os.path.join(path, project,project+'_energy_forms_dict.npy'), allow_pickle=True).item() #loading from adb
    
    relations = ['relationsc:', 'relationsp:', 'relationss:', 'relations1:', 'relations2:','variables:'] #relations type           

    #E_form is aggregating all required lines --> lines between 'energyforms' and 'demand'
    Level_form_dict = {}
    for Reg in adb_files:#looping through regions    
        Level_form_dict[Reg] = {}
        if task == 'ADB':
            dt = os.path.join(path, project, Reg,"data",Reg+'.adb') #path - energy forms
        elif task == 'LDB':
            dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+'.ldb') #path - energy forms

        dt2 = open(dt)#open file
        
        E_form_block = "" #creating block
        found = False #set non starting value
        for line in dt2: #looping through lines
            if found: #if start
                if line.strip() == 'demand:':  #stop condition
                    break
                E_form_block += line.strip() + "\n" #add line
            else: #if found is false
                if line.strip() == 'energyforms:': #starting condition
                    found = True #set as found true

        dt2.close() #close file

        E_forms = E_form_block[:] #str type has no copy
        E_forms = E_forms.split('*') #splitting each energy form
        E_forms = [e for e in E_forms if e != '' if e != '\n'] #excluding lines with no info 
        for n in E_forms: #looping through energy forms
            f1 = ' '.join(n.split('\t')).split('\n') #splitting by tab, then joining then splitting by line
            #fa = ' '.join(n.split('\t')).split('\n')
            f1 = [t.strip() for t in f1 if t!=''] #excluding empty lines
            LVL = f1[0].split(' ') #getting the level
            for com in f1[1:]:  #looping thorugh commoditties
                com = com.strip() #removing empty characteres in begining and in the end of str
                if not com.startswith('#'): #if it is not comment
                    c = com.split(' ') #split by space
                    c = [c for c in c if c!=""]
                    if len(c) > 2:
                        COM = c[-2]
                    else:
                        COM = c[-1]
                    Level_form_dict[Reg][COM+'-'+LVL[-1]] = LVL[0]+'|'+c[0].lower() #building dictionary level|commodity            
    
    main_rel_blk_final = []
    for Reg in adb_files:#looping through regions    

        if task == 'ADB':
            dt = os.path.join(path, project, Reg,"data",Reg+'.adb') #path - energy forms
        elif task == 'LDB':
            dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+'.ldb') #path - energy forms

        df = df_generator() 
        
#GETTING RELATIONS
        main_rel_blk = []
        for r in range(len(relations[:-1])): #looping through type of relations

            dt1 = open(dt) #open file
            rel_block = "" #generating block
            found = False #set atart
            for line in dt1: #looping through lines
                if found: #is start
                    if line.strip() == relations[r+1]: #stop if it reaches the next relation type
                        break
                    rel_block += line.strip() + "\n" #add line
                else: #ow found is false
                    if line.strip() == relations[r]: #check start condition
                        found = True    
            dt1.close() #close file
            
            rel_blk1 = rel_block[:] #str has noattribute copy
            rel_blk = rel_blk1.split('*\n') #spliting each technology
            main_rel_blk.append(rel_blk)
        
        #ADJUSTING TO GET RELATIONS
        main_rel_blk = [mm for mm in main_rel_blk if mm != [""]]
        main_rel_blk = list(itertools.chain(*main_rel_blk))
        main_rel_blk = [bb.split(' ')[0].split('\t')[0] for bb in main_rel_blk if bb.split(' ')[0].split('\t')[0]  != '']
        
        main_rel_blk_final.append(main_rel_blk)
        
    main_rel_blk_final = sorted(list(set(list(itertools.chain(*main_rel_blk_final)))))
        
    if BLUES:
        lc = [item for item in main_rel_blk_final if item.startswith('B') and item.endswith('a') or item.endswith('T')] #GET LAND CONSTRAINTS FOR 2010
        lsc2 = [item for item in main_rel_blk_final if item.startswith('B') and item.endswith('T')] #GET ALL LAND CONSTRAINTS
        lsc2a = [item for item in main_rel_blk_final if item.startswith('B') and item.endswith('a')] #GET LAND CONSTRAINTS for 2010 HISC
        sc2 = [item for item in main_rel_blk_final if not item.startswith('B') and item not in ['BASE', 'LOWC'] ] #removing cons related to land use and BASE and LOWC
        
        sc2 = [i for i in sc2 if not any(x in i for x in ['GHG-Brasil', 'CO2-Brasil', 'tGHG', 'TCO2', 'cGHG']) ] #removing emissions constraints defined
        sc2_c = [c for c in sc2 if c not in sc2_e] #get constraint that are not emissions
        #add constraints related to cement industry
# =============================================================================
#         sc2_c = sc2_c + ['BQkL',	'BQkU',	'BQBL',	'BQBU',	'BQCL',	'BQCU',	'BQKL',	'BQKU',	'BQHL',	'BQHU',	
#                          'BMQ2',	'BMQ3',	'BMQ4',	'BMQ6',	'BMZ1',	'BZ1U',	'BZ2U',	'BZ3U',	'BZ4U']
# =============================================================================
        #sc2_c = [c.split('-')[0] if '-' in c else c for c in sc2_c] #EXCLUDING RELATIONS RELATED TO BRAZIL
       

        sc2 = [i for i in sc2 if not any(x in i for x in ['GHG-Brasil', 'CO2-Brasil', 'tGHG', 'TCO2', 'cGHG']) ] + ['BCCS'] #removing emissions constraints defined
        sc2_c = [c for c in sc2 if c not in sc2_e] #get constraint that are not emissions
        #add constraints related to cement industry
        obcon = [item for item in main_rel_blk_final if item.startswith('B')]
        lcaux = list(set([ll[:3] for ll in lc]))
        obcon2 = [item for item in obcon if not item.startswith(tuple(lcaux))] + ['BCCS']
        sc2_c = sc2_c + obcon2
        sc2_c = [item for item in sc2_c if item not in ['BASE', 'LOWC'] ] #re
        sc2_c = sorted( list(set(sc2_c)) ) #make it cleaner
        
    else:
        sc2_c = main_rel_blk[::]
        
    if task == 'LDB': #loading all relations from adb
        ldb_aux = "_"+ldb
    else:
        ldb_aux = ldb

    #SAVING LIST OF RELATIONS
    with open(os.path.join(path, project, 'IX', ldb, project+'_'+task+'_relation_list'+ldb_aux), 'wb') as fp: #save list
        pickle.dump(sc2_c, fp)
    #EMISSIONS
    with open(os.path.join(path, project, 'IX',ldb, project+'_'+task+'_emission_list'+ldb_aux), 'wb') as fp:#save list
        pickle.dump(sc2_e, fp)
    if BLUES:
        #LAND RELATIONS
        with open(os.path.join(path, project, 'IX',ldb, project+'_'+task+'_land_relation_list'+ldb_aux), 'wb') as fp: #save list
            pickle.dump(lc, fp)          
 
    for Reg in adb_files:#looping through regions    

        if task == 'ADB':
            dt = os.path.join(path, project, Reg,"data",Reg+'.adb') #path - energy forms
        elif task == 'LDB':
            dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+'.ldb') #path - energy forms

        df = df_generator() 
                
        data_file2 = open(dt) #open adb
         
        l = open(dt)  #open adb
        l1 = l.read() #reading it
        l2 = l1.split('\n') #split by line
    #     
        node = [t for t in l2[0:10] if t.upper().startswith('ADB:')] #get nodes
        node = node[0].split(':')[-1].strip(' ') #adjusting format
         
        l.close() #close file
    # =============================================================================
    
    #Get the text between specific lines defined by its content
        tec_block = "" #creating technology block
        found = False #
        for line in data_file2: #looping through adb file
            if found: #if found is True add line
                if line.strip().startswith('endata') or line.strip().startswith('resources:'):  #if line is as defined stop
                    break
                tec_block += line.strip()+'\n' #adding line                
            else: #found is false
                if line.strip() == 'systems:': #if line is systems:
                    found = True #set found
                    #tec_block = "" #creating technology block
        
        data_file2.close()   #close file      
        
        tec_block_1 = tec_block[:] #copying block
        tec1a = tec_block_1.split('*\n') #spliting each technology
        tec1 = [e for e in tec1a if e != '' if e != '\n'] #excluding lines with no info#get rid of last element 
        
        #if BLUES:
            #tec_L =  [s for s in tec1 if s.startswith('Land') or s.startswith('Conv')]#eliminate land and conv after 10
            #tec_L1 = [s for s in tec_L if hasNumbers(s.split(' ')[0].split('\t')[0][-2:])] #getting technologies with number in the last two characters
            #tec_L2 = [s for s in tec_L1 if int(s.split(' ')[0].split('\t')[0][-2:]) >10] #keep only the values higher than 10
                #tec_L10 = [s for s in tec_L1 if int(s.split(' ')[0][-2:]) == 10] #keep only the values higher than 10                
            #tec = [t for t in tec1 if t not in tec_L2] #removing values higher than 10                
                #t_n = [s for s in tec1 if s not in tec] #list of tecs removed
        #else:
        tec = tec1[:]  

        #n = tec[tec.index([t for t in tec if t.startswith('LandUse_emissions_2010')][0])]
        kt = 0 #counting
        for n in tec:    
            #print(n)
            kt += 1 #counting
            f1 = ' '.join(n.split('\t')).split('\n') #adjusting based on tab and lines
            f1 = [t for t in f1 if t!=''] #removing empty lines
            f1 = [t.strip(' ') for t in f1] #removing first character if it is empty
            F0 = f1[0].split(' ') #Getting the name of the technology
            
            tn0 = '_'.join([fn[0].upper()+fn[1:].lower() for fn in F0[0].split("_")]) #[0].upper() + F0[0][1:] #name of the technology
            if BLUES:
                if tn0.startswith('Land_Bal'):
                    tn0 = tn0+'_10'  #auxliary list to adjust Land_Bal tech add '_10' into 'Land_Bal' tec  
                if tn0.startswith('Land') and tn0.endswith('_10') or tn0.startswith('Conv') and tn0.endswith('_10'):
                    tn0 = '_'.join(tn0.split('_')[:-1])  #removing 10 at the end of the name 
                #R_tec[Reg].append(tn0) #getting tecs by region
             
            tn0 = tn0.strip()
            mk = 1 #mode counting
            if any('activity' in s for s in f1): #checking if there is more than 1 mode
                for m in f1: #looping throuhg technology info
                   if 'activity' in m: #getting the activity line
                       mk+=1 #add mode 2 according to activities
                       ind = [i for i, j in enumerate(f1) if 'activity' in j] # getting the position of each activity
                       startk = True #starting element
                       indices = [] #aux list for indices
                       for i in ind: #looping through activity position in the list
                       #getting the interval of each activity in the list by indices
                           if startk: #start
                               indices.append((0,i)) #index from 0 to i which is the first activity
                               startk = False # add one to avoid kk 
                               auxi = i #getting the new first index
                           else: #going one
                               indices.append((auxi,i)) #index from the last one to the next activity stop
                               auxi = i #getting the new first index
                       indices.append((auxi,len(f1)))
                       nf1 = [f1[s:e] for s,e in indices] #splitting the list into lists based on activities indices
            else:#there is no activity
               nf1 = [f1] #if there is no activity get f1 
            
            #GENERAL
            nft = len(nf1) #number of modes           
####################################################################################################################################################################
            #DO
            kname = 0
            while (tn0, node) in [(t0,n0) for t0,n0 in zip(df.index.get_level_values(level=0).tolist(), df.index.get_level_values(level=1).tolist() )]:
                kname+=1
                tn0 = tn0 +'_'+str(kname)
            tp = list(zip(cycle([tn0]), cycle([node]), list(range(1,nft+1)))) #combining the name of the tech and the mode (activity)
            ss = pd.DataFrame(columns=df.columns, index=pd.MultiIndex.from_tuples(tp, names=['tec', 'node','mode'])) #creating series
            df = df.append(ss) #adding series as row
            df = df.sort_index() #by sorting dataframe we improve performance apparently
            print(tn0+' - '+Reg) #print tec name and its region                
            if any(y > 1 for y in [len(ntnt) for ntnt in nf1]): #check if there is ANY activity with at least one parameter (it is more than one because in the list there are always the name of the tec or 'activity')

# =============================================================================
#                           ADDING NEW LIFETIME
# =============================================================================              
                if any(i.startswith('pll') for i in nf1[0]):
                    LFT = [i for i in nf1[0] if i.startswith('pll')] #get pll str
                    LFT = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", LFT[0].split('#')[0]) #get values from pll str
                    LFT = [r_t_z("{:.11f}".format(float(b))) for b in LFT]  #converting to float

                    df.loc[idx[tn0, node, :], idx['lifetime', '1']] = [str(LFT)]*nft    #add into df     

# =============================================================================
#                           ADDING NEW FISRT YEAR
# =============================================================================            
                if any(i.startswith('fyear') for i in nf1[0]):  #if there is any fyear
                    #FYR = [[float(i.split(' ')[-1].split('#')[0]) if i.split(' ')[-1].split('#')[0] != '' else float(i.split(' ')[-2].split('#')[0]) for i in nf1[0] if i.startswith('fyear')]]
                    FYR = [i for i in nf1[0] if i.startswith('fyear')] #get pll str     
                    FYR = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", FYR[0].split('#')[0]) #get values from pll str
                    FYR = [float(b) for b in FYR]  #converting to float     
                    
                    df.loc[idx[tn0, node, :], idx['first_year', '1']] = [str(FYR)]*nft   #add as it is into df
                         
# =============================================================================
#                           ADDING NEW LAST YEAR
# =============================================================================  
                if any(i.startswith('lyear') for i in nf1[0]): #if there is any fyear
                    #FLY = [[float(i.split(' ')[-1].split('#')[0]) if i.split(' ')[-1].split('#')[0] != '' else float(i.split(' ')[-2].split('#')[0]) for i in nf1[0] if i.startswith('lyear')]]                                                     
                    FLY = [i for i in nf1[0] if i.startswith('lyear')] #get pll str     
                    FLY = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", FLY[0].split('#')[0]) #get values from pll str
                    FLY = [float(b) for b in FLY]  #converting to float     
                    
                    df.loc[idx[tn0, node, :], idx['last_year', '1']] = [str(FLY)] *nft   #add as it is into df
                
# =============================================================================
#                           ADDING NEW CAPACITY FACTOR --> PFL CAPACITY FACTOR
# =============================================================================           
                if any(i.startswith('plf') for i in nf1[0]): #checking for plf
                    PLF = [i for i in nf1[0] if i.startswith('plf')] #checking for plf
                    PLF = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", PLF[0].split('#')[0]) #getting all numbers
                    PLF = [r_t_z("{:.11f}".format(float(b))) for b in PLF] #converting str into value    
                
                    df.loc[idx[tn0, node, :], idx['capacity_factor', '1']] = [str(PLF)] * nft #adding as list

# =============================================================================
#                           ADDING NEW INVESTMENT COST
# ============================================================================= 
                if any(i.startswith('inv') for i in nf1[0]): #checking for inv
                    INV = [i for i in nf1[0] if i.startswith('inv')]#checking for inv
                    INV = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", INV[0].split('#')[0])  #getting all numbers
                    INV = [r_t_z("{:.11f}".format(float(b))) for b in INV]  #converting str into value                    
                
                    df.loc[idx[tn0, node, :], idx['investment_cost', '1']] = [str(INV)] * nft #adding as list
   
# =============================================================================
#                           ADDING NEW FIXED COST - FOM
# =============================================================================   
                if any(i.startswith('fom') for i in nf1[0]): #checking for fom
                    FOM = [i for i in nf1[0] if i.startswith('fom')] #checking for fom
                    FOM = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", FOM[0].split('#')[0]) #getting all numbers
                    FOM = [r_t_z("{:.11f}".format(float(b))) for b in FOM] #converting str into value
                
                    df.loc[idx[tn0, node, :], idx['fixed_cost', '1']] = [str(FOM)] * nft #adding as list

# =============================================================================
#                           ADDING NEW CONSTRUCTION TIME - CTIME
# ============================================================================= 
                if any(i.startswith('ctime') for i in nf1[0]): #checking for ctime
                    CTIME = [i for i in nf1[0] if i.startswith('ctime')] #checking for ctime
                    CTIME = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", CTIME[0].split('#')[0])  #getting values
                    CTIME = [r_t_z("{:.11f}".format(float(b))) for b in CTIME] #converting str into value
                
                    df.loc[idx[tn0, node, :], idx['construction_time', '1']] = [str(CTIME)] * nft #adding as list

# =============================================================================
#                           ADDING NEW OPTM - #optm_t as capacity factor info
# =============================================================================                        
                if any(i.startswith('optm') for i in nf1[0]):  #checking for OPTM
                    OPTM = [i for i in nf1[0] if i.startswith('optm')] #checking for OPTM
                    OPTM = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", OPTM[0].split('#')[0]) #getting all numbers
                    OPTM = [r_t_z("{:.11f}".format(float(b))) for b in OPTM] #converting str into value   

                    df.loc[idx[tn0, node, :], idx['optm_t', '1']] = [str(OPTM)] * nft#adding as list

# ============================================================================= 
#                           ADDING NEW MIN UTILIZATION FACTOR - minutil
# =============================================================================                      
                if any(i.startswith('minutil') for i in nf1[0]): #checking for muf  
                    MIN_UTI = [i for i in nf1[0] if i.startswith('minutil')] #checking for muf  
                    MIN_UTI = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", MIN_UTI[0].split('#')[0]) #getting all numbers
                    MIN_UTI = [r_t_z("{:.11f}".format(float(b))) for b in MIN_UTI]       #converting str into value       
                
                    df.loc[idx[tn0, node, :], idx['minimum_utilization_factor', '1']] = [str(MIN_UTI)] * nft #adding as list
               
# =============================================================================
#                           ADDING TOTAL BOUND INSTALLED BDI UP
# ============================================================================= 
                if any(i.startswith('bdi up') for i in nf1[0]): #checking for bdi up  
                    BDI_up = [i for i in nf1[0] if 'bdi up' in i] #checking for bdi up  
                    BDI_up = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", BDI_up[0].split('#')[0]) #getting all float digits
                    BDI_up = [r_t_z("{:.11f}".format(float(b)) ) for b in BDI_up] #converting str into value
                    
                    df.loc[idx[tn0, node, :], idx['bound_total_capacity_up', '1']] = [str(BDI_up)] * nft #adding as list   
 
# =============================================================================
#                           ADDING TOTAL BOUND INSTALLED BDI LO
# =============================================================================                 
                if any(i.startswith('bdi lo') for i in nf1[0]): #checking for bdi lo  
                    BDI_lo = [i for i in nf1[0] if 'bdi lo' in i] #checking for bdi lo  
                    BDI_lo = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", BDI_lo[0].split('#')[0]) #getting all float digits
                    BDI_lo = [r_t_z("{:.11f}".format(float(b)) ) for b in BDI_lo] #converting str into value      
                    
                    df.loc[idx[tn0, node, :], idx['bound_total_capacity_lo', '1']] = [str(BDI_lo)] * nft #adding as list                    
         
# =============================================================================
#                           ADDING NEW BOUND NEW CAPACITY BDC UP
# =============================================================================                 
                if any(i.startswith('bdc up') for i in nf1[0]): #checking for bdc up
                    BDC_up = [i for i in nf1[0] if 'bdc up' in i] #checking for bdc up
                    BDC_up = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", BDC_up[0].split('#')[0]) #getting all float digits
                    BDC_up = [r_t_z("{:.11f}".format(float(b)) ) for b in BDC_up] #converting str into value   
                    BDC_up =  [0] + BDC_up #adjusting bdc in the first because IX and ORIGINAL MSG have different time step approach [str(float(BDC_up[0])/2)]
                    #FOR INSTANCE IX 2015 REPRESENTS 2011-2015; ORIGINAL REPRESENT 2015-2019
                    if len(BDC_up) > len(year_act):
                        BDC_up = BDC_up[:len(year_act)]
                        
                    df.loc[idx[tn0, node, :], idx['bound_new_capacity_up', '1']] = [str(BDC_up)] * nft #adding as list
             
# =============================================================================
#                           ADDING NEW BOUND NEW CAPACITY BDC LO
# =============================================================================                 
                if any(i.startswith('bdc lo') for i in nf1[0]): #checking for bdc lo
                    BDC_lo = [i for i in nf1[0] if 'bdc lo' in i] #checking for bdc up
                    BDC_lo = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", BDC_lo[0].split('#')[0]) #getting all float digits
                    BDC_lo = [r_t_z("{:.11f}".format(float(b)) ) for b in BDC_lo] #converting str into value     
                    BDC_lo = [0] + [r_t_z("{:.11f}".format(float(b)) ) for b in BDC_lo] #adjusting bdc in the first because IX and ORIGINAL MSG have different time step approach
                    #FOR INSTANCE IX 2015 REPRESENTS 2011-2015; ORIGINAL REPRESENT 2015-2019
                    if len(BDC_lo) > len(year_act):
                        BDC_lo = BDC_lo[:len(year_act)]    
                        
                    df.loc[idx[tn0, node, :], idx['bound_new_capacity_lo', '1']] = [str(BDC_lo)] * nft #adding as list
  
# =============================================================================
#                           ADDING NEW HISTORICAL NEW CAPACITY
# =============================================================================                      
                if any(i.startswith('hisc') for i in nf1[0]): #check if there is any hisc
                    HIST = [i for i in nf1[0] if 'hisc' in i] #check if there is any hisc
                    HIST = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", HIST[0].split('hc')[-1].split('#')[0]) #getting all float digits #r"[-+]?d*.d+|d+"
  
                    HIST_y = [int(float(HIST[y])) for y in range(len(HIST)) if y%2==0] #get year and converting them to int
                    HIST_val = [r_t_z("{:.11f}".format(float(HIST[y])) ) for y in range(len(HIST)) if y%2!=0]    #get values and converting them to float
                
                    sums = {} #aux sum dict
                    for key, value in zip(HIST_y,HIST_val): #looping thorugh year and values hisc pair
                        try: #try sum if it already exists
                            sums[key] += float(value) #adding value
                        except KeyError: #ow
                            sums[key] = float(value) #defien value
                    
                    HIST_y = list(sums.keys()) #defining years from sum
                    HIST_val = list(sums.values()) #defining values from sum                    
                    HIST_val = [r_t_z("{:.11f}".format(v)) for v in HIST_val]
                    
                    df.loc[idx[tn0, node, :], idx['historical_new_capacity_years', '1']] = [str(HIST_y)] * nft #adding as list
                    df.loc[idx[tn0, node, :], idx['historical_new_capacity', '1']] = [str(HIST_val)] * nft #adding as list                    


########################################################################################################################                    
                for nne in range( nft ): #looping through tec modes 
                    ff = nf1[nne] #ff is the current mode list of the current tec
                    if len(ff) >1: #check if there is at least one parameter
                        #DO
                        mod2 = nne+1 #mode
                        tn = (tn0, node, mod2)
                        comments = [c for c in ff if c.startswith('#')] #getting comments
                        c_aux = [(c[:3], c.partition('#')[-1]) for c in ff if c.partition('#')[-1] != ''] #getting all lines that contain comments
                        if len(c_aux) >1: #if there are comments
                            c_aux = c_aux[:-1] #remove comment that is already there
                            comments.append(c_aux) #appending other comments into main comment 
                        ff = [f.partition('#')[0] for f in ff] # excluding comments related to the parameters to avoid future issues
                        ff = [ft for ft in ff if ft != ''] #removing empty elements
                        ff = [ft.strip(' ') for ft in ff] # excluding blank space at the end of the line
                    
# =============================================================================
#                           ADDING NEW COMMENTS
# =============================================================================                            
                        if any(len(y) > 0 for y in comments):

                            df.loc[tn, idx['comments', '1']] = str(comments) #add into ldb df #at. can replace loc.                    
# =============================================================================
#                           ADDING NEW INPUTS
# =============================================================================               
                        #INPUT
                        inp = [i for i in ff if i.startswith('inp') or i.startswith('minp')] #getting the inputs
                        if len(inp) > 0: #if there is any input
                            if any(['minp' in s for s in inp]):
                                iv = 0 #counting inputs
                            else:
                                iv = 1
                            for i in inp:#looping through inputs
                                iv+=1 #counting inputs
                                inp1 = [aa for aa in i.split(' ') if aa!='']#splitting the input, the code and te values
                                if i.startswith('minp'): #check if there is a main input
                                    if len(inp1[1].split('-'))>2: #if >2 the node of the input is different from the current node of the tech
                                        df.loc[tn, idx['node_in', '1']] = inp1[1].split('-')[-1] #setting node in  
                                        lcd = Level_form_dict[inp1[1].split('-')[-1] ]
                                    else: #ow same node
                                        df.loc[tn, idx['node_in', '1']] = node #otw the node is the same
                                        lcd = Level_form_dict[ node ]
                                else:
                                    df.loc[tn, idx['node_in', '1']] = node #otw the node is the same
                                    lcd = Level_form_dict[ node ]
                                            
                                ii = ' '.join(str(s) for s in inp1[2:]) #get all inp info into a string separeted by " "
                                
                                if 'E' in ii.split('#')[0].upper():
                                   print(tn0+' - '+Reg + ' HAS E IN INPUT') #print tec name and its region   
                                   inp_val = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?',ii.split('#')[0])
                                else:
                                    inp_val = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", ii.split('#')[0]) #getting all numbers of the input #r"bd[d,.]*b"
                                                         
                                inp_val = [r_t_z("{:.11f}".format(float(b)) ) for b in inp_val] # converting str into value
                                if len(inp_val) > len(vtgs): #check if len inp val is greater than vtgs
                                    inp_val = inp_val[0:len(vtgs)] #adjust it
                               
                                df.loc[tn, idx['input_values', str(iv)]] = str(inp_val) #adding value as list
                                df.loc[tn, idx['input_commodities', str(iv)]] = lcd[inp1[1][:3]].split('|')[-1] #adding commodity
                                df.loc[tn, idx['input_levels', str(iv)]] = lcd[inp1[1][:3]].split('|')[0] #adding level                                   
# =============================================================================
#                           ADDING NEW OUTPUTS
# =============================================================================                     
                        #OUTPUT 
                        out = [i for i in ff if i.startswith('out') or i.startswith('mout')] #getting the inputs
                        if len(out) > 0: #if there is any output
                            if any(['moutp' in s for s in out]): #check if it has main output
                                iv = 0 #counting outputs
                            else:
                                iv = 1               
                            for o in out: #looping through inputs
                                iv+=1 #counting outputs                                
                                out1 = [aa for aa in o.split(' ') if aa!=''] #splitting the output, the code and the values
                                if o.startswith('mout'): #check if there is a main input
                                    if len(out1[1].split('-'))>2: #if >2 the node of the output is different from the current node of the tech
                                        df.loc[tn, idx['node_out', '1']] = out1[1].split('-')[-1] #setting node in
                                        lcd = Level_form_dict[out1[1].split('-')[-1] ]
                                    else:
                                        df.loc[tn, idx['node_out', '1']] = node #otw the node is the same
                                        lcd = Level_form_dict[node ]                                       
                                else:
                                    df.loc[tn, idx['node_out', '1']] = node #otw the node is the same
                                    lcd = Level_form_dict[ node ]
                                            
                                oo = ' '.join(str(s) for s in out1[2:]) #get all inp info into a string separeted by " "
                                
                                if 'E' in oo.split('#')[0].upper():
                                   print(tn0+' - '+Reg + ' HAS E IN OUTPUT') #print tec name and its region   
                                   out_val = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?',oo.split('#')[0])

                                    
                                else: 
                                    out_val = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", oo.split('#')[0]) #getting all numbers of the output
                                                         
                                out_val = [r_t_z("{:.11f}".format(float(b)) ) for b in out_val] # converting str into value
                                
                                if len(out_val) > len(vtgs): #check if len inp val is greater than vtgs
                                    out_val = out_val[0:len(vtgs)] #adjust it
                                    
                                df.loc[tn, idx['output_values', str(iv)]] = str(out_val) #adding value as list
                                df.loc[tn, idx['output_commodities', str(iv)]] = lcd[out1[1][:3]].split('|')[-1] #adding commodity
                                df.loc[tn, idx['output_levels', str(iv)]] = lcd[out1[1][:3]].split('|')[0] #adding level
# =============================================================================
#                           ADDING NEW EMISSIONS
# =============================================================================  
                        emi_fac = [] #emission factor aux list
                        emi_val = [] #emission value aux list
                        emi = [i for i in ff if i.startswith('con1a') or i.startswith('conca')]# ##constraints info or i.startswith('conc')]
                        emi = [i for i in emi if any(substring in i for substring in sc2_e)] 
                        for ee in emi: #looping through constraints
                            emi1 = [aa for aa in ee.split(' ') if aa!=''] #spliting by blank space
                            emi_fac_aux = [e.strip() for e in emi1 if e.split('-')[0] in sc2_e] #check if values are among the emissions previously defined  
                            if len(emi_fac_aux) > 0:  #check if there is any value
                                emi_fac.append(emi_fac_aux[0]) #appending value into main list of emission factors
                                
                                if 'E' in ee.split(emi_fac_aux[0])[-1].split('#')[0].upper():
                                    print(tn0+' - '+Reg + ' HAS E IN EMISSION') #print tec name and its region   
                                    ev_aux = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', ee.split(emi_fac_aux[0])[-1].split('#')[0] )
                                else:
                                    ev_aux = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", ee.split(emi_fac_aux[0])[-1].split('#')[0] ) #getting all numbers of the output emi1[ind+1:] #from ts to the end
                                
                                ev_aux = [r_t_z("{:.11f}".format(float(v)) ) for v in ev_aux] #converting to float   
                                emi_val.append(ev_aux) #appending value into main list of emission values  
                        if len(emi_fac)>0: #if there is any factor                                                        
                            df.loc[tn, idx['emissions', '1']] = str(emi_fac).replace("'","") #add into df
                            df.loc[tn, idx['emission_values', '1']] = str(emi_val) #add into df
                    
# =============================================================================
#                           ADDING NEW ACT CONSTRAINTS
# =============================================================================  
                        con_fac = [] #CON factor aux list
                        con_val = [] #CON value aux list
                        con = [i for i in ff if i.startswith('con1a') or i.startswith('conca')]# or i.startswith('conc')] #constraints info
                        con = [i for i in con if any(substring in i for substring in sc2_c)]
                        for ccon in con: #looping through constraints
                            con1 = [aa for aa in ccon.split(' ') if aa!=''] #spliting by blank space
                            nc_aux = [c.strip() for c in con1 if c.split('-')[0] in sc2_c]
                            if len(nc_aux) > 0: #check if there is any value
                                con_fac.append(nc_aux[0]) #appending value into main list of constraints
                                
                                if  'E' in ccon.split(nc_aux[0])[-1].split('#')[0].upper():
                                    print(tn0+' - '+Reg + ' HAS E IN CON1A') #print tec name and its region   
                                    cv_aux = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', ccon.split(nc_aux[0])[-1].split('#')[0] )
                                else:
                                    cv_aux = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", ccon.split(nc_aux[0])[-1].split('#')[0]) #getting all numbers of the output emi1[ind+1:] #from ts to the end
                                
                                cv_aux = [r_t_z("{:.11f}".format(float(v)) ) for v in cv_aux] #converting to float                                                                           
                                con_val.append(cv_aux) #appending value into main list of constraint values 
                        if len(con_fac)>0:                           
                            df.loc[tn, idx['act_constraints', '1']] = str(con_fac).replace("'","") #add into df
                            df.loc[tn, idx['act_constraint_values', '1']] = str(con_val)    #add into df                       

# =============================================================================
#                           ADDING NEW CAP CONSTRAINTS
# =============================================================================  
                        con_fac = [] #CON factor aux list
                        con_val = [] #CON value aux list
                        con = [i for i in ff if i.startswith('con1c') ]# or i.startswith('conc')] #constraints info
                        con = [i for i in con if any(substring in i for substring in sc2_c)]
                        for ccon in con: #looping through constraints
                            con1 = [aa for aa in ccon.split(' ') if aa!=''] #spliting by blank space
                            nc_aux = [c.strip() for c in con1 if c.split('-')[0] in sc2_c]
                            if len(nc_aux) > 0: #check if there is any value
                                con_fac.append(nc_aux[0]) #appending value into main list of constraints
                                cv_aux = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", ccon.split(nc_aux[0])[-1].split('#')[0]) #getting all numbers of the output emi1[ind+1:] #from ts to the end
                                cv_aux = [r_t_z("{:.11f}".format(float(v)) ) for v in cv_aux] #converting to float                                                                           
                                con_val.append(cv_aux) #appending value into main list of constraint values 
                        if len(con_fac)>0:                           
                            df.loc[tn, idx['cap_constraints', '1']] = str(con_fac).replace("'","") #add into df
                            df.loc[tn, idx['cap_constraint_values', '1']] = str(con_val)    #add into df                       

# =============================================================================
#                           ADDING NEW LAND CONSTRAINTS
# =============================================================================  
                        if BLUES:
                            ll_fac = [] #land con aux list
                            lcon = [i for i in ff if i.startswith('con1a') or i.startswith('con1c')]# or i.startswith('conc')] #constraints info
                            lcon = [i for i in lcon if any(substring in i for substring in lsc2)] #check if there is a total land constraint
                            for llc in lcon: #looping through constraints
                                lcon1 = [aa for aa in llc.split(' ') if aa!=''] #spliting by blank space
                                lc_aux = [l.strip() for l in lcon1 if l in lsc2] #get only land constraints   
                                if len(lc_aux) > 0: #check if there is any value
                                    ll_fac.append(lc_aux[0]) #appending value into main list of land constraints
                            if len(ll_fac)>0:                          
                                df.loc[tn, idx['land_constraint', '1']] = str(ll_fac).replace("'","") #add into df

# =============================================================================
#                           ADDING HISTORICAL LAND CONSTRAINTS
# =============================================================================  
                            if tn0.startswith('Land'):
                                ll_fac2 = [] #land con aux list
                                lcon2 = [i for i in ff if i.startswith('con1a') or i.startswith('con1c')]# or i.startswith('conc')] #constraints info
                                lcon2 = [i for i in lcon2 if any(substring in i for substring in lsc2a)] #check if constraint is in lsc2a constraints that ends with 'a'
                                for llc2 in lcon2: #looping through constraints
                                    lcon3 = [aa for aa in llc2.split(' ') if aa!=''] #spliting by blank space
                                    hlcf_aux = [hl for hl in lcon3 if hl in lsc2a] #get only land constraints   
                                    if len(hlcf_aux) > 0: #check if there is any value
                                        ll_fac2.append(hlcf_aux[0]) #appending value into main list of land constraints
                                if len(ll_fac2)>0:                          
                                    df.loc[tn, idx['historical_land_constraint', '1']] = str(ll_fac2).replace("'","") #add into df
  
                                 
# =============================================================================
#                           ADDING NEW VARIABLE COST - VOM
# ============================================================================= 
                        #VARIABLE COST - VOM #VOM varies according to the mode/activity
                        if any(i.startswith('vom') for i in ff): #checking for vom
                            VOM = [i for i in ff if i.startswith('vom')] #checking for vom
                            VOM = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", VOM[0].split('#')[0]) #getting all numbers
                            VOM = [r_t_z("{:.11f}".format(float(b)) ) for b in VOM]  #converting str into value
                            df.loc[tn, idx['variable_cost', '1']] = str(VOM)  #adding as list


# =============================================================================
#                           ADDING BOUND ACTIVITY BDA UP
# =============================================================================                 
                        #BOUND ACTIVITY BDA
                        if any(i.startswith('bda up') for i in ff): #checking for bda up  
                            BDA_up = [i for i in ff if 'bda up' in i] #checking for bda up  
                            BDA_up = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", BDA_up[0].split('#')[0]) #getting all float digits
                            BDA_up_val = [r_t_z("{:.11f}".format(float(b)) ) for b in BDA_up] #converting str into value
                            df.loc[tn, idx['bound_activity_up', '1']] = str(BDA_up_val) #adding as list

# =============================================================================
#                           ADDING BOUND ACTIVITY BDA LO
# =============================================================================                
                        #BDA LO
                        if any(i.startswith('bda lo') for i in ff): #checking for bda lo
                            BDA_lo = [i for i in ff if 'bda lo' in i] #checking for bda lo
                            BDA_lo = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", BDA_lo[0].split('#')[0]) #getting all float digits
                            BDA_lo_val = [r_t_z("{:.11f}".format(float(b)) ) for b in BDA_lo] #converting str into value
                            df.loc[tn, idx['bound_activity_lo', '1']] = str(BDA_lo_val) #adding as list
                          
# =============================================================================
#                           ADDING BOUND ACTIVITY fixed
# =============================================================================                 
                        #BOUND ACTIVITY BDA
                        if any(i.startswith('bda fx') for i in ff): #checking for bda fx
                            BDA = [i for i in ff if 'bda fx' in i] #checking for bda fx
                            BDA = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", BDA[0].split('#')[0]) #getting all float digits
                            BDA_val_up = [r_t_z("{:.11f}".format(float(b)+0) ) for b in BDA] #converting str into value
                            BDA_val_lo = [r_t_z("{:.11f}".format(float(b)-0) ) for b in BDA] #converting str into value
                            
                            df.loc[tn, idx['bound_activity_up', '1']] = str(BDA_val_up) #adding as list
                            df.loc[tn, idx['bound_activity_lo', '1']] = str(BDA_val_lo) #adding as list

################################################################################################################################################ 
        df.drop(('0', '0', '0'), inplace=True)        
        df_f = pd.concat([df_f, df])#, join='inner')

    df_f.drop(('0', '0', '0'), inplace=True)     
    dfc = df_generator()
    df_f = df_f[dfc.columns]
    dt_final =  datetime.now() - startTime #finish time
    print('Done.\n')
    print('Runtime:    '+ str(dt_final)+'\n')

    
    if task == 'LDB':
        df_f.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_'+ldb+'_original_df'+'.csv'), sep=';') #exporting ldb_df
        np.save(os.path.join(path, project,'IX',ldb,project+'_'+task+'_energy_forms_dict_'+ldb+'.npy'),Level_form_dict) #save level|commodity dictionary
    else:
        df_f.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_original_'+'df'+'.csv'), sep=';') #exporting adb_df
        np.save(os.path.join(path, project,'IX',ldb,project+'_'+task+'_energy_forms_dict.npy'),Level_form_dict) #save level|commodity dictionary   
    return df_f        

#%%
def resource_conv(tpp, adb_files, ldb=''):
    """
    get resource commodities, values and grades from adb/ldb files
    if it is ldb it replaces old values for the new ones
    """


    print('\nworking on resources')
    task, path, project = tpp #getting data
    #os.chdir(path)   # set working drive 
    if ldb == "":
        resource_df = pd.DataFrame(columns=['node', 'commodity', 'grade', 'value']) #aux resource df
        dd1 = np.load(os.path.join(path,project,'IX', ldb,project+'_'+task+'_energy_forms_dict.npy'),allow_pickle='TRUE').item()  #load egergy form dictionary        
        for Reg in adb_files: #looping through regions
            dd = dd1[Reg] 

            dt = os.path.join(path, project, Reg,"data", Reg+".adb")  #adb path
            data_file2 = open(dt) #open file       
            l = open(dt) #open file
            l1 = l.read() #read it
            l2 = l1.split('\n') #spliting by line
            
            node = [t for t in l2[0:10] if t.upper().startswith('ADB:')] #get node
            node = node[0].split(':')[-1].strip(' ')  #adjust it
             
            l.close()#close file
        # =============================================================================
        
        #Get the text between specific lines defined by its content
            res_block = "" #creating resource block
            found = False     #set to not start
            for line in data_file2: #looping throguh lines
                if found: #if start 
                    if line.strip() == 'endata': #if attend this stop
                        break
                    res_block += line.strip() + "\n" #add line                    
                else: # foud is false
                    if line.strip() == 'resources:': #look for this condition
                        found = True #set as true
                        #tec_block = ""
            data_file2.close()    #close file     
            
            res1 = res_block[:] #copying it
            res = res1.split('*\n') #spliting each technology
            res = [d for d in res if d != ""]  #get free from empty elements
            
            for tt in res: #looping through resources
                tt1 =  tt.split('\n') #splitting by line
                #generating df
                aux = pd.DataFrame({'node': node, #node
                                    'commodity': dd[tt1[0].split(' ')[-1]].split('|')[-1], #commodity from energy form
                                    'grade': tt1[2].split(' ')[-1], #grade
                                    'value': int(tt1[3].split('\t')[-1]) #value
                                    }, index=[0])
                resource_df = pd.concat([resource_df, aux], ignore_index=True) #concat
        resource_df.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_resource_'+'df'+'.csv'), index=False, sep=';') #exporting to csv            

    elif ldb != '': #if ldb 
        #load egergy form dictionary 
        try: #load resource df from adb
            resource_df = pd.read_csv(os.path.join(path,project,'IX', project+'_ADB_resource_'+'df'+'.csv'), dtype='object', sep=';')
        except ValueError:
            raise Exception("There is no resource df in the current path:"+ os.path.join(path,project, project+'_resource_'+'df'+'.csv')+"\n\ntry to set ldb = '' or check your pathes") 
        dd1 = np.load(os.path.join(path, project,'IX',ldb, project+'_'+task+'_energy_forms_dict_'+ldb+'.npy'),allow_pickle='TRUE').item()  
        new_rr = []        
        for Reg in adb_files: #looping through regions
            dd = dd1[Reg]

            dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+".ldb") #adb path
            data_file2 = open(dt)            #open file
            l = open(dt)#open file
            l1 = l.read()   #read lines
            l2 = l1.split('\n') #split by line

            node = [t for t in l2[0:10] if t.upper().startswith('ADB:')] #node
            node = node[0].split(':')[-1].strip(' ') #node
             
            l.close() #close file
        # =============================================================================
        
        #Get the text between specific lines defined by its content
            res_block = "" #creating resource block
            found = False     #set to not start
            for line in data_file2: #looping throguh lines
                if found: #if start 
                    if line.strip() == 'endata': #if attend this stop
                        break
                    res_block += line.strip() + "\n"  #add line                    
                else: # foud is false
                    if line.strip() == 'resources:': #look for this condition
                        found = True #set as true
                        #tec_block = ""
            data_file2.close()    #close file     
            
            res1 = res_block[:] #copying it
            res = res1.split('*\n') #spliting each technology
            res = [d for d in res if d != ""]  #get free from empty elements
            
            for rr in res: #looping through resources
                tt1 =  rr.split('\n') #splitting by line
                tt1 = [t.strip().strip('\t') for t in tt1 if t != ""] 
                new_rr.append([node, dd[tt1[0].split(' ')[-1]].split('|')[-1]])
                if [node, dd[tt1[0].split(' ')[-1]].split('|')[-1]] in resource_df[['node','commodity']].values.tolist(): #check if commodity and node have been already defined

                    if any(r.startswith('grade') for r in tt1): #if garade it there
                        resource_df['grade'] = [r.split(' ')[-1] for r in tt1 if r.startswith('grade')][0] #replace by new grade
                    if any(r.startswith('volume') for r in tt1): #check if there is any volume
                        resource_df['value'] = [int(r.split('\t')[-1]) for r in tt1 if r.startswith('volume')][0] #replace value by volume
                else: #ow it is necessary to add new row into df
                    aux = pd.DataFrame({'node': node,
                                        'commodity': dd[tt1[0].split(' ')[-1]].split('|')[-1],
                                        'grade': tt1[2].split(' ')[-1],
                                        'value': int(tt1[3].split('\t')[-1])
                                        }, index=[0])
                    resource_df = pd.concat([resource_df, aux], ignore_index=True) #concat into df

        for old_rr in resource_df[['node','commodity']].values.tolist(): #looping through resource commoditties and nodes
            if old_rr not in new_rr: #check if old rr is in new rr
                resource_df = resource_df.loc[~((resource_df['commodity'] == old_rr[1]) & (resource_df['node'] == old_rr[0]))] #removing value that is not in ldb
        
        resource_df.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_resource_'+'df_'+ldb+'.csv'), index=False, sep=';') #exporting to csv                      
    return resource_df
#%%
# =============================================================================
# #DEMAND
# =============================================================================
def demand_conv(tpp, adb_files, year_act, count_y, ldb):
    """
        This function gets the info from the message adb file
        and organize it in tables and dictionaries
        a dictionary with all technologies constraints variables of each adb file (each region) is exported
        as well a csv file with all information of all technologies for each region
        The function return a dictionary with all technologies and their parameters in a list for each region
        It is necessary to define the path in which your projects are stored and the name of the project
    """
    print('\nworking on demand')
    task, path, project = tpp #getting data
    if task == 'ADB':
        ldb_aux = ldb
    else:
        ldb_aux = '_'+ldb
        
    #os.chdir(path)# set working drive   
    
    demand_df = pd.DataFrame(columns=['node', 'commodity', 'level', 'year_act', 'time', 'value', 'comment']) #aux demand df
    Level_form_dict1 = np.load(os.path.join(path,project,'IX',ldb,project+'_'+task+'_energy_forms_dict'+ldb_aux+'.npy'),allow_pickle='TRUE').item()  #lo
    for Reg in adb_files:#looping through adb files
        Level_form_dict = Level_form_dict1[Reg]  #load egergy form dictionary 

        dt = os.path.join(path, project, Reg,"data", Reg+".adb") #path file
        data_file2 = open(dt) #open file
        l = open(dt)#open file
        l1 = l.read() #read it
        l2 = l1.split('\n') #split by line
         
        node = [t for t in l2[0:10] if t.upper().startswith('ADB:')] #node
        node = node[0].split(':')[-1].strip(' ') #node
    
        l.close() #close file
        
        #Get the text between specific lines defined by its content
        demand_block = "" #generating block
        found = False #false        
        for line in data_file2: #looping through line
            if found: #if found is true               
                if line.strip() == 'loadcurve:': #if this condition is attended  
                    break 
                demand_block += line.strip() + "\n"  #ad line
            else: #ow found is false
                if line.strip() == 'demand:': # if line is this
                    found = True #set tture
   
        data_file2.close()         #close file
        
        demand1 = demand_block[:] #copying it
        demand = demand1.split('\n') #spliting each technology
        demand = [d for d in demand if d != ""]  #drop empty elements

        k = 0 #counting
        for n in demand: #looping through demand
            k+=1 #counting
            if n.startswith('#'): #check if its a comment
                comme = ('_').join(n.split('\t')) #add as coment
            else: #ow
                if len(n.split(' ')[0].split("\t")[0])>3: # if it does not follow the pattern X-D 
                    raise Exception('check Energy forms or demand energy forms of:\t '+"Region: "+Reg+"\n"+n) #raise exception
                else: #ow it is right
                    n1 = n.split('ts')[-1]
                    DEM = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", n1) #getting all float digits
                    DEM = DEM[count_y:] #starting from first year modelled
                    if len(DEM) == 0:
                        DEM = [0]
                    while len(DEM) < len(year_act): #check if length matches
                        DEM = DEM+[DEM[-1]] #make them match
                    DEM = [float(d) for d in DEM] #converting to float
                    
                    com_ = Level_form_dict[n.split(' ')[0].split("\t")[0]].split('|')[-1] #get commodity
                    lvl_ = Level_form_dict[n.split(' ')[0].split("\t")[0]].split('|')[0] #get level
                    tt_ = 'year' #set time as year
            if k%2==0: #if k is even it means that we have already gotten all necessary info
                for yy,dd in zip(year_act, DEM): #looping through year all and DEM
                    aux = expand_grid_name(['node', 'commodity', 'level', 'year_act', 'time', 'value', 'comment'], #df columns
                                           [node], [com_], [lvl_], [yy], [tt_], [dd], [comme]) #df values
                        
                    aux = aux[demand_df.columns] #ordering columns
                    demand_df = pd.concat([demand_df, aux], ignore_index=True) #concatenating df

    demand_df.to_csv(os.path.join(path,project,'IX',ldb,project+'_'+task+'_demand_df'+ldb_aux+'.csv'), index=False, sep=';')    #export to csv
            
    return demand_df
            
#%%
def capacity_factor_conversion( tpp, adb_files, ldr_file, time1, ldb=''):
    """
        This function gets the locad_curve info from the message adb.ldr file
        and organize it in tables representing the capacity factor for all techs with load_curve
        The function return the cf table

    """
    print('\nworking on capacity factor')
    task, path, project = tpp #getting data    
    if ldr_file:
        if time1 not in ['months', 'month-hours', 'original']: #check option defined
            raise Exception("time should be defined as 'months', 'month-hours' or 'original'")
    else:
        if time1 != 'original': #check option defined
            raise Exception("if you are not getting load region information from an ldr file, time must be equal 'original'")        
    
    #os.chdir(path)     # set working drive
    
    regdit = {} #aux dict
    for Reg in adb_files: #looping through regions
        if ldr_file: #if tis using ldr file to get ldr data
            if ldb == '': #check if it is adb or ldb
                dt = os.path.join(path, project,Reg,"data", Reg+"_adb.ldr") #ADB       
            else: #OW IT IS LDB
                dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+".ldr")  #LDB
        else: #if you are getting it from adb/ldb files
            if ldb == '': #check if it is adb or ldb
                dt = os.path.join(path, project,Reg,"data", Reg+".adb") #ADB       
            else: #OW IT IS LDB
                dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+".ldb")  #LDB
        try:
            data_file2 = open(dt) #open file
        except ValueError: #raise exception if you cannot load demand df from adb
            raise Exception("There is no ldr file or adb/ldb file in the current path:\t "+ dt)                 
        
        #Get the text between specific lines defined by its content
        lc_block = "" #cf block
        found = False #set starting
       
        if ldr_file:
            stop = ""
            start = 'loadcurves:'
        else:
            stop = "relationsc:"
            start = 'loadcurve:'            
        for line in data_file2: #looping through lines
            if found: #start
                if line.strip() == stop: #stop in this case
                    break
                lc_block += line.strip() + "\n"  #add line
            else: #if found is false
                if line.strip() == start: #start condition
                    #print('Y')
                    found = True #set as true

        data_file2.close()   #close file 
        
        lc = lc_block[:] #str has no attribute copy   

        # =============================================================================
        lc1 = lc.split('\n') #split by line
        if ldr_file:
            lc2 = [s for s in lc1 if s.startswith('systems')] #get systems names that represents the tec with ldr
        else:
            lc2 = [s.split(' ')[0] for s in lc1 if s.startswith('systems')] #get systems names that represents the tec with ldr            
        lc2 = lc2+["_end_"] #add an end value
        lc = lc+'_end_'  #add an end value
        
        load_tecs = {}
        for k in range(len(lc2)-1):
            load_tecs[lc2[k]] = find_between(lc, lc2[k], lc2[k+1])
        
        m_load_tec = {} #aux dict 2
        for lcn in lc2[:-1]: #looping through name of ldr tecs
            #print(lcn)

            ll = [l for l in load_tecs[lcn].split('\n') if l != ""] #split by line and get rid of empty ones
            if ldr_file:
                ldr = [] #aux list 
                if len(ll[1].split(' ')) >12 or len(ll[1].split(' ')) != len([l for l in ll if l.startswith('1.0')]):
                    raise Exception ("Check your ldr data for:\t "+ dt) #RAISE EXCEPTION IF LEN OF MAJOR VALUE IS BIGGER THAN 12 OR IF LENGTH OF MAJOR TEMPORAL VALUE DOES NOT MATCH NUMBER OF MINOR TEMPORAL VALUE
                ll_aux = [l for l in ll[1:] if not l.startswith('1.0')] #EXCLUDING 1.0000
                ll_aux_M_o = [float(l) for l in ll_aux[0].split(' ')] #MAJOR TEMPORAL VALUE
                ll_aux_m_o = [[float(a) for a in lla.split(' ')] for lla in ll_aux[1:] ] #MINOR TEMPORAL VALUE   
                ll_aux_m = {}
                
                if time1 == 'month-hours': #check if its month hours
                    if len(ll_aux_M_o) <=12 and 12%len(ll_aux_M_o) ==0: # check if it is possible to convert original data from major temporal value to 12 months
                        na = int(12/len(ll_aux_M_o)) #aux value to extend to 12 values
                        ll_aux_M = list(itertools.chain(*[[l/na]*na for l in ll_aux_M_o])) #extending for 12 values
                    else: #ow raise exception bc it is not possible to take advantage of original data
                        raise Exception("It is not possible to do a month-hours (24h for each month in a total of 288 time steps) time from this database\ncheck major ldr levels")
                    for llm in range(len(ll_aux_m_o)): #looping through length of minor temporal resolution
                        if len(ll_aux_m_o[llm]) <=24 and 24%len(ll_aux_m_o[llm]) ==0: # check if it is possible to convert original data form major temporal value to 24 hours
                            na = int(24/len(ll_aux_m_o[llm])) #aux value to extend to 24 values
                            ll_aux_m[llm] = list(itertools.chain(*[[l/na]*na for l in ll_aux_m_o[llm]])) #extending for 24 values
                        else: #ow raise exception bc it is not possible to take advantage of original data
                            raise Exception("It is not possible to do a month-hours (24h for each month in a total of 288 time steps) time from this database\ncheck minor ldr levels")

                #if time1 != "months":    # if it is not month                                              
                    #ldr = [] #aux list 
                    fc = int(len(ll_aux_M)/len(ll_aux_M_o))
                    for n in range(len(ll_aux_M_o)): #LOOPING THROUGH LEN OF MAJOR TEMPORAL VALUE            
                        ll_aux_m1 = np.array(ll_aux_m[n]) * ll_aux_M[n] / round(sum(ll_aux_m[n]),5) #adjusting minor temporal resolution values to match MAJOR RESOLUTION VALUE
                        m1_aux = list(ll_aux_m1) * fc # adjusting to match number of MAJOR values
                        ll_aux_m1 = m1_aux[:]
                        ldr.append(ll_aux_m1) #appending into ldr
                    ldr = [round(l,8) for l in list(itertools.chain(*ldr))] #adjusting ldr by flatteninfg the list
                    
                elif time1 == "months": # if it is month
                    if len(ll_aux_M) <12 and 12%len(ll_aux_M) ==0: # check if it is possible to convert original data form major temporal value to 12 months
                        na = int(12/len(ll_aux_M)) #aux value to extend to 12 values
                        ll_aux_M = list(itertools.chain(*[[l/na]*na for l in ll_aux_M]))   #extending for 12 values
                        ldr = ll_aux_M.copy() #copying major ldr to ldr
                        time = [31/365, 28/365, 31/365, 30/365, 31/365, 30/365, 31/365, 31/365, 30/365, 31/365, 30/365, 31/365] #time for month
                    else:#raise exception in case it is not possible to convert original data to month
                        raise Exception("It is not possible to do a months (JAN to DEC with 12 time steps) time from this database\ncheck major ldr levels")
                
                else:# if it is original
                    for n in range(len(ll_aux_M_o)):
                        ll_aux_m1 = np.array(ll_aux_m_o[n]) * ll_aux_M_o[n] / round(sum(ll_aux_m_o[n]),5) #adjusting minor temporal resolution values to match MAJOR RESOLUTION VALUE
                        ldr.append(ll_aux_m1) #appending into ldr
                    ldr = [round(l,8) for l in list(itertools.chain(*ldr))] #adjusting ldr by flatteninfg the list
                        
                lcn2 = "_".join([fn[0].upper()+fn[1:].lower() for fn in lcn.split('.')[1].split("_")])#lcn.split('.')[1] #adjusting name
                lcn2 = lcn2.strip()
                m_load_tec[lcn2] = ldr #get all ldr into a dictionary
            else:
                ll = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", "".join(ll))
                ldr  = [float(l) for l in ll]
                lcn2 = "_".join([fn[0].upper()+fn[1:].lower() for fn in lcn.split('.')[1].split("_")])#lcn.split('.')[1] #adjusting name
                lcn2 = lcn2.strip()
                m_load_tec[lcn2] = ldr
                
        if time1 != 'original':
            time = [31/365, 28/365, 31/365, 30/365, 31/365, 30/365, 31/365, 31/365, 30/365, 31/365, 30/365, 31/365] #time for month
            if time1 == 'month-hours':
                time = list(itertools.chain(*[[t/24]*24 for t in time])) #time for month hours
                time = [round(t,8) for t in time] #round time values
        else:
            if ldr_file: #check if ldr file is true and get time from it
                data_file2 = open(dt) #open file              
            
                #Get the text between specific lines defined by its content
                time_block = "" #cf block
                found = False #set starting
               
                for line in data_file2: #looping through lines
                    if found: #start
                        if line.strip() == 'loadcurves:': #stop in this case
                            break
                        if line.strip().startswith("day"):
                            time_block += line.strip() + "\n"  #add line
                        if line.strip().startswith("length"):
                            time_block += line.strip() + "\n"  #add line                        
                    else: #if found is false
                        if line.strip() == 'loadregions:': #start condition
                            #print('Y')
                            found = True #set as true
        
                data_file2.close()   #close file 
                
                tb1 = time_block[:] #str has no attribute copy      
                tb = [t for t in tb1.split('\n') if t != ""] #spliting by line and removing empty ones
                tb_day = [t.split('day')[-1] for t in tb if t.startswith('day')] #get the ones related to the day
                tb_len = [re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", t.split('length')[-1]) for t in tb if t.startswith('length')] #get the ones related to the length
    
                tb_day = [float(t.strip()) for t in tb_day] #day values adjusted
                tb_len = [[float(t) for t in tt] for tt in tb_len] #length values adjusted
                
                tb_d_t = list(np.array(tb_day)/365) #day values representation of the year
                
                time = [list(np.array(l) *d/ sum(l)) for d,l in zip(tb_d_t, tb_len)] #time adjustment for the number of periods in the database
                time = list(itertools.chain(*time)) #make list flat
                time = [round(t,8) for t in time]#rounding it
                
            else: #get time from adb/ldb database
                data_file2 = open(dt) #open file              
            
                #Get the text between specific lines defined by its content
                time_block = "" #cf block
                found = False #set starting
               
                for line in data_file2: #looping through lines
                    if found: #start
                        if line.strip() == 'energyforms:': #stop in this case
                            break
                        time_block += line.strip() + "\n"  #add line                        
                    else: #if found is false
                        if line.startswith('length'): #start condition
                            time_block += line #add line
                            found = True #set as true
        
                data_file2.close()   #close file 
                
                time = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", time_block) #time from adb/ldb database
                time = [float(t) for t in time]
                
        def NormalizeData(data):#normalize to the max value
            return (data ) / (np.max(data)) 
        
# =============================================================================
#         def NormalizeData_mean(data):
#             return (data ) / (np.mean(data))
# =============================================================================
        
        #regdit[Reg+'_cf'] = {}
        regdit[Reg+'_cf_norm'] = {} #aux df
        #regdit[Reg+'_cf_mean'] = {}
        for t in m_load_tec.keys(): #looping through keys that represent each tec name
            #regdit[Reg+'_cf'][t] = np.array(m_load_tec[t])/np.array(time)
            cf = np.array(m_load_tec[t])/np.array(time) #get capacity factor associated with this ldr by dividing each ldr value by time representation
            #regdit[Reg+'_cf_norm'][t] = NormalizeData(regdit[Reg+'_cf'][t])
            t = t[0].upper() + t[1:] #adjusting name
            regdit[Reg+'_cf_norm'][t] = NormalizeData(cf) #normalize data according to maximum value to get it right as capacity factor
            #regdit[Reg+'_cf_mean'][t] = NormalizeData_mean(regdit[Reg+'_cf'][t])
            ldr_info = regdit
    #np.save(os.path.join(path,project,'ldr_info_'+ldb), regdit) #saving dictionary

    return ldr_info
#%%
def tec_link(BLUES_tec_df, tec_list, tpp, ldb=''):
    """
    tec link is a function that get the link between the technologies by listing their inputs and outputs
    """
    startTime = datetime.now() 
    print('\nworking on tec link chain')
    if ldb == "":
        ldb_aux = ldb
    else:
        ldb_aux = '_'+ldb
    task, path, project = tpp #getting data            
    #os.chdir(path) # set working drive
    #idx = pd.IndexSlice #pandas slice
    i_check_list = [] # tecs that commoditties and level do not have same length
    o_check_list = [] # tecs that commoditties and level do not have same length

    #link_t = {} #aux dict
    #link_t =  {key: {'input_for_tec':None, 'output_from_tec':None} for key in tec_list} #aux dict
    
    dft = BLUES_tec_df[['input_levels', 'input_commodities', 'output_levels', 'output_commodities']].copy() #GET ONLY INPUT OUTPUT COLUMNS
    #dft = dft.fillna('nan')  #REPLACING NAN
    dfcc = df_generator() 
    dfcc_c = [c for c in dfcc.columns.tolist() if c in dft.columns.tolist() ] 
    #dft.sort_index( axis=1, level = 1, inplace=True)
    dft = dft[dfcc_c] #ORDERING COLUMNS

    #GET INPUT TECS
    def input_f(t2):
        if hasNumbers(t2[-4:]): #checking it all the last 4 characteres are number --> if it is True it means that they are existing technologies installed before 2008
            t = t2[:-5] #removing the last 4 characteres to get the parameters
        else:
            t = t2 #ow it is not necessary        
        # GET INPUT TECS
        #i_ldr_lvl_com2 = [] #input auxiliary list
        # GET t commodities and levels
        i_ll = dft[dft.index.get_level_values(0)==t]['input_levels'].values.tolist() #getting levels
        i_ll2 = list(itertools.chain(*i_ll)) #flattening the list 
        i_ll2 = [ll for ll in i_ll2 if ll != 'nan'] #excluding nan values
        
        i_cc = dft[dft.index.get_level_values(0)==t]['input_commodities'].values.tolist() #getting commodities
        i_cc2 = list(itertools.chain(*i_cc)) #flattening the list
        i_cc2 = [cc for cc in i_cc2 if cc != 'nan'] #removing nan
        
        if len(i_ll2) == len(i_cc2): #check if levels and commodities have same length
            i_llcc = [(vl, vc) for vl,vc in zip(i_ll2, i_cc2)] #list of pair lvl com
            i_llcc = list(set(i_llcc)) #setting list
    
        else: #they must have smae length
            print('check')
            i_check_list.append(t) #append into i list
        return i_llcc    #return level-commodty input of each tech
    
    aai = list(map(input_f, tec_list)) #map applies the function for each elemnt in the list

    #GET  TECS THAT HAVE OUTPUT THE INPUT LVL-COM    
    def input_f2(i_llcc):
        inp_t2_f = []
        if len(i_llcc) >0:
            for lc2 in i_llcc:
                #print(lc2)
                #append tecs that have as output level and com ay of the lvl com pair defined
                i_new_tecs2 = dft[(dft['output_levels'] == lc2[0]) & (dft['output_commodities'] == lc2[1])] 
                i_new_tecs2.dropna(how='all', inplace=True) #removing nan rows from datafram to get only the ones with data
                inp_t = i_new_tecs2.index.get_level_values(level=0).tolist() #getting tecs
                inp_t2 = list(set(inp_t)) #removing duplicate tecs
                inp_t2_f.append(inp_t2)
        
        inp_t2_f = list(unique(list(itertools.chain(*inp_t2_f))))
        return inp_t2_f #return tecs that have as output the level-commodity pair defined
    
    ii = list(map(input_f2, aai)) #map applies the function for each elemnt in the list
    
    #GET OUTPUT TECS
    def output_f(t2):
        if hasNumbers(t2[-4:]): #checking it all the last 4 characteres are number --> if it is True it means that they are existing technologies installed before 2008
            t = t2[:-5] #removing the last 4 characteres to get the parameters
        else:
            t = t2 #ow it is not necessary             
        # GET OUTPUT TECS
        #i_ldr_lvl_com2 = [] #input auxiliary list
        # GET t commodities and levels
        o_ll = dft[dft.index.get_level_values(0)==t]['output_levels'].values.tolist() #getting levels
        o_ll2 = list(itertools.chain(*o_ll)) #flattening the list 
        o_ll2 = [ll for ll in o_ll2 if ll != 'nan'] #excluding nan values
        
        o_cc = dft[dft.index.get_level_values(0)==t]['output_commodities'].values.tolist() #getting commodities
        o_cc2 = list(itertools.chain(*o_cc)) #flattening the list
        o_cc2 = [cc for cc in o_cc2 if cc != 'nan'] #removing nan
        
        if len(o_ll2) == len(o_cc2): #check if levels and commodities have same length
            o_llcc = [(vl, vc) for vl,vc in zip(o_ll2, o_cc2)] #list of pir lvl com
            o_llcc = list(set(o_llcc)) #setting list
    
        else: #they must have smae length
            print('check')
            o_check_list.append(t) #append into i list
        return o_llcc   #return level-commodty output of each tech
    
    aao = list(map(output_f, tec_list)) #map applies the function for each elemnt in the list

    #GET  TECS THAT HAVE INPUT THE OUTPUT LVL-COM
    def output_f2(o_llcc):
        out_t2_f = []
        if len(o_llcc) >0: 
            for lc2 in o_llcc:
                #print(lc2)
                #append tecs that have as output level and com ay of the lvl com pair defined
                o_new_tecs2 = dft[(dft['input_levels'] == lc2[0]) & (dft['input_commodities'] == lc2[1])] 
                o_new_tecs2.dropna(how='all', inplace=True) #removing nan rows from datafram to get only the ones with data
                out_t = o_new_tecs2.index.get_level_values(level=0).tolist() #get tecs
                out_t2 = list(set(out_t)) #removing duplicate tecs
                out_t2_f.append(out_t2)

        out_t2_f = list(unique(list(itertools.chain(*out_t2_f))))    
        return out_t2_f #return tecs that have as input the level-commodity pair defined
    
    oo = list(map(output_f2, aao)) #map applies the function for each elemnt in the list
    
    aa_ii_oo = aai + aao
    aa_ii_oo2 = list(itertools.chain(*aa_ii_oo))
    aa_ii_oo3 =  list(set([ tuple(sorted(t)) for t in aa_ii_oo2 ]))
    aa_ii_oo4 = sorted(aa_ii_oo3, key=lambda tup: tup[0])
    aa_ii_oo5 = [[att] for att in aa_ii_oo4]
    ii2 = list(map(input_f2, aa_ii_oo5)) #map applies the function for each elemnt in the list
    oo2 = list(map(output_f2, aa_ii_oo5)) #map applies the function for each elemnt in the list
    
    link_t = {key:{'0get_input_from_tec': i, '1gen_output_for_tec':o} for key,i,o in zip(tec_list, ii, oo)} #link of technologies
    lvl_com_t = {key:{'tec_use_as_input': i, 'tec_gen_output':o} for key,i,o in zip(aa_ii_oo4,  oo2, ii2)} #link of technologies
    

    np.save(os.path.join(path,project,'IX',ldb,project+'_'+task+'_link_t_info'+ldb_aux), link_t) #save it 
    np.save(os.path.join(path,project,'IX',ldb,project+'_'+task+'_lvl_com_t_info'+ldb_aux), lvl_com_t) #save it
    
    dt_final =  datetime.now() - startTime #finish time
    print('Done.')
    print('Runtime:    '+ str(dt_final))
    return link_t, lvl_com_t
#%%
def ldr_tecs(ldr_info, BLUES_tec_df, tpp, adb_files, link_t, ldb):
    print('\nworking on ldr tecs')
    task, path, project = tpp #getting data    
    ldr_tec = [] #aux list
    for l in ldr_info.keys(): #looping through ldr info for each region
        for l1 in ldr_info[l].keys(): #looping inside for each technology with ldr
            if l1 not in ldr_tec: #check if it is not in ldr tec
                ldr_tec.append(l1) #append in case it is not
     
    ldr_tec_0 = ldr_tec.copy() #make a copy

    #############################################################################################################
    ldr_list = [] #aux ldr list

    for Reg in adb_files: #looping through regions
        
        if task == 'ADB':
            dt = os.path.join(path, project, Reg,"data",Reg+'.adb') #path - energy forms
        elif task == 'LDB':
            dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+'.ldb') #path - energy forms
        
        dt1 = open(dt) #reading adb
        
        E_form_block = "" #creating block
        found = False #starting
        for line in dt1: #looping through lines
            if found: #when found is True add line
                if line.strip() == 'demand:':  #stop line
                    break
                E_form_block += line.strip() + "\n"  #addgin line
            else: #found is false
                if line.strip() == 'energyforms:': #starting line
                    found = True #set found = True

        dt1.close() #close file

        E_forms1 = E_form_block[:] #copying E form
        E_forms = E_forms1.split('*') #splitting each energy form
        E_forms = [e for e in E_forms if e != '' if e != '\n'] #removing empty lines 
                   
        for n in E_forms: #looping through energy forms
            f1 = ''.join(n.split('\t')).split('\n') #splitting by tab, then joining then splitting by line
            f1 = [t for t in f1 if t!=''] #excluding empty lines
            LVL = f1[0].split(' ') #getting the commodity
            for l in f1[1:]: #looping through levels
                l = l.strip(' ') #removing empty characteres in the beggining and in the end of str
                if not l.startswith('#'): #if its not a comment
                    COM = l.split(' ') #split by empty characteres
                    if len(COM) >2 and COM[-1] == 'l': #get commodities that are listed as load factor
                        ldr_list.append((LVL[0],COM[0].lower())) #append ldr commodities into ldr list

    ldr_list = list(set(ldr_list)) #remove duplicates
    ldr_tecs = [] #aux ldr tec list
    for lc in ldr_list: #looping through ldr list of com-lvl
        lvl = lc[0] #level
        com = lc[1] #com
        kilen = 1 + len(BLUES_tec_df['input_commodities'].columns.values.tolist())
        kolen = 1 + len(BLUES_tec_df['output_commodities'].columns.values.tolist())
        #append tecs that have as input/output level and com any of the lvl com pair defined
        [ldr_tecs.append(BLUES_tec_df[(BLUES_tec_df[idx['input_levels', str(k)]]== lvl) & (BLUES_tec_df[idx['input_commodities', str(k)]]== com)].index.get_level_values('tec').values.tolist()) for k in list(range(1,kilen))]    
        [ldr_tecs.append(BLUES_tec_df[(BLUES_tec_df[idx['output_levels', str(k)]]== lvl) & (BLUES_tec_df[idx['output_commodities', str(k)]]== com)].index.get_level_values('tec').values.tolist()) for k in list(range(1,kolen))]            
    
    ldr_tecs = list(itertools.chain(*ldr_tecs))    #flattening list
    ldr_tecs = sorted(list(set(ldr_tecs)))     #remove duplicates
    
    check_cf_ldr = [a for a in ldr_tec_0 if a not in ldr_tecs] #check if ldr info from capacity factor matches ldr tecs
    if len(check_cf_ldr)>0: #in case not 
        raise Exception("Check your data because the following technologies have ldr info that do not match commodity-level ldr info\n"+str(check_cf_ldr))

    ldr_tec_1 = ldr_tecs.copy() #copy of ldr info
    #ldr_tec_dict = {k: v for k, v in link_t.items() if k in ldr_tec_o} #getting link of technologies for ldr tecs
    
    #ldr_d =  [ l for l in ldr_tec_1 if l not in ldr_tec_o]        #check ldr tecs that were added
    #ldr_tec_dict1 = {k: v for k, v in link_t.items() if k in ldr_tec_1}  #getting link of technologies for new ldr tecs
    
    int_ldr_tec = [] #TECS THAT HAVE defined time in and out as time (1-12, 1-288....) 
    final_ldr = [] #output will be year
    first_ldr = [] #input will be year
    ldr_tec_1_aux = ldr_tec_1.copy() #copy of ldr info
    
    while len(ldr_tec_1_aux) >0: #check if aux ldr info is >0
        for ldr in ldr_tec_1_aux: #looping through ldrs tecs
            ldr_tec_1_aux = [t for t in ldr_tec_1_aux if t != ldr] #keep ldr tecs that are different from the current ldr in the loop    
            
            ldr_out = link_t[ldr]['1gen_output_for_tec'] #get output of ldr tec     link_t = {key:{'0get_input_from_tec': i, '1gen_output_for_tec':o} for key,i,o in zip(tec_list, ii, oo)} #link of technologies
            aux_out = any(elem in ldr_tec_1 for elem in ldr_out) #check if there is any output tec in ldr 
            
            ldr_inp = link_t[ldr]['0get_input_from_tec'] #get input of ldr tec
            aux_inp = any(elem in ldr_tec_1 for elem in ldr_inp) #check if there is any input tec in ldr 
                   
            if not aux_out: #it means that there is no output in LDR
                final_ldr.append(ldr) #So this is a final ldr tec
    
        # =============================================================================
        # There are two different approaches for ldr in older versions
        # =============================================================================
                        
            if not aux_inp: #it means that there is no input in LDR
                first_ldr.append(ldr) #So ldr is a first ldr technology

    #ldr_tec_dict1 = {k: v for k, v in link_t.items() if k in ldr_tec_1} #ldr tec link dictionary of input and output tecs
    #ldr_tec_dict_final = {k: v for k, v in link_t.items() if k in final_ldr}     #final ldr tec link dictionary of input and output tecs
    #ldr_tec_dict_first = {k: v for k, v in link_t.items() if k in first_ldr}  #first ldr tec link dictionary of input and output tecs
    
    mid_ldr = [m for m in ldr_tec_1 if m not in final_ldr and m not in first_ldr and m not in int_ldr_tec]    #get mid ldr, the ones that have as input and output time as defined
    #ldr_tec_dict_mid = {k: v for k, v in link_t.items() if k in mid_ldr} #mid ldr tec link dictionary of input and output tecs
    
    # =============================================================================
    # FOR CHECKING THE MID LDR TECS
    # mid_inp = [ldr_tec_dict_mid[m]['input_tec'] for m in mid_ldr] #list of input tec of the mid ldr tecs
    # mid_inp = list(itertools.chain(*mid_inp)) #flattening the list
    # [m for m in mid_inp if m not in ldr_tec_1] # check if any of the values is not a ldr tec
    # all(elem in ldr_tec_1 for elem in mid_inp)  # check if any of the values is not a ldr tec
    # 
    # mid_out = [ldr_tec_dict_mid[m]['output_tec'] for m in mid_ldr] #list of output tec of the mid ldr tecs
    # mid_out = list(itertools.chain(*mid_out)) #flattening the list
    # [m for m in mid_out if m not in ldr_tec_1] # check if any of the values is not a ldr tec
    # all(elem in ldr_tec_1 for elem in mid_out)        # check if any of the values is not a ldr tec     
    # =============================================================================
    return final_ldr, first_ldr, mid_ldr, ldr_tecs, ldr_list
  
#%%
def Relation_constraints(tpp, adb_files,  ldb='', land=False, emission=False):
    """
    This script goes over all relations defined in the adb/ldb database
    relations might be splitted into land, emission and others which is the default
    """

    if ldb == "":
        ldb_aux = ldb
    else:
        ldb_aux = '_'+ldb
        
    if land: #if land load land relation list
        print('\nworking on land constraints')
    elif emission: #if emission load emission relation list
        print('\nworking on emission constraints')
    else: #ow load other relation list
        print('\nworking on constraints')

    task, path, project = tpp
    #os.chdir(path)     # set working drive
    relations = ['relationsc:', 'relationsp:', 'relationss:', 'relations1:', 'relations2:','variables:'] #relations type
    if task == 'ADB':   #check ldb
        rel_df = pd.DataFrame(columns=['node', 'relation', 'condition', 'cost', 'unit_group', 'ldr', 'lower', 'upper', 'type', 'comment']) #aux df
       
        for Reg in adb_files: #looping through regions
            for r in range(len(relations[:-1])): #looping through type of relations
                
                dt = os.path.join(path, project, Reg,"data", Reg+".adb") #adb path    
                dt = open(dt) #open file
                rel_block = "" #generating block
                found = False #set atart
                for line in dt: #looping through lines
                    if found: #is start
                        if line.strip() == relations[r+1]: #stop if it reaches the next relation type
                            break
                        rel_block += line.strip() + "\n"  #add line
                    else: #ow found is false
                        if line.strip() == relations[r]: #check start condition
                            found = True    
                dt.close() #close file
                
                rel_blk1 = rel_block[:] #str has noattribute copy
                rel_blk = rel_blk1.split('*\n') #spliting each technology
                #rel_blk = rel_blk[:-1]
                #rel_b_name = [r.split(' ')[0].split('\t')[0] for r in rel_blk] #rel name os blocks
                
                if land: #if land load land relation list
                    with open(os.path.join(path, project,'IX',project+'_'+task+'_land_relation_list'), 'rb') as fp:
                        rel_ = pickle.load(fp)  
                elif emission: #if emission load emission relation list
                    with open(os.path.join(path, project,'IX',project+'_'+task+'_emission_list'), 'rb') as fp:
                        rel_ = pickle.load(fp)     
                else: #ow load other relation list
                    with open(os.path.join(path, project,'IX', project+'_'+task+'_relation_list'), 'rb') as fp:
                        rel_ = pickle.load(fp)
                        
                rel_1 = sorted(list(set([r.split('-')[0] for r in rel_]))) #get all rels minus the node info of them
                rel = [r for r in rel_blk if r.startswith(tuple(rel_1))] #keep relation of the block that are listed in technologies
                for nn in rel: #looping through relations
                    f1 = ' '.join(nn.split('\t')).split('\n') #adjusting based on tab and lines
                    f1 = [t for t in f1 if t!=''] #removing empty lines
                    f1 = [t[1:] if t[0]==' ' else t for t in f1] #removing first character if it is empty   
                    f1 = [f.strip() for f in f1] #removing empty characteres
                    rel_name = f1[0].split(' ')[0] #get relation name
                    rel_condition = f1[0].split(' ')[-1] #get relation name
                    if any(r.startswith('cost') for r in f1): #look for cost info into relation block
                        cost = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('cost')][0]) #getting float digits 'd*.?d+'
                    else: #else set cost as nan
                        cost = [['nan']]

                    ug = re.findall('(?:: )(.*)', [f for f in f1 if f.startswith('units group')][0])[0].split(',') #unit groups info
                    ug = str([[ u.strip() for u in ug]]) #removing empty lines in the beggining and in the end
                    ldr = [f for f in f1 if f.startswith('for_ldr')][0].split(' ')[-1] #get ldr info
                    if any(r.startswith('upper') for r in f1): #check if there upper info
                        up = str([re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('upper')][0])]) #getting all float digits
                    else: #ow set nan
                        up = str(['nan']) #
                    if any(r.startswith('lower') for r in f1): #check if there lower info
                        lw = str([re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('lower')][0])])  #getting all float digit #'-?d*.?d+'
                    else: #ow set nan
                        lw = str(['nan'])
                    if any(r.startswith('type') for r in f1): #look for type info
                        tp = [f for f in f1 if f.startswith('type')][0].split(' ')[-1]
                    else: #ow set nan
                        tp = ['nan']
                    if any(r.startswith('#') for r in f1): #look for comments
                        cm = f1[-1].split('#')[-1]#[1:]
                    else: #ow set nan
                        cm = ['nan']
                    
                    #generating dataframe
                    aux = pd.DataFrame({'node':Reg,
                                        'relation':rel_name,
                                        'condition':rel_condition,
                                        'cost': str(cost),
                                        'unit_group': ug,
                                        'ldr': ldr,
                                        'lower':lw,
                                        'upper':up,
                                        'type': tp,
                                        'relation_type':relations[r],
                                        'comment':cm}, index=[0])
                        
                    rel_df = pd.concat([rel_df, aux], ignore_index=True) #concatenating into rel_df

    elif task == 'LDB': #check if ldb supposes to be defined
        try: #try load rel_df from adb
            if land:
                rel_df = pd.read_csv(os.path.join(path,project,'IX', project+'_ADB_land_relation_df.csv'),sep=';', dtype='object')
            elif emission:
                rel_df = pd.read_csv(os.path.join(path,project,'IX', project+'_ADB_emission_relation_df.csv'),sep=';', dtype='object')        
            else:
                rel_df = pd.read_csv(os.path.join(path,project,'IX', project+'_ADB_relation_df.csv'), sep=';', dtype='object')
        except ValueError: #raise exception if you cannot load demand df from adb
            raise Exception("There is no rel df (or emission|land reld df) in the current path:"+ os.path.join(path,project,project)+"\n\ntry to set ldb = '' or check your pathes") 
        
        new_rr = []
        for Reg in adb_files: #looping through regions
            for r in range(len(relations[:-1])):#looping through relation
    
                dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+".ldb")  #path
                dt = open(dt)#open file
                rel_block = "" #genreating block
                found = False#set starting
                for line in dt: #looping through lines
                    if found:   #check start condition
                        if line.strip() == relations[r+1]:  #stop loop
                            break
                        rel_block += line.strip() + "\n"  #add line
                    else: #found is false
                        if line.strip() == relations[r]: #check condition to set start
                            found = True #start

                #print(rel_block)
                dt.close() #close file
                
                rel_blk1 = rel_block[:] #str has no attribute copy
                rel_blk = rel_blk1.split('*\n') #spliting each technology
                #rel_blk = rel_blk[:-1]
                #rel_b_name = [r.split(' ')[0].split('\t')[0] for r in rel_blk] #rel name os blocks
                
                if land: #load relation
                    with open(os.path.join(path, project,'IX',ldb,project+'_'+task+'_land_relation_list_'+ldb), 'rb') as fp:
                        rel_ = pickle.load(fp)
                elif emission:
                    with open(os.path.join(path, project,'IX',ldb,project+'_'+task+'_emission_list_'+ldb), 'rb') as fp:
                        rel_ = pickle.load(fp) 
                else:
                    with open(os.path.join(path, project,'IX',ldb,project+'_'+task+'_relation_list_'+ldb), 'rb') as fp:
                        rel_ = pickle.load(fp)
                        
                rel_1 = sorted(list(set([r.split('-')[0] for r in rel_]))) #get all rels minus the node info of them
                #[r for r in rel_1 if r not in rel_b_name] #for checking if there are missing relations
                rel = [r for r in rel_blk if r.startswith(tuple(rel_1))] #defining relation that were found in tec parameters
                for nn in rel: #looping throuhg relations
                    f1 = ' '.join(nn.split('\t')).split('\n') #adjusting based on tab and lines
                    f1 = [t for t in f1 if t!=''] #removing empty lines
                    f1 = [t[1:] if t[0]==' ' else t for t in f1] #removing first character if it is empty   
                    f1 = [f.strip() for f in f1] #removing empty characteres
                    rel_name = f1[0].split(' ')[0] #get relation name
                    rel_condition = f1[0].split(' ')[-1] #get relation name
                    new_rr.append([Reg, rel_name])
                    
                    if [Reg, rel_name] in rel_df[['node', 'relation']].values.tolist(): #check if relation-node pair is already in the adb rel df 
                        
                        if any(r.startswith('cost') for r in f1):#look for cost info
                            cost = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('cost')][0]) #getting float digits 'd*.?d+'
                            rel_df.loc[(rel_df['node'] == Reg ) & (rel_df['relation'] == rel_name), 'cost'] = str(cost) #set new value
                        if any(r.startswith('units') for r in f1):#look for units group
                            ug = re.findall('(?:: )(.*)', [f for f in f1 if f.startswith('units group')][0])[0].split(',') #unit groups info
                            ug = [[ u.strip() for u in ug]] #removing empty lines in the beggining and in the end
                            rel_df.loc[(rel_df['node'] == Reg ) & (rel_df['relation'] == rel_name), 'unit_group'] = str(ug) #set new value
                        if any(r.startswith('for_ldr') for r in f1):
                            ldr = [f for f in f1 if f.startswith('for_ldr')][0].split(' ')[-1] #get ldr info
                            rel_df.loc[(rel_df['node'] == Reg ) & (rel_df['relation'] == rel_name), 'ldr'] = str(ldr) #set new value
                        if any(r.startswith('upper') for r in f1): #look for upper info
                            up = str([re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('upper')][0])]) #getting all float digits
                            rel_df.loc[(rel_df['node'] == Reg ) & (rel_df['relation'] == rel_name), 'upper'] = str(up) #set new value
                        if any(r.startswith('lower') for r in f1):
                            lw = str([re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('lower')][0])]) #getting all float digits
                            rel_df.loc[(rel_df['node'] == Reg ) & (rel_df['relation'] == rel_name), 'lower'] = str(lw)
                        if any(r.startswith('#') for r in f1):#look for comments
                            cm = f1[-1].split('#')[-1]#[1:] #get comments
                            rel_df.loc[(rel_df['node'] == Reg ) & (rel_df['relation'] == rel_name), 'comment'] = str(cm) #set new value
                        if any(r.startswith('type') for r in f1): #lookf for type info
                            tp = [f for f in f1 if f.startswith('type')][0].split(' ')[-1] #get type info
                            rel_df.loc[(rel_df['node'] == Reg ) & (rel_df['relation'] == rel_name), 'type'] = str(tp) #set new value
                    else: #relation-node pair is not in the rel_df from adb
                        if any(r.startswith('cost') for r in f1): #look for cost info into relation block
                            cost = re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('cost')][0]) #getting float digits 'd*.?d+'
                        else: #else set cost as nan
                            cost = [['nan']]
    
                        ug = re.findall('(?:: )(.*)', [f for f in f1 if f.startswith('units group')][0])[0].split(',') #unit groups info
                        ug = str([[ u.strip() for u in ug]]) #removing empty lines in the beggining and in the end
                        ldr = [f for f in f1 if f.startswith('for_ldr')][0].split(' ')[-1] #get ldr info
                        if any(r.startswith('upper') for r in f1): #check if there upper info
                            up = str([re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('upper')][0])]) #getting all float digits
                        else: #ow set nan
                            up = str(['nan']) #
                        if any(r.startswith('lower') for r in f1): #check if there lower info
                            lw = str([re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", [f for f in f1 if f.startswith('lower')][0])])  #getting all float digit #'-?d*.?d+'
                        else: #ow set nan
                            lw = str(['nan'])
                        if any(r.startswith('type') for r in f1): #look for type info
                            tp = [f for f in f1 if f.startswith('type')][0].split(' ')[-1] #getting type
                        else: #ow set nan
                            tp = ['nan']
                        if any(r.startswith('#') for r in f1): #look for comments
                            cm = f1[-1].split('#')[-1]#[1:] #getting comments
                        else: #ow set nan
                            cm = ['nan']
                        
                        #generating df
                        aux = pd.DataFrame({'node':Reg,
                                            'relation':rel_name,
                                            'condition':rel_condition,
                                            'cost': str(cost),
                                            'unit_group': ug,
                                            'ldr': ldr,
                                            'lower':lw,
                                            'upper':up,
                                            'type': tp,
                                            'relation_type':relations[r],
                                            'comment':cm}, index=[0])
                                
                        rel_df = pd.concat([rel_df, aux], ignore_index=True) #concatenating df

        for old_rr in rel_df[['node','relation']].values.tolist(): #looping through resource commoditties and nodes
            if old_rr not in new_rr: #check if old rr is in new rr
                rel_df = rel_df.loc[~((rel_df['relation'] == old_rr[1]) & (rel_df['node'] == old_rr[0]))] #removing value that is not in ldb

    
    if land: #exporting data according to the relation input
        rel_df.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_land_relation_df'+ldb_aux+'.csv'), index=False, sep=';')
    elif emission:
        rel_df.to_csv(os.path.join(path,project,'IX', ldb,project+'_'+task+'_emission_relation_df'+ldb_aux+'.csv'), index=False, sep=';')        
    else:
        rel_df.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_relation_df'+ldb_aux+'.csv'), index=False, sep=';')

    return rel_df             
#%%
def task_function(tpp, adb_files, vtgs, year_act, BLUES, sc2_e, ldb = '', new_df = False, **kwargs): 
    """
    main function that set if the user wants to get data from adb or ldb
    the user should also define if the data will be generated from scratch 
    or will use an existing one previously created
    """
    task, path, project = tpp
    print('\nworking on '+task+' task')
    if task not in ['ADB', 'LDB']: #check if task is correctly set
        raise Exception("Task should be 'ADB' or 'LDB'")
    if task == 'ADB': #if its adb
        if new_df: #if it is new df
            adb_df = conv_data(tpp, adb_files, ldb, vtgs, year_act, BLUES, sc2_e, df_generator)# call funtion to get data from adb
        else: #ow
            try: #read csv that was created earlier
                adb_df = pd.read_csv(os.path.join(path,project,'IX',ldb,project+'_'+task+"_original_df.csv"), index_col=[0,1,2], header=[0,1], dtype='object', sep=';')# Works without sep
            except:
                raise Exception('There is no '+project+'_ADB_df.csv in '+ os.path.join(path,project )+'\n\ncheck your path/project or try new_df=True')   
        adb_df.sort_index(inplace=True) #sorting index
        adb_df = adb_df.fillna('nan')#replace nan value by 'nan'
        BLUES_tec_df = adb_df.copy() #making a copy
        return BLUES_tec_df #return adb df
    
    elif task == 'LDB': #if its ldb
        try: #read csv that was created earlier - this is the modified df for adb not the original
            adb_df = pd.read_csv(os.path.join(path,project,'IX',project+'_ADB_original_df.csv'), index_col=[0,1,2], header=[0,1], dtype='object',sep=';')#Works without sep
        except:
            raise Exception('There is no '+project+'_ADB_df.csv in '+ os.path.join(path,project ) +"\n\ncheck your path/project or try task = 'ADB'")        
        adb_df.sort_index(inplace=True) #sorting index
        adb_df = adb_df.fillna('nan')
        if new_df:     #if it is new ldb df
            ldb_df = conv_data(tpp, adb_files, ldb, vtgs, year_act, BLUES, sc2_e, df_generator)# call function to ldb data
        else: #ow
            try: #try to read previously created ldb df
                ldb_df = pd.read_csv(os.path.join(path,project,'IX',ldb,project+'_'+task+'_'+ldb+"_original_df.csv"), index_col=[0,1,2], header=[0,1], dtype='object',sep=';')#, Works without sep
                #adb_ldb_df.columns = pd.MultiIndex.from_tuples(adb_ldb_df.columns.values
            except:
                raise Exception('There is no '+project+'_LDB_'+ldb+'_df.csv in '+ os.path.join(path,project,ldb ) +"\n\ncheck your path/project or try new_df=True")

        #read adb_ldb df      
  
        ldb_df.sort_index(inplace=True) #sorting index
        ldb_df = ldb_df.fillna('nan')    #replace nan value by 'nan'

        t_ldb = ldb_df.index.get_level_values(level=0).tolist() #ldb technologies
        
        t_adb = adb_df.index.get_level_values(level=0).tolist() #adb technologies
        
        new_ldb_tec = [t for t in t_ldb if t not in t_adb] #new ldb tecs
        
        t_ldb_old = [t for t in t_ldb if t not in new_ldb_tec] #ldb tecs that were already in adb tecs
        adb_updt_df = adb_df.loc[idx[t_ldb_old, :, :]].copy() #df with only ldb tecs that are not new

        new_ldb_df = ldb_df.loc[idx[new_ldb_tec, :, :]] #df with only new ldb        
        adb_ldb_df = ldb_df.loc[idx[t_adb, :, :]].copy() #df with only adb tecs
        

        for col in ldb_df.columns.tolist(): #updating columns according to the adb
            adb_ldb_df[col] = np.where(adb_ldb_df[col] == 'nan', adb_updt_df[col], adb_ldb_df[col])        

        adb_ldb_df = pd.concat([adb_ldb_df, new_ldb_df])
        adb_ldb_df.sort_index(inplace=True) #sorting index
        adb_ldb_df = adb_ldb_df.fillna('nan')#replace nan value by 'nan'           
        
        adb_df.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_adb_'+ldb+'_original_df'+'.csv'), sep=';') #exporting adb_ldb_df

        adb_ldb_df.sort_index(inplace=True) #sorting index
        adb_ldb_df = adb_ldb_df.fillna('nan')    #replace nan value by 'nan'        
        t_adb_ldb = adb_ldb_df.index.get_level_values(level=0).tolist()        
        
        cc = [tt for tt in t_ldb if tt not in t_adb_ldb] #CHECK IF ALL LDB TECS ARE IN ADB_LDB DF
        if len(cc) >0:
            print('There are tecs that were not added correctly, Check: '+str(cc))

        #f_tecs = [ t for t in t_adb_ldb if t not in t_ldb] #eliminating tecs that are not in LDB BUT THERE ARE IN ADB
        
        return adb_df, adb_ldb_df, ldb_df

#%%
def relation_values(tpp,rel_df, year_act, count_y, ldb = "", emission_values=False):
    """
    # =============================================================================
    # Getting all relation upper and lower parameters
    # =============================================================================
    """
    task, path, project = tpp
    if emission_values:
        print('\nworking on emission values')
    else:  
        print('\nworking on relation values')
    if ldb != "":
        ldb_aux = "_"+ldb
    else:
        ldb_aux = ldb
    relation_upper = pd.DataFrame(columns=['relation', 'node_rel','year_rel','value']) #creating dataframe
    relation_lower = pd.DataFrame(columns=['relation', 'node_rel','year_rel','value']) #creating dataframe
    
    #removing rows that have no relevant info    
    rel_df = rel_df[~( (rel_df['upper'] == str(['nan']) ) & (rel_df['lower'] == str(['nan']) ) & (rel_df['cost'] == str([['nan']]) ) )]
    
    rel_ = rel_df['relation'][~rel_df['relation'].str.endswith('a')].values.tolist() #getting all relations and removing all land relations that ends with a
    rel_ = sorted(list(set(rel_))) #sorting and setting list
    
    for r in rel_: #looping thorugh relations
        nod = rel_df[rel_df['relation'] == r]['node'] #getting node
        aux_upper = pd.DataFrame(columns=['relation', 'node_rel','year_rel','value']) #creating dataframe
        aux_lower = pd.DataFrame(columns=['relation', 'node_rel','year_rel','value']) #creating dataframe
        for n in nod: #lopping thorugh nodes
            
            vru = [float(rr.strip().strip("'")) if rr.strip().strip("'") != 'nan' else 'nan' for rr in rel_df[(rel_df['relation'] == r) & (rel_df['node'] == n)]['upper'].values.tolist()[0].strip("[]").split(',')] #getting upper values
            vrl = [float(rr.strip().strip("'")) if rr.strip().strip("'") != 'nan' else 'nan' for rr in rel_df[(rel_df['relation'] == r) & (rel_df['node'] == n)]['lower'].values.tolist()[0].strip("[]").split(',')] #getting lower values
            
            if 'nan' not in list(set(vru)) and len(list(set(vru))) == 1: #if upper values is not nan and len == 1
                vru = vru #setting upper value
                aux = expand_grid_name(['relation', 'node_rel','year_rel','value'], #df columns
                                       [r], [n], year_act, vru)  #generating dataframe
                aux_upper = pd.concat([aux_upper, aux], ignore_index=True) #concatenating dataframe
            elif 'nan' not in list(set(vru)) and len(list(set(vru))) > 1: #if upper value is not nan and len > 1
                vru = vru[count_y:] #setting upper value
                while len(vru) < len(year_act):  #check if their length match
                    vru = vru + [vru[-1]] #extending value to match number of years
                for v,y in zip(vru, year_act): #looping thorugh years and values
                    aux = expand_grid_name(['relation', 'node_rel','year_rel','value'],
                                           [r], [n], [y], [v]) #generating dataframe
                    aux_upper = pd.concat([aux_upper, aux], ignore_index=True) #concatenating dataframe
            
            elif 'nan' in list(set(vru)) and len(list(set(vru))) > 1: #check if it is allright
                print('check - there are nan and float values - '+r)
                    
            if 'nan' not in list(set(vrl)) and len(list(set(vrl))) == 1: #if lower values is not nan and len == 1
                vrl = vrl #setting lower value
                aux = expand_grid_name(['relation', 'node_rel','year_rel','value'],
                                       [r], [n], year_act, vrl) #generating dataframe
                aux_lower = pd.concat([aux_lower, aux], ignore_index=True) #concatenating dataframe
            elif 'nan' not in list(set(vrl)) and len(list(set(vrl))) > 1: #if lower value is not nan and len > 1
                vrl = vrl[count_y:] #setting lower value
                while len(vrl) < len(year_act): #check if their length match
                    vrl = vrl + [vrl[-1]] #extending value to match number of years
                for v,y in zip(vrl, year_act): #looping thorugh years and values
                    aux = expand_grid_name(['relation', 'node_rel','year_rel','value'], #df columns
                                           [r], [n], [y], [v]) #generating dataframe
                    aux_lower = pd.concat([aux_lower, aux], ignore_index=True) #concatenating dataframe
            elif 'nan' in list(set(vrl)) and len(list(set(vrl))) > 1: #check if it allright
                print('check - there are nan and float values - '+r)

        if list(set(rel_df[rel_df['relation'] == r]['relation_type']))[0] ==  'relationsc:':
            aux_upper['year_rel'] = 'cumulative' 
            aux_lower['year_rel'] = 'cumulative'         
        relation_upper = pd.concat([relation_upper, aux_upper], ignore_index=True)
        relation_lower = pd.concat([relation_lower, aux_lower], ignore_index=True)
        
    relation_upper.drop_duplicates(keep='first', inplace=True)
    relation_lower.drop_duplicates(keep='first', inplace=True)  
    relation_upper.reset_index(inplace=True, drop=True)
    relation_lower.reset_index(inplace=True, drop=True)       
    if emission_values:
        relation_upper.rename(columns=dict(zip(relation_upper.columns.tolist(), [ 'type_emission','node', 'type_year', 'value'])), inplace = True)
        #bound_emission	node | type_emission | type_tec | type_year --> adjust emission upper
        relation_upper.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_emission_bound'+ldb_aux+'_df'+'.csv'), index=False, sep=';') #exporting upper dataframe to csv
        return relation_upper
    else:    
        relation_upper.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_relation_upper'+ldb_aux+'_df'+'.csv'), index=False, sep=';') #exporting upper dataframe to csv
        relation_lower.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_relation_lower'+ldb_aux+'_df'+'.csv'), index=False, sep=';') #exporting lower dataframe to csv
        return relation_upper, relation_lower

#%%
def relation_cost_calculaion(tpp,rel_df, year_act, count_y, ldb = "", emission_values=False):
    """
    # =============================================================================
    # Getting all relation upper and lower parameters
    # =============================================================================
    """
    if emission_values:
        print('\nworking on emission cost')
    else:  
        print('\nworking on relation cost')
    task, path, project = tpp
    if ldb != "":
        ldb_aux = "_"+ldb
    else:
        ldb_aux = ldb
    relation_cost = pd.DataFrame(columns=['relation', 'node_rel','year_rel','value']) #creating dataframe

    
    rel_ = rel_df['relation'][~rel_df['relation'].str.endswith('a')].values.tolist() #getting all relations and removing all land relations that end with a
    rel_= sorted(list(set(rel_))) #sorting and setting list
    
    for r in rel_: #looping thorugh relations
        nod = rel_df[rel_df['relation'] == r]['node'].values.tolist() #getting node
        for n in nod: #lopping thorugh nodes           
            cost_v = [float(rr.strip().strip("'")) if rr.strip().strip("'") != 'nan' else 'nan' for rr in rel_df[(rel_df['relation'] == r) & (rel_df['node'] == n)]['cost'].values.tolist()[0].strip("[]").split(',')] #getting upper values
            cost_v = [f for f in cost_v if str(f) != 'nan']
            if len(cost_v) >0 : #if upper values is not nan and len == 1
                if len(cost_v) == 1:
                    vru = cost_v #setting upper value
                    aux = expand_grid_name(['relation', 'node_rel','year_rel','value'], #df columns
                                           [r], [n], year_act, vru)  #generating dataframe
                    relation_cost = pd.concat([relation_cost, aux], ignore_index=True) #concatenating dataframe
                else: #if upper value is not nan and len > 1
                    vru = cost_v
                    vru = vru[count_y:] #setting upper value
                    while len(vru) < len(year_act):  #check if their length match
                        vru = vru + [vru[-1]] #extending value to match number of years
                    for v,y in zip(vru, year_act): #looping thorugh years and values
                        aux = expand_grid_name(['relation', 'node_rel','year_rel','value'],
                                               [r], [n], [y], [v]) #generating dataframe
                        relation_cost = pd.concat([relation_cost, aux], ignore_index=True) #concatenating dataframe


    if emission_values:
        relation_cost.rename(columns=dict(zip(relation_cost.columns.tolist(), ['type_emission', 'node' , 'year_act', 'value'])), inplace = True)
        #tax_emission	node | type_emission | type_tec | type_year    
        relation_cost.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_emission_cost'+ldb_aux+'_df'+'.csv'), index=False, sep=';') #exporting upper dataframe to csv
        return relation_cost
    else:    
        relation_cost.to_csv(os.path.join(path,project,'IX',ldb, project+'_'+task+'_relation_cost'+ldb_aux+'_df'+'.csv'), index=False, sep=';') #exporting upper dataframe to csv
        return relation_cost    
#%%
def old_tecs(BLUES_tec_df, first_hist_new_cap, tecso, tpp, BLUES=True):
    """
    Get old technologies: The ones that are defined as existing and cannot be expanded
    """
    print('\nworking on old_tecs')
    task, path, project = tpp 
    #creating new tecs for historical new capacity before 2008
    tn_old = [] #new techs according to the year
    tno = [] #tno is listing techs that have any historical new capacity old earlier than 2008
    tn_old_out = [] #removing tec bc all hisc are earlier than 2008
    for tn in tecso: #looping through tec list
        hisc_y = BLUES_tec_df.loc[tn, ('historical_new_capacity_years')].values[0][0].strip('[]').split(',') #look for historical new capacity values
        hisc_aux = [int(yy) for yy in hisc_y if yy != 'nan'] #get only the non nan values
        hisc_y = [int(yy) for yy in hisc_y if float(yy) <= 2000] #first_hist_new_cap]  #get values lesser or equal 2005

        if BLUES:
            if 2005 in (hisc_y):  #if 2005 is among the values
                #print(tn)
                #print(hisc_y)
                hisc_y = [yy for yy in hisc_y if yy<2000 or yy>2004] #get all values that are NOT between 2000 and 2004, all these values will be aggregated into 2005 technologies
        if len(hisc_aux) >0: # if there are hisc aux values
            if hisc_aux == hisc_y: #if values of hisc aux and hisc y are equal
                if BLUES:
                    if 'existing' in tn.lower() or 'exist' in tn.lower(): #check if it is an existing technology
                        tn_old_out.append(tn) #if yes append to old tecs that will be removed
        if len(hisc_y) >0: #if there are hisc values
            for yy in hisc_y: #looping through them
                tn_old.append(tn+'_'+str(yy)) #append year in the end of the name
                tno.append(tn) #append original tec name

    if BLUES:            
        return tn_old, tn_old_out
    else:
        return tn_old
#%%
# =============================================================================
# #LAND USE TECHNOLOGIES ADJUSTMENT
# =============================================================================
def land_adj(BLUES_tec_df, tecs_1, land_rel_df, first_land_tec, main_d_period):
    """
    This functions gets the total land availabulity by node and land class from upper and lower value in tthe 
    """
    print('\nworking on land tec adjustments')
    land_rel_dic = {} #creating dictionary
    land_nodes = list(set(BLUES_tec_df.loc[first_land_tec].index.get_level_values('node').values.tolist())) #nodes of the first land tec
    for Reg in land_nodes: #get nodes of the first land tec
    #GET THE UPPER AND LOWER VALUE OF EACH LAND CLASS AND SUMMING THEM TO CHECK IF THEY MATCH THE FIRST LAND BDA
        tup = land_rel_df[(land_rel_df['node'] == Reg) & (land_rel_df['relation'].str.endswith('a'))]['upper'].values.tolist() #getting the upper value of year 'a' (a means 2010 in the old framework) land relation 
        ttup = [t.strip("[]").split(',')[0].strip("'") for t in tup] #getting only the first value
        #ttup = list(itertools.chain(*tup)) #flattening the list
        ttup = [float(at) for at in ttup]#converting to value
        land_rel_dic[Reg+'_upper'] = sum(ttup) #summing the upper value to get the total land availability in the node (sum of all land classes)
    
        tlw = land_rel_df[(land_rel_df['node'] == Reg) & (land_rel_df['relation'].str.endswith('a'))]['lower'].values.tolist() #getting the lower value of 'a' (a means 2010 in the old framework) land relation 
        ttlw = [t.strip("[]").split(',')[0].strip("'") for t in tlw] #getting only the first value
        ttlw = [float(at) for at in ttlw]    #converting to value 
        land_rel_dic[Reg+'_lower'] = sum(ttlw) #summing the lower value to get the total land availability in the node
        
        badup = float(BLUES_tec_df.loc[idx[first_land_tec, Reg, :], ('bound_activity_up')].values[0][0].strip('[]').strip("'")) #getting bda up of the fisrt land technology
        badlw = float(BLUES_tec_df.loc[idx[first_land_tec, Reg, :], ('bound_activity_lo')].values[0][0].strip('[]').strip("'")) #getting bda lo of the fisrt land technology  
        
        if badup - badlw < 2: #check if they are the same, they must be because it provides the total area available
            land_rel_dic[Reg+'_'+first_land_tec] = badup #getting the bda value to compare it with the sum previously done
        else: #if lower and upper value are not equal check it
            print('check_:  '+Reg)
            break
        
        land_classes = sorted((list(set([lt[-2:] for lt in land_rel_df.relation.values.tolist() if lt.endswith('a')])))) #get all land classes for 2010
        for lt in land_classes: #looping through the different classes of land tec (A,B, C...G)
        #GET THE UPPER AND LOWER VALUE OF EACH LAND CLASS
            tup = land_rel_df[(land_rel_df['node'] == Reg) & (land_rel_df['relation'].str.endswith(lt))]['upper'].values.tolist() #getting the upper value of the land class
            ttup = [t.strip("[]").split(',')[0].strip("'") for t in tup] #getting only the first value
            ttup = [float(at) for at in ttup] #converting to value
            land_rel_dic[Reg+'_'+lt[0]+'_upper'] = sum(ttup) #summing the upper value to get the total availability of each land class (A, B...F,G)
        
            tlw = land_rel_df[(land_rel_df['node'] == Reg) & (land_rel_df['relation'].str.endswith(lt))]['lower'].values.tolist() #getting the loweer value of the land class
            ttlw = [t.strip("[]").split(',')[0].strip("'") for t in tlw] #getting only the first value
            ttlw = [float(at) for at in ttlw]    #converting to value
            land_rel_dic[Reg+'_'+lt[0]+'_lower'] = sum(ttlw) #summing the lower value to get the total availability of each land class (A, B...F,G)
                
    
    #tecs_L1 = [s for s in tecs_1 if s.startswith('Land_') and 'Bal_Total' not in s]      #getting primary land level tecs
    tecs_L_conv = sorted(list(set([s for s in tecs_1 if s.startswith('Conv_') ])))  #getting conversion tecs
    
    BLUES_tec_df.replace('Prim_Land1', 'Prim_Land', inplace=True) #adjusting primary land level
    BLUES_tec_df.replace('Prim_Land2', 'Prim_Land', inplace=True) #adjusting primary land level
    BLUES_tec_df.replace('Sec_Land2', 'Sec_Land', inplace=True) #adjusting seconday land level
    BLUES_tec_df.replace('Sec_Land1', 'Sec_Land', inplace=True) #adjusting seconday land level

    # replacing and adjusting conv tecs input level of primary land tecs
    kilen = list(range(1,1+len(BLUES_tec_df['input_commodities'].columns.values.tolist()) ))
    kolen = list(range(1,1+len(BLUES_tec_df['output_commodities'].columns.values.tolist()) ))
    for k in kilen:
        BLUES_tec_df.loc[idx[:, :, :], idx['input_levels', str(k)]] = np.where(np.logical_and(BLUES_tec_df.loc[idx[:, :, :], idx['input_levels', str(k)]].str.startswith('Prim_Land'), 
                                                                                                    BLUES_tec_df.loc[idx[:, :, :], idx['input_commodities', str(k)]].str.startswith('land') ) ,
                                                                                              'Land_Resource',
                                                                                              BLUES_tec_df.loc[idx[:, :, :], idx['input_levels', str(k)]])
        
    for k in kolen:
        # replacing and adjusting conv tecs output level from secondary to primary
        BLUES_tec_df.loc[idx[tecs_L_conv, :, :], idx['output_levels', str(k)]] = np.where(BLUES_tec_df.loc[idx[tecs_L_conv, :, :], idx['output_levels', str(k)]].str.startswith('Sec_Land'), 
                                                                                              'Prim_Land',
                                                                                              BLUES_tec_df.loc[idx[tecs_L_conv, :, :], idx['output_levels', str(k)]])

    
    land_classes = sorted([l[0] for l in land_classes]) #GETTING CLASSES OF LAND
    land_aux = ["Land_Available_"+l for l in land_classes] #CREATING LAND GENERATION TECHNOLOGY FOR EACH CLASS
    for ll in land_aux: #looping through land classes

        tp = list(zip(cycle([ll]), land_nodes, cycle(list(range(1,2))))) #combining the name of the tech, node, and the mode (activity)
        s = pd.DataFrame(columns=BLUES_tec_df.columns, index=pd.MultiIndex.from_tuples(tp, names=['tec', 'node','mode'])) #creating series
            #BLUES_tec_df = BLUES_tec_df.append(pd.Series(name=tn))
        BLUES_tec_df = BLUES_tec_df.append(s) #adding series as row
        BLUES_tec_df = BLUES_tec_df.sort_index(level=0)        
        for tn in tp:#looping thorugh land nodes
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['node_in','1']] = tn[1]
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['node_out','1']] = tn[1]
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['time','1']] = 'year'            
            #add lft =5
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['lifetime','1']] = str([main_d_period]) #lifetime = 5
            # muf =1
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['minimum_utilization_factor','1']] = '[1]' #muf must be 1 to produce the correct amount of land every time step

            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['investment_cost','1']] = '[0.1]' #inv c
            #add cf = 1
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['capacity_factor','1']] = '[1]' #CF IS ONE TO MATC MUF
            #add bdi up
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['bound_total_capacity_up','1']] = str([land_rel_dic['_'.join([tn[1],tn[0].split('_')[-1], 'upper'])]+0.5]) #GET BDC FROM LAND CONSTRAINTS
            #bdi lo
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['bound_total_capacity_lo','1']] = str([max(0, land_rel_dic['_'.join([tn[1],tn[0].split('_')[-1], 'lower'])]-0.5)]) #GET BDC FROM LAND CONSTRAINTS
            #output
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['output_commodities','1']] = '_'.join([tn[0].split('_')[0], tn[0].split('_')[-1]]).lower() #GET COMMODITY Land_A, Land_B....Land_G
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['output_levels','1']] = 'Land_Resource' #LEVEL
            BLUES_tec_df.at[idx[tn[0], tn[1], tn[2]], idx['output_values','1']] = '[1]'#VALUE
                
    new_t = BLUES_tec_df.index.get_level_values(level=0).tolist()
    #remove land balance total -->
    tecs_2 = [t for t in new_t if t not in ["Land_Bal_Total", 'Deficit_Land']] #REMOVING BALANCE AND DEFICIT
    BLUES_tec_df = BLUES_tec_df.loc[idx[tecs_2, :, :]] # removing "Land_Bal_Total" technology
    BLUES_tec_df = BLUES_tec_df.sort_index() #SORTING INDEX
    
    #tecs = BLUES_tec_df.index.get_level_values(level=0).tolist() #dataframe tecs
    tec_land1 = [s for s in tecs_2 if s.startswith('Land_')] #get all primary land tecs it does not include the land generation tecs
    tec_land1_new = ["_".join([tt[0].upper() + tt[1:] if tt == 'dbl' else tt for tt in s.split("_")]) for s in tec_land1] #new name of primary land tec adjusting dbl to Dbl

    
    dct_land1 = dict(zip(tec_land1, tec_land1_new)) #dictionary mapping the relation between old technology names and the new ones
    map_level(BLUES_tec_df, dct_land1, level=0) #replacing old tec names for the new ones based on the dictionary   
    
    BLUES_tec_df = BLUES_tec_df.sort_index() #sorting df
    tecs_2 = BLUES_tec_df.index.get_level_values(level=0).tolist() #get tecs
    tecs_2 = sorted(list(set(tecs_2)))


    BLUES_tec_df = BLUES_tec_df.fillna('nan')
    return BLUES_tec_df, tec_land1_new, tecs_2, land_nodes
#%%
def land_tec_link(link_t, land_generation ):
    """
    returning only land technologies chain that will be used to build land ix pattern
    """
    print('\nworking on land tec link')
    #land_rel_df['relation'][land_rel_df['relation'].str.endswith('T')].values.tolist() #get relations that ends with T they are the TOTAL BALANCE constraints
    land_chain = {} # dictionary to track land tecs chain
    for lc in land_generation:
        nt = 3 # number of levels after the first land technology   
        land_chain_aux = {} #aux dictionary
        t = 0 # counting
        land_chain_aux[str(t)] = [lc] #first land technology as list in the aux dict
        #land_chain[str(t)] = [first_land_tec] #first land technology as list in the land chain dict
        while t <= nt:  #while t does not reach nt do the following
            land_chain_aux[str(t+1)] = [] # creating list for the auxiliary level t+1
            for lt1 in land_chain_aux[str(t)]: # looping through aux dict
                lt1b = "_".join([tt[0].upper() + tt[1:] if tt == 'dbl' else tt for tt in lt1.split("_")])
                lto = link_t[lt1]['1gen_output_for_tec'] # getting output tecs of tec lt1 link_t = {key:{'0get_input_from_tec': i, '1gen_output_for_tec':o} for key,i,o in zip(tec_list, ii, oo)} #link of technologies
                lti = link_t[lt1]['0get_input_from_tec'] # getting input tecs of tec lt1  
                land_chain[str(t)+'_'+lt1b] = {} # creating dict to add outp and inp tecs of tec lt1
                land_chain[str(t)+'_'+lt1b]['output'] = lto #adding output tecs
                land_chain[str(t)+'_'+lt1b]['input'] = lti #adding input tecs
                land_chain_aux[str(t+1)].append(lto) #adding output tecs in the aux for the next level
            land_chain_aux[str(t+1)] = list(itertools.chain(*land_chain_aux[str(t+1)])) #flattening the list
            t+=1 #add one to go to the next level
    return land_chain
#%%
def land_ix_adj(land_nodes, tec_land1_new, land_chain):
    """
    # =============================================================================
    # #ADJUSTING LAND TO THE BLUES PATTERN
    # =============================================================================
    """
    print('\nworking on land tec ix conversion')
    # Check types of land to keep the balance according to the defined methodology
    #tec_land2 = list(set(["_".join(t.split("_")[:-1]) for t in tec_land1_new])) #types of land #forest, savanna, 
    class_land_tec = sorted(list(set(["_".join(t.split("_")[1:]) for t in tec_land1_new if 'Available' not in t.split("_")[1:]]))) #types of land for each class (A, B, C...)
    cat_land_dict = {} #aux dictionary
    for t in class_land_tec: #looping through type_class_lands
        tt = '1_Land_'+t  #adjusting name to match land chain
        #t = t[0].upper()+t[1:] #adjusting name
        cat_land_dict[t] = list(set([tl for tl in land_chain[tt]['output'] if not tl.startswith(('Conv', 'Oi_', 'Agr_')) ]))[0] #getting tecs for each land chain output of the class_type_land
    
    #type land tec and the technologies responsible for convert from primary to secondary
    cat_land_tec_df = pd.DataFrame(list(cat_land_dict.items()), columns=['type_land_tec', 'technology']) #converting dict to df to link cat_land_tec and its technologies
    
    tec_land1_set = [t for t in sorted(list(set(tec_land1_new))) if 'Available' not in t.split("_")[1:]] #get unique values for primary land tec

    #primary technologies that are impacted by type_land_tec
    map_land_rel_df = pd.DataFrame(columns=['node', 'technology', 'type_land_tec'])     #aux df
    for n in land_nodes:#looping through land nodes
        #generating df
        aux = pd.DataFrame({'node':n,
                            'technology':tec_land1_set,
                            'type_land_tec':class_land_tec})
        map_land_rel_df = pd.concat([map_land_rel_df, aux], ignore_index=True) #concatenating
    
    sec_land_tec = cat_land_tec_df['technology'].values.tolist() #secondary land tecs
    land_tec = tec_land1_set + sec_land_tec #all land tecs
    
    return cat_land_tec_df, map_land_rel_df, sec_land_tec, land_tec, class_land_tec
#%%
def land_yr_tec_ajd(BLUES_tec_df, tecs, tec_L, year_act, count_y):
    """
    This function takes all parameter information from the different land technologies (Land and Conv)
    for the different years they are defined in the BLUES model and organize it according to the year of each technology
    """
    print("\nworking on Land/Conv year tec adjustments")
    # =============================================================================
    land_data = BLUES_tec_df.copy()    #copying BLUES data without comment and historical land constraint    
    land_data = land_data.loc[idx[tec_L, :,:]] #subsetting based onland and Conv tecs

    tec_la = sorted(list(set( land_data.index.get_level_values(level=0).tolist() ) )) #getting a list of all land conv tecs
    tec_la_aux = sorted( list( set( ["_".join(t.split("_")[:-1]) if hasNumbers(t[-2:]) else t for t in tec_la] )) ) #getting land conv tecs not considering year elements
    for lca1 in tec_la_aux: #looping thorugh each land conv tec not considering the years
        #print(lca1)
        tec_la_aux2 = [t for t in tec_la if t.startswith(lca1)] #get all year tecs
        for nn in list( set (land_data.loc[idx[tec_la_aux2, :,:]].index.get_level_values(level=1).tolist() )): #looping thorugh nodes
            land_data_aux = land_data.loc[idx[tec_la_aux2, nn,:]].copy() #subsetting by tec and node
            land_data_aux = land_data_aux.loc[:, ~(land_data_aux == land_data_aux.iloc[0]).all()] #removing columns with equal values        
            #land_data_aux[[col for col in land_data_aux if not land_data_aux[col].nunique()==1]] #removing columns with equal values
            for lca in tec_la_aux2:  #looping thorugh tecs
                #ADJUSTING  EMISSIONS
                if "emissions" in land_data_aux.columns and 'emission_values' in land_data_aux.columns: #CHECK IF BOTH EMISSIONS AND EMISSION VALUES ARE AMONG COLUMNS
                    org = land_data_aux.loc[idx[lca, nn, :]]['emission_values'].values.tolist() #get the original value for emission values
                    # sort a list based on other list order
                    original_order = sort_list(list(range(len(land_data_aux.loc[idx[lca, nn, :]]['emissions'].values.tolist()[0][0].strip('[]').split(',') ) )),
                                                  [v.strip() for v in land_data_aux.loc[idx[lca, nn, :]]['emissions'].values.tolist()[0][0].strip('[]').split(',')] )   
                     
                    # sorting emissions #replace old by sorted one
                    land_data_aux.loc[idx[lca, nn, :], idx['emissions', '1']] = str(sorted( [v.strip() for v in land_data_aux.loc[idx[lca, nn, :]]['emissions'].values.tolist()[0][0].strip('[]').split(',')]  ) )
                     
                    ev = land_data_aux.loc[idx[lca, nn, :]]['emission_values'].values.tolist()[0][0].split('], [') #adjusting values, ordering them to match emissions names order
                    ev = [[r_t_z("{:.11f}".format(float(ne) ) ) for ne in re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", e)] for e in ev] #get numbers
                    evs = [ev[i] for i in original_order]  #sorting acccording to the rule
                    land_data_aux.loc[idx[lca, nn, :], idx['emission_values', '1']] = str(evs) #replace values     
                    mod = land_data_aux.loc[idx[lca, nn, :]]['emission_values'].values.tolist() #modified values
                    
                    if mod != org: #check if modified is equal to original #if they are not equal probably the emissions are different
                        print("There is an issue in the different land technologies for each year")
                        print(lca+' - '+nn)
                        sys.exit() #breaking script 
                if "emissions" in land_data_aux.columns:
                    land_data_aux = land_data_aux.drop(columns=['emissions'], level=0) #removing EMISSION COLUMN
            #land_data_aux = land_data_aux.loc[:, ~(land_data_aux == land_data_aux.iloc[0]).all()] #removing EMISSION COLUMN
        
            for cc in land_data_aux.columns: #looping through remaining columns
                tec_la_aux3 = tec_la_aux2[:] #include first tec (tec for 2010) in the sum
                #GETTING VALUES AND ADJUSTING THEM
                EVA = land_data_aux.loc[idx[tec_la_aux2[0], nn, :], cc].values.tolist()[0].split('], [') #.strip('[]').split(',')])
                EVA = [[float(ne) for ne in re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", e)] for e in EVA] #GET NUMBERS
                for ena in range(len(EVA)): #LOOPING THROUGH EACH LIST OF VALUES
                    EVA[ena] = EVA[ena][count_y:] #adjusting according to first year
                    while len(EVA[ena]) < len(year_act):#ADJUSTING LENGTH OF LIST
                        EVA[ena] = EVA[ena] + [0] #ADD ZEROS UNTIL MATCH LEN YEAR ACT
                    if len(EVA[ena]) > len(year_act): #adjusting length if it is bigger than year act
                        EVA[ena] = EVA[ena][:len(year_act)]
                if len(EVA) ==1: #IF LENGTH IS EQUAL 1
                    EVA = EVA[0] #GET ONLY FIRST LIST
                new_v = np.zeros_like(EVA) #NEW V ACCORDING TO EVA
                
                for lcaa in tec_la_aux3:   #LOOPING THROUGH technologies
                    if hasNumbers(lcaa[-2:]): #check tec has number in the end
                        yrl = int('20'+lcaa[-2:]) #converting number into years
                    else: #if ther is no number it is 2010
                        yrl = 2010
                    #GETTING VALUES AND ADJUSTING THEM
                    EV = land_data_aux.loc[idx[lcaa, nn, :], cc].values.tolist()[0].split('], [') #.strip('[]').split(',')])
                    EV = [[float(ne) for ne in re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", e) if float(ne) !=0] for e in EV] #GET NUMBERS
                    
                    for en in range(len(EV)): #LOOPING THROUGH EACH LIST OF VALUES
                        EV[en] = EV[en] #adjusting according to first year
                        while len(EV[en]) < len(year_act[:year_act.index(yrl)+1]): #ADJUSTING LENGTH OF LIST
                            EV[en] = [0] + EV[en]  #ADD ZEROS UNTIL MATCH LEN of current tec year
                        while len(EV[en]) < len(year_act): #ADJUSTING LENGTH OF LIST
                            EV[en] =  EV[en] + [0] #ADD ZEROS UNTIL MATCH LEN YEAR ACT                        
                        if len(EV[en]) > len(year_act): #check if length is greater than number of years
                            EV[en] = EV[en][count_y:len(year_act)]
                    if len(EV) ==1:  #IF LENGTH IS EQUAL 1
                        EV = EV[0] #GET ONLY FIRST LIST
                    new_v = new_v + np.array(EV) #summing values into new array value
                
                new_v = new_v[count_y:]
                if cc == ('bound_new_capacity_up', '1'):
                    new_v = [0] + list(new_v)[:-1] #it is not allowed construction in the first year -->bc of differenece between Original and IX years
                    new_v = np.array(new_v)
                    #print(new_v)
                    
                new_v = [str(new_v.tolist())] #adjust format to match original format
                #land_data_aux.loc[idx[tec_la_aux2[0], nn, :], cc].values.tolist()
                land_data_aux.loc[idx[tec_la_aux2[0], nn, :], cc] = new_v #replace old value by new value
            
            #renaming bound_new_capacity 2 bound_total_capacity
            new_c = [c.replace('new', 'new') if 'new' in c else c for c in land_data_aux.columns.levels[0]] #new column value
            land_data_aux.rename(columns=dict(list(zip(land_data_aux.columns.levels[0].tolist(), new_c))), level=0, inplace=True)  #rename column

            check_col = ['input_commodities', 'input_levels', 'output_commodities', 'output_levels']
            if  len([i for i in check_col if i in land_data_aux.columns.levels[0]]) >0:
                print("\nThe following columns should not be in land data_aux: " + str(check_col)  )
                sys.exit() #breaking script
                
            #land_data_aux.loc[idx[tec_la_aux2[0], nn, :]]
            land_data.update(land_data_aux.loc[idx[tec_la_aux2[0], nn, :]]) #UPDATE COLUMNS  OF LAND_DATA AUX INTO LAND DATA DF FOR LAND/CONV TEC
    
    return land_data
#%%
def bal_eq(tpp, adb_files, ldb):
    print('\nworking on balance equalty commodity-levels')
    task, path, project = tpp

    #############################################################################################################
    bal_eq_list = [] #aux ldr list

    for Reg in adb_files: #looping through regions
        
        if task == 'ADB':
            dt = os.path.join(path, project, Reg,"data",Reg+'.adb') #path - energy forms
        elif task == 'LDB':
            dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+'.ldb') #path - energy forms
        
        dt1 = open(dt) #reading adb
        
        E_form_block = "" #creating block
        found = False #starting
        for line in dt1: #looping through lines
            if found: #when found is True add line
                if line.strip() == 'demand:':  #stop line
                    break
                E_form_block += line.strip() + "\n"  #addgin line
            else: #found is false
                if line.strip() == 'energyforms:': #starting line
                    found = True #set found = True

        dt1.close() #close file

        E_forms1 = E_form_block[:] #copying E form
        E_forms = E_forms1.split('*') #splitting each energy form
        E_forms = [e for e in E_forms if e != '' if e != '\n'] #removing empty lines 
                   
        for n in E_forms: #looping through energy forms
            f1 = ''.join(n.split('\t')).split('\n') #splitting by tab, then joining then splitting by line
            f1 = [t for t in f1 if t!=''] #excluding empty lines
            LVL = f1[0].split(' ') #getting the level
            for com in f1[1:]: #looping through levels
                com = com.strip(' ') #removing empty characteres in the beggining and in the end of str
                if not com.startswith('#'): #if its not a comment
                    COM = com.split(' ') #split by empty characteres
                    COM = [c for c in COM if c!=""]                    
                    if len(COM) >2 and COM[-1] == 'f': #get commodities that are listed as load factor
                        bal_eq_list.append((LVL[0],COM[0].lower())) #append ldr commodities into ldr list
    
    bal_eq_df = pd.DataFrame(bal_eq_list, columns=['level', 'commodity',])
    return bal_eq_df
#%%
def lvl_com_check(tpp, adb_files, ldb):
    """
    Check level commodity values
    """
    print('\nworking on level-commodity check')
    task, path, project = tpp

    if task == 'ADB':
        lvl_com_cm_df = pd.DataFrame(columns=['level', 'commodity','code',  'comment', 'status'])
    elif task == 'LDB':
        lvl_com_cm_df = pd.read_csv(os.path.join(path, project,'IX', project+'_ADB_level_commodity'+'.csv'), sep=';')  
        #a = pd.read_csv(os.path.join(path, project, project+'_ADB_level_commodity'+'.csv'), sep=';')  
        #all(a==lvl_com_cm_df)
#E_form is aggregating all required lines --> lines between 'energyforms' and 'demand'
    for Reg in adb_files:#looping through regions    
        if task == 'ADB':
            dt = os.path.join(path, project, Reg,"data",Reg+'.adb') #path - energy forms
        elif task == 'LDB':
            dt = os.path.join(path, project, Reg,"data",Reg+'_'+ldb+'.ldb') #path - energy forms

        dt2 = open(dt)#open file
        
        E_form_block = "" #creating block
        found = False #set non starting value
        for line in dt2: #looping through lines
            if found: #if start
                if line.strip() == 'demand:':  #stop condition
                    break
                E_form_block += line.strip() + "\n"  #add line
            else: #if found is false
                if line.strip() == 'energyforms:': #starting condition
                    found = True #set as found true

        dt2.close() #close file
        
        Level_form_dict = {}
        E_forms1 = E_form_block[:] #str type has no copy
        E_forms = E_forms1.split('*') #splitting each energy form
        E_forms = [e for e in E_forms if e != '' if e != '\n'] #excluding lines with no info
        lcc = []        
        for n in E_forms: #looping through energy forms
            f1 = ' '.join(n.split('\t')).split('\n') #splitting by tab, then joining then splitting by line
            #fa = ' '.join(n.split('\t')).split('\n')
            f1 = [t.strip() for t in f1 if t!=''] #excluding empty lines
            LVL = f1[0].split(' ') #getting the level
            for com in f1[1:]:  #looping thorugh commoditties
                com = com.strip() #removing empty characteres in begining and in the end of str
                if not com.startswith('#'): #if it is not comment
                    c = com.split(' ') #split by space
                    c = [c for c in c if c!=""]
                    cm_ = 'nan'
                    if f1.index(com)+1 != len(f1):
                        if f1[f1.index(com)+1].startswith('#'):
                              cm_ =  f1[f1.index(com)+1]
                    if len(c) > 2:
                        COM = c[-2]
                    else:
                        COM = c[-1]
                    Level_form_dict[COM+'-'+LVL[-1]] = LVL[0]+'|'+c[0].lower() #building dictionary level|commodity
                    
                    tp_lcc = (LVL[0], c[0].lower(),COM[-1]+'-'+LVL[-1], cm_)
                    lcc.append(tp_lcc)
                lcc = unique(lcc)
                
        lcc_df = pd.DataFrame(lcc, columns=['level', 'commodity', 'code', 'comment'])
        
        dt2 = open(dt)#open file
        
        E_form_block = "" #creating block
        found = False #set non starting value
        for line in dt2: #looping through lines
            if found: #if start
                if line.strip() == 'endata':  #stop condition
                    break
                E_form_block += line #add line
            else: #if found is false
                if line.strip() == 'demand:': #starting condition
                    found = True #set as found true
    
        dt2.close() #close file        

        lcc_df['status'] = ['Used' if any(x in E_form_block for x in  [code+' ', code+'\n', code+'\t', code+'-'])  else 'Not_Used' for code in lcc_df['code'].values.tolist()]

        
        lvl_com_cm_df = pd.concat([lvl_com_cm_df, lcc_df], ignore_index=True)
        
        lvl_com_cm_df['status_v'] = np.where(lvl_com_cm_df['status']=='Used', 99, 0)
        
        lvl_com_cm_df.sort_values(['level', 'commodity', 'status_v'], inplace=True)
        lvl_com_cm_df.drop_duplicates(['level', 'commodity'],  keep='last', inplace=True, ignore_index=True)    
        lvl_com_cm_df = lvl_com_cm_df[['level', 'commodity','code',  'comment', 'status']]
        
        
    if task == 'LDB': 
        lvl_com_cm_df.to_csv(os.path.join(path, project, 'IX', ldb, project+'_'+task+'_level_commodity'+'_'+ldb+'.csv'), sep=';', index=False)
    else:
        lvl_com_cm_df.to_csv(os.path.join(path, project, 'IX', ldb, project+'_'+task+'_level_commodity'+'.csv'), sep=';', index=False)  
    
    return lvl_com_cm_df

