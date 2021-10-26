"""
Created on Thu Aug 15 16:55:01 2019

@author: Fabio Diuana
"""
"""
Writing all dataframes for the technologies in the model
"""
from datetime import datetime
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt 
from itertools import cycle
import itertools
import string
import pickle
import os
import sys
import statistics as stat
idx = pd.IndexSlice #pandas option to slice list
from BLUES_conversion_r170_functions import * #import from BLUES conversion functios
# =============================================================================
# #functions that I am using in the main function
# sys.path.append('D:FabioDiuana_Pdrive_IIASAModelsNEST-BLUESInput_data_scripts')
# sys.path.append('D:/Fabio/Diuana_Pdrive_IIASA/Models/NEST-BLUES/Input_data_scripts')
# from xfuncs import *
# =============================================================================
#%%
def build_adb_dfs(tpp, ldb, BLUES_tec_df, tecs_, all_time_info, old_tec_list, ldr, ldr_info, land_rel_df, rel_df, BLUES=True, tec_dict = {}):
    startTime = datetime.now() #getting starting time
    warning_list = [] #listing issues
    not_def_rel = [] #listing tecs that have no unit group defined
    final_ldr, first_ldr, int_ldr, ldr_tecsv, ldr_list, ldr_original = ldr #LDR_TECSV comprises all ldr tecs
    years, count_y1, time, first_d_period, main_d_period, start_year = all_time_info #, first_hist_new_cap 
    task, path, project = tpp 
    tn_old, tn_old_out, tn_old_out_y = old_tec_list  # it is not being used anymore
        
    #for tn1 in tecs_: #looping through technologies
    def aux__build_tec(tn1):

        #############print(tn1)
        #if hasNumbers(tn1[-4:]): #checking it all the last 4 characteres are number --> if it is True it means that they are existing technologies installed before 2008
        #    tn = tn1[:-5] #removing the last 4 characteres to get the parameters
        #else:
        tn = tn1 #ow it is not necessary
            
        a = BLUES_tec_df[BLUES_tec_df.index.get_level_values(0)==tn] #getting all parameters of the current tec; subesetting the maind dataframe
        
        #add all empty parameters dataframe:
        vtg_year_df = pd.DataFrame(columns=['node',  'vintage', 'year_act'])
        vtg_year_time_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'time'])
        lft_df = pd.DataFrame(columns=['node',  'vintage', 'value'])
        ctime_df = pd.DataFrame(columns=['node',  'vintage', 'value'])
        inpt_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode', 'node_in', 'commodity', 'level', 'time', 'time_in', 'value'])
        output_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode', 'node_out', 'commodity', 'level', 'time', 'time_out', 'value'])
        main_output_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode', 'time', 'value'])
        main_input_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode', 'time', 'value'])
        capacity_factor_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'time', 'value'])
        var_cost_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode', 'time', 'value'])
        inv_cost_df = pd.DataFrame(columns=['node',  'vintage', 'value'])
        fix_cost_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'value'])
        muf_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'value'])
        bdi_lo_df = pd.DataFrame(columns=['node',  'year_act', 'value'])
        bdi_up_df = pd.DataFrame(columns=['node',  'year_act', 'value'])
        bdc_lo_df = pd.DataFrame(columns=['node',  'vintage',  'value'])
        bdc_up_df = pd.DataFrame(columns=['node',  'vintage',  'value'])
        bda_lo_df = pd.DataFrame(columns=['node',  'year_act', 'mode', 'time',  'value'])
        bda_up_df = pd.DataFrame(columns=['node',  'year_act', 'mode', 'time', 'value'])
        emission_factor_df = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode',  'emission', 'value'])
        historical_new_capacity_df = pd.DataFrame(columns=['node',  'vintage',  'value'])
        relation_activity_df = pd.DataFrame(columns=['relation','node_rel','year_rel','node_loc','year_act','mode', 'value']) #creating dataframe
        relation_capacity_df = pd.DataFrame(columns=['relation','node_rel','year_rel', 'value']) #creating dataframe
          
        #NODE 
        nodes = list(set(BLUES_tec_df.loc[tn].index.get_level_values('node').values.tolist())) #getting nodes of the technology
        
        for nod in nodes:    #looping through nodes
            
            print(tn1 + ' - '+nod )
            vtgs = years[:] #setting vintages
            year_act = years[1:][:] #setting modelled years
            count_y = count_y1
            #MODE
            mod = BLUES_tec_df.loc[idx[tn, [nod], :]].index.get_level_values('mode').values.tolist() #getting modes of the technology
            
            #CHECK IF HISTORICAL LAND OR HISC ARE DIFFERENT FROM NAN
            if BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_land_constraint')].values.tolist()[0][0] != 'nan' or BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_new_capacity_years')].values.tolist()[0][0] != 'nan': #if land constraint is not nan
                #CHECK IF HISC IS DIFFERENT FROM NAN
                if BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_new_capacity_years')].values.tolist()[0][0] != 'nan':
                    hisc_b1 = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_new_capacity')].values.tolist())] #getting historical values
                    hisc_b1 = [float(ttv.strip().strip("'")) for ttv in hisc_b1[0][0].strip('[]').split(',') if ttv != 'nan'] #removing nan and get values as float
                    if sum(hisc_b1) > 0: #IF HISC SUM GREATER THAN 0
                        vtgs = years[:] #setting vintages  
                    else:
                        #vtgs = year_act[:]
                        vtgs = years[:]
                else: #IF HISC IS NAN IT MEANS THERE HISTORICAL LAND
                    vtgs = years[:]
                    
# =============================================================================
#             elif BLUES_tec_df.loc[idx[tn, [nod], :], ('investment_cost')].values.tolist()[0][0] != 'nan':
#                 vtgs = year_act[:]
# =============================================================================
                
            else: # IT DOES NOT NEED FIRST VTG
                vtgs = year_act[:]
    ##########################################################################################################################      
            #LIFETIME: #check lyear and adjust it
            lft_aux = pd.DataFrame(columns=['node',  'vintage',  'value']) #creating aux dataframe for lifetime
            lft_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('lifetime')].values.tolist())] #getting a list of values for all modes and removing the equal ones
            if len(lft_b) == 1: # if it is TRUE it means that the values are the same for all modes
                if str(lft_b[0][0]) == 'nan': #if the value is 'nan'                      
                        #lft_a = 'nan'
                        lft_a = [25]
                        
                        
                else: # if its not nan
                    lft_a = [float_round(float(v.strip().strip("'")),0, ceil ) for v in lft_b[0][0].strip('[]').split(',')] #converting string list into list of values
                    #lft_a = [lfv if lfv>=1 else 1 for lfv in lft_a] #replacing lifetime smaller than 1

                if lft_a != 'nan':                         
                    if len(lft_a) == 1: #it means the values are the same for the entire period
                        aux = expand_grid_name(['node', 'vintage',  'value'],   #datafamere collumns
                                                           [nod], vtgs,  lft_a ) #generating dataframe  from these values
                        aux = aux[lft_aux.columns] #adjusting to the correct collums order
                        lft_aux = pd.concat([lft_aux, aux], ignore_index=True) #concatenating aux dataframe into empty df
                    else:           # it means the values vary according to the years 
                        lft_a = [float_round(float(sv), 0, ceil) for sv in lft_a if str(sv) !='nan'] #converting sv that is not nan into float
                        # check if I need to repeat the values for the 2 last years or if it is not considering the 2 first years
                        while len(lft_a) < len(vtgs): #making the quantity of lfta a values equal to the number of vintage years
                            lft_a = lft_a + [lft_a[-1]] #always adding the last value to make them equal
                        
                        for vin,v in zip(vtgs, lft_a): #looping through vintages and lft values
                            aux = expand_grid_name(['node', 'vintage',  'value'],  #datafamere collumns
                                                           [nod], [vin], [v] ) #generating dataframe  from these values
                            aux = aux[lft_aux.columns] #adjusting to the correct collums order
                            lft_aux = pd.concat([lft_aux, aux], ignore_index=True) #concatenating aux dataframe into empty df
            
            else: # the values are different for each mode           
                #print("CHECK IT bc for lft values should be the same for all modes")
                warning_list.append("CHECK IT bc for lft values should be the same for all modes:  "+str(tn)+' - '+str(nod)) # it should not happen, add info into list of warning
            
            lft_aux = lft_aux[lft_df.columns]    #adjusting to the correct collums order

            if lft_a !='nan':
                lft_v = lft_a[0] #getting first value as the reference value
            else:
                lft_v = 1 #first_d_period
                
            if BLUES_tec_df.loc[idx[tn, [nod], :], ('last_year')].values.tolist()[0][0] != 'nan': # if last year is not nan , it means the tec should be decomissioned in this year
                lft_aux['value'] = float(BLUES_tec_df.loc[idx[tn, [nod], :], ('last_year')].values.tolist()[0][0].strip('[]')) - vtgs[0] #adjusting the lifetime value to its match with the last year
                lft_v = [float( BLUES_tec_df.loc[idx[tn, [nod], :], ('last_year')].values.tolist()[0][0].strip('[]') ) - v for v in vtgs]#[0] #doing the same for the value taht will be used to build the auxiliary time dataframe
                lft_v = [lf for lf in lft_v if lf>0]
                lv = float(BLUES_tec_df.loc[idx[tn, [nod], :], ('last_year')].values.tolist()[0][0].strip('[]'))
                if lv not in years:
                    lv = years[(min([int(abs(lv-y)) for y in years]))]
                    
                year_act = year_act[: year_act.index(lv)+1] #ADJUSTING YEAR ALL OF THIS TECHNOLOGY in order to IT DOES NOT GO BEYOND THE LAST YEAR
                #vtgs = vtgs[:vtgs.index(lv)+1]
                vtgs = vtgs[:len(lft_v)]
    
            if len([lft_a]) >1: # if lfta contains more than one list of values it is wrong
                warning_list.append('lfta has more than one list of values:  '+str(tn)+' - '+str(nod))


            if BLUES_tec_df.loc[idx[tn, [nod], :], ('first_year')].values.tolist()[0][0] != 'nan': #check if first year is not nan, in this case the tec cannnot be comissioned before this year
                fv = float(BLUES_tec_df.loc[idx[tn, [nod], :], ('first_year')].values.tolist()[0][0].strip('[]')) #getting the value, from string to value
                #if fv not in vtgs: # checking if first year is among vintage year options
                fv = [v for v in vtgs if v>fv][0] # in case not , it is taking the first year later than the first year (this approach matches Original BLUES)
                   
                vtgs = vtgs[vtgs.index(fv):] #adjusting vintage period according to the first year for this tec
                year_act = year_act[year_act.index(fv):] #adjusting year all period according to the first year for this tec
                count_y = years.index(fv) -1 #adjusting the count_y value         
                
                #vtg_year = vtg_year[vtg_year.vintage >= fv]
                #vtg_year_time = vtg_year_time[vtg_year_time.vintage >= fv]
                #year_time = year_time[year_time.year_act.isin(vtg_year.year_act.values.tolist())]
                #vtg_time = vtg_time[vtg_time.vintage >= fv]     

            vtg_year_aux = expand_grid_name(['vintage','year_act'], vtgs,year_act) #generating df with all vintages and years
            vtg_year_aux = vtg_year_aux[vtg_year_aux.vintage <= vtg_year_aux.year_act] #subseting for years>vintages
            vtg_year_aux = vtg_year_aux.reset_index().drop(columns='index') #removing index
            vtg_year_aux['vintage_aux'] = [start_year[y] for y in vtg_year_aux['vintage']]
            vtg_year_aux['year_act_aux'] = [start_year[y] for y in vtg_year_aux['year_act']]

            lft_aux = lft_aux[lft_aux['vintage'].isin(vtgs)]# GETTING TECHNICALL LIFETIME DF SET BY VTGS
            
            if isinstance(lft_v, list):#len(lft_v) >1: #ADJUSTING LIFETIME VALUE ACCORDING TO DIFFERENT VTGS WHEN LIFTEIME VARIES OVER TIME
                if len(lft_v) >1:
                    conditions = [lft_aux['vintage'] ==  v for v in vtgs]
                    choices = lft_v
                    lft_aux['value'] = np.select(conditions, choices, default="ERRO")
                    lft_aux['value'] = lft_aux['value'].astype(float)
                    vtg_year = vtg_year_aux.copy()
            else:
                #generating df that correlates vintages and years according to the lifetime # by removing values that are lesser than lifetime
                vtg_year = vtg_year_aux.drop(vtg_year_aux[np.logical_or(vtg_year_aux.year_act_aux-vtg_year_aux.vintage_aux >= lft_v, # it was >= earlier
                                                                vtg_year_aux.year_act-vtg_year_aux.vintage < 0)].index).reset_index().drop(columns='index')

            vtg_year = vtg_year[['vintage','year_act']]
            if BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_land_constraint')].values.tolist()[0][0] != 'nan' or BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_new_capacity_years')].values.tolist()[0][0] != 'nan': #if land constraint is not nan
                if lft_v == first_d_period:
                    aux_vtg_year = pd.DataFrame({'vintage':[vtgs[0]],
                                                 'year_act':[year_act[0]] })
    
                    vtg_year = pd.concat([aux_vtg_year, vtg_year], ignore_index=True)
                
            
            #generating df that correlates vintages, years and time according to the lifetime
            vtg_year_time = expand_grid_name(vtg_year.columns.tolist()+['time'],
                                        #Line above is to get the column names and line below is to get the values
                                         vtg_year.values, time)

            lft_df = pd.concat([lft_df, lft_aux], ignore_index = True) #concatenating aux dataframe into main empty df 
              
            vtgs = sorted(list(set(vtg_year['vintage'].values.tolist() ) ) )
            #year_time = vtg_year_time[vtg_year_time.vintage==vtg_year_time.year_act].copy()[['year_act', 'time']] #generating year time df according to the lifetime
            #vtg_time = vtg_year_time[vtg_year_time.vintage==vtg_year_time.year_act].copy()[['vintage', 'time']] #generating vintage time df according to the lifetime
            year_time = vtg_year_time[['year_act', 'time']].copy() #generating year time df according to the lifetime
            year_time.drop_duplicates(keep='first', inplace=True, ignore_index=True)
            vtg_time = vtg_year_time[['vintage', 'time']].copy() #generating vintage time df according to the lifetime
            vtg_time.drop_duplicates(keep='first', inplace=True, ignore_index=True)
            
            year_act = sorted(list(set(vtg_year.year_act.values.tolist()[:])))
            if tn not in ldr_tecsv: #if technology is not a ldr tec
                vtg_year_time['time'] = 'year' #getting year as time
                vtg_year_time.drop_duplicates(keep = 'first', inplace=True) #new version has ignore_index|| #dropping duplicates since there will be 12 year time
                vtg_year_time.reset_index(drop=True, inplace=True)#reseting and dropping index
                
                #doing the same for other time auxiliary df
                year_time['time'] = 'year'
                year_time.drop_duplicates(keep = 'first', inplace=True) #new version has ignore_index
                year_time.reset_index(drop=True, inplace=True)
                
                #doing the same for other time auxiliary df
                vtg_time['time'] = 'year'
                vtg_time.drop_duplicates(keep = 'first', inplace=True) #new version has ignore_index
                vtg_time.reset_index(drop=True, inplace=True)
                
            vtg_year2 = vtg_year.copy()    
            vtg_year_time2 = vtg_year_time.copy()
            vtg_year2['node'] = nod
            vtg_year_time2['node'] = nod
            vtg_year2 = vtg_year2[['node', 'vintage', 'year_act']]
            vtg_year_time2 = vtg_year_time2[['node', 'vintage', 'year_act', 'time']]
            vtg_year_df = pd.concat([vtg_year_df, vtg_year2], ignore_index=True)
            vtg_year_time_df = pd.concat([vtg_year_time_df, vtg_year_time2], ignore_index=True)
    ##########################################################################################################################  
    # =============================================================================
    #     #CONSTRUCTION TIME
    # =============================================================================
            
            ctime_aux = pd.DataFrame(columns=['node',  'vintage',  'value']) #creating aux dataframe for ctime
            ctime_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('construction_time')].values.tolist())] #getting a list of values for all modes and removing the equal ones
            if len(ctime_b) == 1: # if it is TRUE it means that the values are the same for all modes
                if str(ctime_b[0][0]) != 'nan': #if its nan
                    ctime_a = [float(v.strip().strip("'")) for v in ctime_b[0][0].strip('[]').split(',')] #coverting string list ('[a,b,c...z]') into list of float values
                    
                    if len(ctime_a) == 1: #it means the values are the same for the entire period
                        aux = expand_grid_name(['node', 'vintage',  'value'],  #df columns
                                                           [nod], vtgs,  ctime_a ) #df values
                        aux = aux[ctime_aux.columns] #oredring columns
                        ctime_aux = pd.concat([ctime_aux, aux], ignore_index=True) #concatenating into aux ctime df
                    else:           # it means the values vary according to the years 
                        ctime_a = [float(sv.strip("'")) for sv in ctime_a if str(sv) !='nan'] #removing nan values
                        # check if I need to repeat the values for the 2 last years or if it is not considering the 2 first years
                        while len(ctime_a) < len(vtgs): #checking if vintages and values length match
                            ctime_a = ctime_a + [ctime_a[-1]] #making them match
                        
                        for vin,v in zip(vtgs, ctime_a): #looping through vintages and values
                            aux = expand_grid_name(['node', 'vintage',  'value'],  #df columns
                                                           [nod], [vin], [v] ) #df values
                            aux = aux[ctime_aux.columns] #ordering columns
                            ctime_aux = pd.concat([ctime_aux, aux], ignore_index=True) #concatenating aux df for each vintage-value pair
            
            else: # the values are different for each mode           
                #print("CHECK IT bc for ctime values should be the same for all modes")
                warning_list.append("CHECK IT bc for ctime values should be the same for all modes:  "+str(tn)+' - '+str(nod)) #it should not happen add into warning list
                
            ctime_aux = ctime_aux[ctime_df.columns]    #oredering columns
            ctime_df = pd.concat([ctime_df, ctime_aux], ignore_index = True) #concatenating into final df
    ##########################################################################################################################      
    # =============================================================================
    #     #MAIN INPUT VALUES
    # =============================================================================
            main_in_val_aux = BLUES_tec_df.loc[idx[tn, [nod], :], idx['input_values', '1']].values.tolist() #getting main output value
            main_in_val1 = [[vv.strip('[]')] for vv in main_in_val_aux if vv.strip('[]') != 'nan' ] #adjust it
            if 'nan' not in main_in_val1:
                main_in_val2 = [[float(vv.strip("'")) for vv in v[0].split(', ')] for v in main_in_val1] #converting to value
                main_in_val = [] #IT IS USED TO ADJUST OTHER VALUES
                mmi_mode = 0 
                for main_in_val_adj in main_in_val2:
                    mmi_mode+=1
                    #mmi_mode = main_in_val2.index(main_in_val_adj) + 1 #getting mode
                    if len(main_in_val_adj) > count_y: #if there is more than one value
                        main_in_val_adj = main_in_val_adj[count_y:] #adjust based on count y
                    while len(main_in_val_adj) < len(year_act): #make length of years and values match
                        main_in_val_adj = main_in_val_adj + [main_in_val_adj[-1]] #repeating the last one
                        
                    #main_in_val = main_in_val_adj[:len(year_act)]
                    main_in_val.append(main_in_val_adj[:len(year_act)]) #append into main_val
                    
                    for mvalmi, y in zip(main_in_val_adj, year_act): #generating main val df
                        main_in_val_df = expand_grid_name(['node', 'year_act','mode', 'value'], #dataframe columns
                                                                              [nod], [y], [mmi_mode], [mvalmi]).merge(      #dataframe values
                                                                                                              vtg_year_time, how='left', on='year_act')   #merging with time auxiliary df 
                        main_in_val_df = main_in_val_df[main_input_df.columns]
                        main_input_df = pd.concat([main_input_df, main_in_val_df], ignore_index=True)    
    ##########################################################################################################################  
    # =============================================================================
    #     #INPUT
    # =============================================================================
            com_aux = BLUES_tec_df.loc[idx[tn, [nod], :], 'input_commodities'].values.tolist() #getting all commoditity values
            lvl_aux = BLUES_tec_df.loc[idx[tn, [nod], :], 'input_levels'].values.tolist() #getting all level values
            com = [[c for c in cc if c!= 'nan'] for cc in com_aux] #removing nan values and keeping each mode in a sublist
            lvl = [[l for l in ll if l  != 'nan'] for ll in lvl_aux] #removing nan values and keeping each mode in a sublist
    
            if len(com) != len(lvl): #checking if number of commoditties and levels match in the mode aspect
                warning_list.append('OUT_check why len com is different from len lvl:  '+str(tn)+' - '+str(nod)) # if not adding warn into warning list
                #print('OUT_check why len com is different from len lvl:\t'+tn)
                #break
            else:  #go ahead
                val_aux = BLUES_tec_df.loc[idx[tn, [nod], :], 'input_values'].values.tolist()
                val = [[v.strip('[]') for v in vv if v != 'nan'] for vv in val_aux]
                if len(val) != len(lvl): #checking if number of commoditties and levels match in the mode aspect
                    warning_list.append('OUT_check why len val is different from len lvl/com:  '+str(tn)+' - '+str(nod)) # if not adding warn into warning list
                    #print('OUT_check why len val is different from len lvl/com:\t'+tn)
                else: #go ahead
                    aux_i = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode',  'commodity', 'level', 'time', 'value']) #auxiliary empty dataframe
        # =================================================================================================
                    for m,c1,l1,mv1 in zip(mod, com, lvl, val): # looping through mode, commodity, level and values
                        #aux_m = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode',  'commodity', 'level', 'time', 'value']) #auxiliary empty dataframe
                        if len(mv1) != len(l1) != len(c1): # check if the number of values, levels and commoditties in each list matches  in the mode level
                            warning_list.append('OUT_check len of each list from val, com and lvl:  '+str(tn)+' - '+str(nod)) #if it is not -- adding info into warning list
                            #print('OUT_check len of each list from val, com and lvl:\t'+tn)
                        else: #if all values length match
                            for c, l, mv in zip(c1,l1,mv1): #looping through each values of commodity, levl and value i the mode subset
                                #aux_1 = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode',  'commodity', 'level', 'time', 'value']) #auxiliary empty dataframe
                                mval1 = [float(v.strip().strip("'")) for v in mv.strip('[]').split(',')] #converting value or list of vlues to float or list of float
                                if len(mval1) == 1: #it means the values does not vary over yars
                                    aux_B_i = expand_grid_name(['node', 'vintage','mode', 'commodity', 'level', 'value'], #dataframe columns
                                                                 [nod], vtgs, [m], [c], [l], mval1).merge(              #dataframe values
                                                                                                 vtg_year_time, how='left', on='vintage') #merging with time auxiliary df
                                    aux_B_i = aux_B_i[aux_i.columns] #adjusting order of columns
                                    aux_i = pd.concat([aux_i, aux_B_i], ignore_index=True) #concatenating into the first auxiliary df
                                else: # in this case the values vary according to the years
                                    if len(mval1) > count_y:
                                        mval1 = mval1[count_y:] #removing value according to the first year_act
                                    while len(mval1) < len(year_act): # adjustin the length of the list values
                                        mval1 = mval1 + [mval1[-1]] # making it matches the length of modelled years
                                    for mval, y in zip(mval1, year_act): #looping thorugh years and values
                                        aux_B_i = expand_grid_name(['node', 'year_act','mode', 'commodity', 'level', 'value'], #dataframe columns
                                                                              [nod], [y], [m], [c], [l], [mval]).merge(      #dataframe values
                                                                                                              vtg_year_time, how='left', on='year_act')   #merging with time auxiliary df                               
    
    # =============================================================================
    #                                     conditions =  [ (aux_B['vintage'] == vtgs[ky]) for ky in list(range(mvlen))]
    #                                     choices = [float(mval_aux[ks].strip(' ')) for ks in range(mvlen)]
    #                                     aux_B['value'] = np.select(conditions, choices, default=mval)
    # =============================================================================
                                        aux_B_i = aux_B_i[aux_i.columns] #adjusting order of columns
                                        aux_i = pd.concat([aux_i, aux_B_i], ignore_index=True)     #concatenating into the first auxiliary df                     

# =============================================================================
#                                 aux_1 = aux_1[aux_m.columns] #adjusting order of columns
#                                 aux_m = pd.concat([aux_m, aux_1], ignore_index=True) #concatenating into the second auxiliary df 
#                         aux_m = aux_m[aux.columns]
#                         aux = pd.concat([aux, aux_m], ignore_index=True)
# =============================================================================
                                        
                    if list(itertools.chain(*BLUES_tec_df.loc[idx[tn, [nod], :], 'node_in'].values.tolist())) == [nod]: #checking if node out is equal to the current node
                        aux_i['node_in'] = aux_i.node #if yes they are equal
                    else:   #in case they are different
                        conditions = [np.logical_and(aux_i['node'] == nod, aux_i['mode']== m1 ) for m1 in mod] #condition to match current node and mode
                        choices = list(itertools.chain(*BLUES_tec_df.loc[idx[tn, [nod], :], 'node_in'].values.tolist())) #choices that get the node out from database dataframe
                        aux_i['node_in'] = np.select(conditions, choices, default='ERRO') #replacing according to conditions and choices
                    if tn not in ldr_tecsv: # it is not a ldr tec
                        aux_i['time_in'] = aux_i.time #time in is equal time which is 'year' in this case
                    else: #it is a ldr tec
                        if tn in first_ldr: #if technology is not amongst the ldr tecs or is IN first first ldr 
                            aux_i['time_in'] = 'year' #time in is year
                        else: #else time out is time for ldr commodity-level with ldr
                            conditions = [((aux_i['level'] == lc[0]) & (aux_i['commodity'] == lc[1]))  for lc in ldr_list] #get all ldr commodity level
                            choices = [aux_i['time'].values.tolist() for i in range(len(conditions))] #choice according to computed value from ldr values adjusted by cf/optm
                            aux_i['time_in'] = np.select(conditions, choices, default='year') #replacing value according to conditions/choices                    

                    aux_i = aux_i[inpt_df.columns] #adjusting order of columns
                    inpt_df = pd.concat([inpt_df, aux_i], ignore_index=True)  #concatenating into the final df 
    ##########################################################################################################################      
    # =============================================================================
    #     #MAIN OUTPUT VALUES
    # =============================================================================
            main_val_aux = BLUES_tec_df.loc[idx[tn, [nod], :], idx['output_values', '1']].values.tolist() #getting main output value
            main_val1 = [[vv.strip('[]')] for vv in main_val_aux if vv.strip('[]') != 'nan' ] #adjust it
            if 'nan' not in main_val1:
                main_val2 = [[float(vv.strip("'")) for vv in v[0].split(', ')] for v in main_val1] #converting to value
                main_val = [] #IT IS USED TO ADJUST OTHER VALUES
                mmo_mode = 0 
                for main_val_adj in main_val2:
                    mmo_mode+=1
                    #mmo_mode = main_val2.index(main_val_adj) + 1 #getting mode
                    if len(main_val_adj) > count_y: #if there is more than one value
                        main_val_adj = main_val_adj[count_y:] #adjust based on count y
                    while len(main_val_adj) < len(year_act): #make length of years and values match
                        main_val_adj = main_val_adj + [main_val_adj[-1]] #repeating the last one
                        
                    #main_val_adj = main_val_adj[:len(year_act)]
                    main_val.append(main_val_adj[:len(year_act)] ) #append into main_val
                    
                    for mvalmo, y in zip(main_val_adj, year_act): #generating main val df
                        main_val_df = expand_grid_name(['node', 'year_act','mode', 'value'], #dataframe columns
                                                                              [nod], [y], [mmo_mode], [mvalmo]).merge(      #dataframe values
                                                                                                              vtg_year_time, how='left', on='year_act')   #merging with time auxiliary df 
                        main_val_df = main_val_df[main_output_df.columns]
                        main_output_df = pd.concat([main_output_df, main_val_df], ignore_index=True)
    # =============================================================================
    #     #OUTPUT
    # =============================================================================
            com_aux = BLUES_tec_df.loc[idx[tn, [nod], :], 'output_commodities'].values.tolist() #getting all commoditity values
            lvl_aux = BLUES_tec_df.loc[idx[tn, [nod], :], 'output_levels'].values.tolist() #getting all level values
            com = [[c for c in cc if c!= 'nan'] for cc in com_aux] #removing nan values and keeping each mode in a sublist
            lvl = [[l for l in ll if l  != 'nan'] for ll in lvl_aux] #removing nan values and keeping each mode in a sublist
    
            if len(com) != len(lvl): #checking if number of commoditties and levels match in the mode aspect
                warning_list.append('OUT_check why len com is different from len lvl:  '+str(tn)+' - '+str(nod)) # if not adding warn into warning list
                print('OUT_check why len com is different from len lvl:\t'+tn)
                #break
            else:  #go ahead
                
                val_aux = BLUES_tec_df.loc[idx[tn, [nod], :], 'output_values'].values.tolist() #geting output values
                val = [[v.strip('[]') for v in vv if v != 'nan'] for vv in val_aux] #adjusting them

                if len(val) != len(lvl): #checking if number of commoditties and levels match in the mode aspect
                    warning_list.append('OUT_check why len val is different from len lvl/com:  '+str(tn)+' - '+str(nod)) # if not adding warn into warning list
                    #print('OUT_check why len val is different from len lvl/com:\t'+tn)
                else: #go ahead
        # =================================================================================================
                    aux_o = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode',  'commodity', 'level', 'time', 'value']) #auxiliary empty dataframe
                    for m,c1,l1,mv1 in zip(mod, com, lvl, val): # looping through mode, commodity, level and values
                        #aux_m = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode',  'commodity', 'level', 'time', 'value']) #auxiliary empty dataframe
                        if len(mv1) != len(l1) != len(c1): # check if the number of values, levels and commoditties in each list matches  in the mode level
                            warning_list.append('OUT_check len of each list from val, com and lvl:  '+str(tn)+' - '+str(nod)) #if it is not -- adding info into warning list
                            #print('OUT_check len of each list from val, com and lvl:\t'+tn)
                        else: #if all values length match
                            for c, l, mv in zip(c1,l1,mv1): #looping through each values of commodity, levl and value i the mode subset
                                #aux_1 = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode',  'commodity', 'level', 'time', 'value']) #auxiliary empty dataframe
                                mval1 = [float(v.strip().strip("'")) for v in mv.strip('[]').split(',')] #converting value or list of vlues to float or list of float #strip remove first and last character of string
                                if len(mval1) == 1: #it means the values does not vary over years
                                    aux_B_o = expand_grid_name(['node', 'vintage','mode', 'commodity', 'level', 'value'], #dataframe columns
                                                                 [nod], vtgs, [m], [c], [l], mval1).merge(              #dataframe values
                                                                                                 vtg_year_time, how='left', on='vintage') #merging with time auxiliary df
                                    aux_B_o = aux_B_o[aux_o.columns] #adjusting order of columns
                                    aux_o = pd.concat([aux_o, aux_B_o], ignore_index=True) #concatenating into the first auxiliary df
                                else: # in this case the values vary according to the years
                                    if len(mval1) > count_y:
                                        mval1 = mval1[count_y:] #removing value according to the first year_act
                                    while len(mval1) < len(year_act): # adjustin the length of the list values
                                        mval1 = mval1 + [mval1[-1]] # making it matches the length of modelled years
                                    for mval, y in zip(mval1, year_act): #looping thorugh years and values
                                        aux_B_o = expand_grid_name(['node', 'year_act','mode', 'commodity', 'level', 'value'], #dataframe columns
                                                                              [nod], [y], [m], [c], [l], [mval]).merge(      #dataframe values
                                                                                                              vtg_year_time, how='left', on='year_act')   #merging with time auxiliary df                               
    
    # =============================================================================
    #                                     conditions =  [ (aux_B['vintage'] == vtgs[ky]) for ky in list(range(mvlen))]
    #                                     choices = [float(mval_aux[ks].strip(' ')) for ks in range(mvlen)]
    #                                     aux_B['value'] = np.select(conditions, choices, default=mval)
    # =============================================================================
                                        aux_B_o = aux_B_o[aux_o.columns] #adjusting order of columns
                                        aux_o = pd.concat([aux_o, aux_B_o], ignore_index=True)     #concatenating into the first auxiliary df                     
                                
# =============================================================================
#                                 aux_1 = aux_1[aux_m.columns] #adjusting order of columns
#                                 aux_m = pd.concat([aux_m, aux_1], ignore_index=True) #concatenating into the second auxiliary df 
#                             aux_m = aux_m[aux.columns]
#                             aux = pd.concat([aux, aux_m], ignore_index=True)
# =============================================================================
                            
                    if list(itertools.chain(*BLUES_tec_df.loc[idx[tn, [nod], :], 'node_out'].values.tolist())) == [nod]: #checking if node out is equal to the current node
                        aux_o['node_out'] = aux_o.node #if yes they are equal
                    else:  #in case they are different
                        conditions = [np.logical_and(aux_o['node'] == nod, aux_o['mode'] == m1 ) for m1 in mod] #condition to match current node and mode
                        choices = list(itertools.chain(*BLUES_tec_df.loc[idx[tn, [nod], :], 'node_out'].values.tolist())) #choices that get the node out from database dataframe
                        aux_o['node_out'] = np.select(conditions, choices, default='ERRO') #replacing according to conditions and choices
                    if tn not in ldr_tecsv: # it is not a ldr tec
                        aux_o['time_out'] = aux_o.time #time in is equal time which is 'year' in this case
                    else: #it is a ldr tec
                        if tn in final_ldr: #if technology is a final ldr the time out is 'year'
                            aux_o['time_out'] = 'year' #time in is year
                        else: #else time out is time for ldr commodity-level with ldr
                           
                            conditions = [((aux_o['level'] == lc[0]) & (aux_o['commodity'] == lc[1]))  for lc in ldr_list] #get all ldr commodity level
                            choices = [aux_o['time'].values.tolist() for i in range(len(conditions))] #choice according to computed value from ldr values adjusted by cf/optm
                            aux_o['time_out'] = np.select(conditions, choices, default='year') #replacing value according to conditions/choices                    
                            
                    aux_o = aux_o[output_df.columns] #adjusting order of columns
                    output_df = pd.concat([output_df, aux_o], ignore_index=True)   #concatenating into the final df   

                    #duplicateRowsDF = output_df[output_df.duplicated()]
                    #duplicateRowsDF = inpt_df[inpt_df.duplicated()]            
    ########################################################################################################################## 
    # =============================================================================
    # #CAPACITY FACTOR
    # =============================================================================
            #aux = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'time', 'value']) #creatigng auxiliary df
            aux_cf_B = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'time', 'value']) #auxiliary df
            #check if any value related to capacity factor is nan 
            if str(BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0]) == 'nan' and str(BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0]) == 'nan' and tn not in list(ldr_info[nod+'_cf_norm'].keys()):
                pass # both values are nan --> pass
# =============================================================================
#             
#                 CF_v = 1 #if both are nan setting cpacity factor as 1
#                 aux_CF = expand_grid_name(['node', 'vintage', 'value'], #dataframe columns
#                                              [nod], vtgs, [CF_v]).merge( #df values
#                                                      vtg_year_time, how='left', on='vintage') #merging into aux time df
#                 aux_CF = aux_CF[aux_B.columns] #oredering columns
#                 aux_B = pd.concat([aux_B, aux_CF], ignore_index=True) #concatenating dfs
            else: # at least one of the parameters are not nan
# =============================================================================
                cf_dic = {} #auxiliary df dictionary
                if tn in list(ldr_info[nod+'_cf_norm'].keys()): #check if tn has an explicit ldr representation from ldr file
                    cf_t = ldr_info[nod+'_cf_norm'][tn] #if yes getting the normalized value
                    if str(BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0]) == 'nan' and str(BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0]) == 'nan':
                        cf_cf = [cf_t.mean()] #if there is no optm or plf set the average value of ldr info as de ldv average cf value
                        #check if capacity facotr parameter is nan and optm is not nan
                    elif str(BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0]) == 'nan' and str(BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0]) != 'nan':
                        #in this case cf will be represented by optm
                        cf_cf = [ float(cfv.strip().strip("'")) for cfv in  BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0].strip('[]').split(',')]
                    #otherwise
                    elif str(BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0]) == 'nan' and str(BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0]) != 'nan':
                        #cf will be reprsented by capcaity factor parameter
                        cf_cf = [ float(cfv.strip().strip("'")) for cfv in  BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0].strip('[]').split(',')]
                    
                    else:#in case both them are not nan we are using the highest one
                        optm_v = [ float(cfv.strip().strip("'")) for cfv in  BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0].strip('[]').split(',')] #getting otpm value
                        cf_v = [ float(cfv.strip().strip("'")) for cfv in  BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0].strip('[]').split(',')] #getting capacity factor value
                        #cf_cf =  [(g + h) / 2 for g, h in zip(cf_v, optm_v)]
                        if stat.mean(cf_v) > stat.mean(optm_v): #check if mean of cf i gerater than optm average
                            cf_cf = cf_v #in this case cf is the chosen one
                        else: #otherwise
                            cf_cf = optm_v #optm is the one that ll be used
                    if len(cf_cf) == 1:  #if there is only one cf value it means it does not vary over time
                        cf_cf_v_aux = float(cf_cf[0]) #getting the value as float
                        del cf_cf #deleting cf cf variable
                        cf_t_v = cf_t * cf_cf_v_aux / cf_t.mean() #adjusting the averaege of ldr capacity factor according to the cf/optm value
                        
                        cf_t_v = list(cf_t_v)
                        cf_t_v1 = [1 if cc>1 else float_round(cc, 2, ceil) for cc in cf_t_v]
                        aux_CF = expand_grid_name(['node', 'vintage', 'value'], #df columns
                                                 [nod], vtgs, [9999999]).merge( #df values
                                                         vtg_year_time, how='left', on='vintage')  #merging into time aux df
                        
                        #conditions for each time (monthly in this case)
    # =============================================================================
                         #adjust it for different time options #check if ldr ingo has to be adjusted as well
    # =============================================================================
                        conditions = [(aux_CF['time']==n) for n in time] #defining conditions as the time
                        choices = cf_t_v1 #choice according to computed value from ldr values adjusted by cf/optm
                        aux_CF['value'] = np.select(conditions, choices, default=9999) #replacing value according to conditions/choices
                        aux_CF = aux_CF[aux_cf_B.columns] #oredering columns
                        aux_cf_B = pd.concat([aux_cf_B, aux_CF], ignore_index=True) #concatenating dfs
    
                    else: #values vary over time
                        if len(cf_cf) > count_y:
                            cf_cf = cf_cf[count_y:] #removing value according to the first year_act
                        while len(cf_cf) < len(year_act): #check if cf values length and year length match
                            cf_cf = cf_cf + [cf_cf[-1]] #making them match by repeating the last values 
                        for y in year_act: # loopin through years
                            cf_auxn = cf_t * cf_cf[year_act.index(y)] / cf_t.mean() #auxiliary cf dictionary for year y is the ldr value pattern adjusted according to the cf/optm value for this year
                            cf_auxn = list(cf_auxn)
                            cf_dic[y] = [1 if cc>1 else float_round(cc, 2, ceil) for cc in cf_auxn]
                            aux_CF = expand_grid_name(['node', 'year_act', 'value'], #df columns
                                                 [nod], [y], [9999999]).merge( #df values
                                                         vtg_year_time, how='left', on='year_act') #merging with time aux df
    
                        #conditions for each time (monthly in this case)
    # =============================================================================
                         #adjust it for different time options #check if ldr ingo has to be adjusted as well
    # =============================================================================
                            conditions = [(aux_CF['time']==n) for n in time] #defining conditions as the time
                            choices = cf_dic[y] #choice according to computed value from ldr values adjusted by cf/optm
                            aux_CF['value'] = np.select(conditions, choices, default=9999) #replacing value according to conditions/choices
                            aux_CF = aux_CF[aux_cf_B.columns] #oredering columns
                            aux_cf_B = pd.concat([aux_cf_B, aux_CF], ignore_index=True) #concatenating dfs
    
                else: #in case it has no explicit ldr representation
                    #check if optm AND capacity factor are not nan
                    if str(BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0]) != 'nan' and str(BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0])!= 'nan': 
                        optm_v = [ float(cfv.strip().strip("'")) for cfv in  BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0].strip('[]').split(',')] #getting otpm values as list of float values
                        cf_v = [ float(cfv.strip().strip("'")) for cfv in  BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0].strip('[]').split(',')] #getting capacity factor value as list of flaot values
                        if len(cf_v) > 1 and len(optm_v) > 1: #if both length are greater than one
                            if len(cf_v) == len(optm_v): #if both have same length
                                #CF =  [(g + h) / 2 for g, h in zip(cf_v, optm_v)] #When the tec has both optm and pll (cf ) I am applying the mean between them. CHECK IF TIS OK
                                if stat.mean(cf_v) > stat.mean(optm_v): # in this case check if capacity factor is higher than optm
                                    CF = cf_v #use capacity factor
                                else: #otherwise
                                    CF = optm_v #use optm
                            elif len(cf_v) > len(optm_v): #if length of capacity factor is greater than optm
                                CF = cf_v #use capacity factor
                            else: #otherwise
                                CF = optm_v #use optm
                            del cf_v, optm_v #del both
                        elif len(cf_v) > 1 or len(optm_v) > 1: # if length of one of them is gerater than one
                            if stat.mean(cf_v) > stat.mean(optm_v): # in this case check if capacity factor is higher than optm
                                CF = cf_v #use cf
                            else: #otherwise
                                CF = optm_v #use otpm
                        else: #in this case both have length one
                            optm_v = float(optm_v[0]) #get value as float
                            cf_v = float(cf_v[0]) #get value as float
                            CF = [max(cf_v, optm_v)] #get the biggest one
                            del cf_v, optm_v # del variables
                    
                    elif str(BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0]) != 'nan': #check if optm is not nan
                        optm_v = [ float(cfv.strip().strip("'")) for cfv in  BLUES_tec_df.loc[idx[tn, [nod], :], 'optm_t'].values[0][0].strip('[]').split(',')] #get optm values
                        CF = optm_v #use optm
                        del optm_v #del variable
                    elif str(BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0])!= 'nan': #check if capacity factor is not nan
                        cf_v = [ float(cfv.strip().strip("'")) for cfv in  BLUES_tec_df.loc[idx[tn, [nod], :], 'capacity_factor'].values[0][0].strip('[]').split(',')] #get cf values
                        CF = cf_v #use cf
                        del cf_v #del variable
                    if len(CF) >count_y: # len CF IS larger than one remove vintage value
                        CF = CF[count_y:] #removing value according to the first year_act
                        for y in year_act: # looping through vintages
                            while len(CF) < len(year_act): # check if defined CF length matches year all length
                                CF = CF + [CF[-1]] #make them match
                            CF_v = CF[year_act.index(y)] #get the value that correponds to the specific year y
                            aux_CF = expand_grid_name(['node', 'year_act', 'value'],  #df columns
                                                 [nod], [y], [CF_v]).merge( #df values
                                                         vtg_year_time, how='left', on='year_act') #merging into time aux df
                            aux_CF = aux_CF[aux_cf_B.columns] #oredering columns
                            aux_cf_B = pd.concat([aux_cf_B, aux_CF], ignore_index=True) #concatenating into aux df
                    else: #if it  CF length is equal 1 get everything
                        aux_cf_B = expand_grid_name(['node', 'year_act', 'value'],  #df columns
                                             [nod], year_act, CF).merge( #df values
                                                     vtg_year_time, how='left', on='year_act') #merging into time aux df

                        
            aux_cf_B = aux_cf_B[capacity_factor_df.columns] #ordering columns 
            capacity_factor_df = pd.concat([capacity_factor_df, aux_cf_B], ignore_index=True) #concatenating into final df

    ##########################################################################################################################         
    # =============================================================================
    #         #MIN UTILIZATION FACTOR
    # =============================================================================
            
            muf_aux = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'value']) #df aux
            if tn == 'Elect_Final': #if technology is among final ldr tecs
                
                muf = aux_cf_B.copy() #capacity_factor_df.copy() #minimum utilization factor must be euqal capacity factor average over time
                muf['value'] = [float(t) for t in muf.value.tolist()] #to guarantee that all values are as float
                muf_aux = muf.groupby(muf_df.columns.tolist()[:-1], as_index=False)['value'].mean()       #computing mean and excluding time column
                #muf_aux['value'] = muf_aux['value'].values -0.02 #relaxing MUF of ldr technologies 
                muf_aux['value'] = [ float_round(vv, 2, floor) for vv in muf_aux['value']] 
                
# =============================================================================
#             elif tn in int_ldr: #if tn technology is an intermediary ldr tec
#                 muf = capacity_factor_df[capacity_factor_df['node'] == nod].copy() #copying capacity factor
#                 muf_y = muf[muf.time=='year']  #subsetting df for time value equal year
#                 muf_t = muf[muf.time!='year'] #subsetting df for time values different of year
#     
#                 muf_y['value'] = [float(t) for t in muf_y.value.tolist()] #make sure all values are float
#                 muf_t['value'] = [float(t) for t in muf_t.value.tolist()]         #make sure all values are float
#                 muf_aux_y = muf_t.groupby(muf_df.columns.tolist()[:-1], as_index=False)['value'].mean()  #get mean and remove time column
#                 muf_aux_t = muf_t.groupby(muf_df.columns.tolist()[:-1], as_index=False)['value'].mean()  #get mean and remove time column
#                 if muf_aux_y.equals(muf_aux_t): #check if both datafarames , for year time and higher resolution time are equal
#                     muf_aux = muf_aux_y.copy() #if they are equal get one of them
#                 else: #otherwise
#                     warning_list.append('check capacity factor - int_ldr_tec '+str(tn)) #append message into warning list
# =============================================================================
                
            else: #if tec is not among the previous case
                muf_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('minimum_utilization_factor')].values.tolist())] #get values for each mode
                if len(muf_b) == 1: # if it is TRUE it means that the values are the same for all modes
                    if str(muf_b[0][0]) != 'nan': #if it is not nan
                        muf_a = [float(v.strip().strip("'")) for v in muf_b[0][0].strip('[]').split(',')] #get values list as float
                        if len(muf_a) == 1: #it means the values are the same for the entire period
                            aux = expand_grid_name(['node', 'vintage',  'value'], #df columns
                                                               [nod], vtgs,  muf_a ).merge( #df values
                                                                  vtg_year, how='left', on='vintage') #merging into aux time df
                            aux = aux[muf_aux.columns] #oredering columns
                            muf_aux = pd.concat([muf_aux, aux], ignore_index=True) #concatenating into aux df
                        else:           # it means the values vary according to the years 
                            muf_a = [float(sv) for sv in muf_a if str(sv) !='nan'] #removing nan values
                            # check if I need to repeat the values for the 2 last years or if it is not considering the 2 first years
                            if len(muf_a) > count_y:
                                muf_a = muf_a[count_y:] #removing values according to the first year all value
                            while len(muf_a) < len(year_act): #check if length of muf and year all are matching
                                muf_a = muf_a + [muf_a[-1]] #make them match
                            
                            for y,v in zip(year_act, muf_a): #through year and value
                                aux = expand_grid_name(['node', 'year_act',  'value'], #df columns
                                                               [nod], [y], [v] ).merge( #df values
                                                                  vtg_year, how='left', on='year_act') #merging into aux time df
                                aux = aux[muf_aux.columns] #oredering colum
                                muf_aux = pd.concat([muf_aux, aux], ignore_index=True) #concat into aux df
    
                else: # the values are different for each mode           
                    #print("CHECK IT bc for muf values should be the same for all modes")
                    warning_list.append("CHECK IT bc for muf values should be the same for all modes:  "+str(tn)+' - '+str(nod)) # it should not happen I believe || add warn into warning list
            muf_aux = muf_aux[muf_df.columns]#ordering columns
            muf_df = pd.concat([muf_df, muf_aux], ignore_index = True) #concat into final df
    
    ##########################################################################################################################      
    # =============================================================================
    #         #EMISSION FACTOR:
    # =============================================================================
            aux_B = pd.DataFrame(columns=['node',  'year_act',  'mode',  'emission','value']) #creating aux df for emissions
            emi_aux =  [a.strip('[]').split(',') for a in [aa[0] for aa in BLUES_tec_df.loc[idx[tn, [nod], :], 'emissions'].values.tolist()]]
            emi_aux = [[ef.strip().strip("'") for ef in ef2] for ef2 in emi_aux] #get emission type
            emi_aux = [[ef for ef in ef2 if ef != 'nan'] for ef2 in emi_aux]  #removing nans
            if len(emi_aux) >0:#if there is any emission
                emi_val = [ef[0].strip().strip("'") for ef in BLUES_tec_df.loc[idx[tn, [nod], :], 'emission_values'].values] #get emission value strip remove the first and the last character of a string
                emi_val = [ef.strip('[]').split('], [') for ef in emi_val] #split values according to each mode
                emi_val = [[re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", l2) for l2 in l1] for l1 in emi_val] #get number
                emi_val = [[[float(l3.strip("'")) for l3 in l2] for l2 in l1 ] for l1 in emi_val] #converting values into float according to the correct lists for each emission and mode
                for emi1, ev1, m in zip(emi_aux, emi_val, mod): #looping through emission values and modes
                    #mode_out_val = main_val[m-1]
                    
                    #adjusting value based on the main input
# =============================================================================
#                     if len(main_in_val) > 0:    
#                         #mode_out_val = main_val[m-1]
#                         mode_inp_val = main_in_val[m-1]
#                         #ef_io = list(np.array(mode_out_val)/np.array(mode_inp_val) )
#                         ef_io = list(np.array([1] )/np.array(mode_inp_val) )
#                     else:
#                         #ef_io = main_val[m-1]       
# =============================================================================
                    ef_io = list(np.array([1] * len(year_act) ) )
                        
                    for em, ev in zip(emi1, ev1):# looping thorugh emission type and its emission values
                        em = em.split('-')[0]
                        if len(ev) == 1: # if len ev >1 it menas values vary over time
                            ev = ev * len(year_act)
                            
                        if len(ev) > count_y:
                            ev = ev[count_y:] #removing value according to the first year_act
                        while len(ev) < len(year_act): #check if values length between year_act and emission match 
                            ev =  ev + [ev[-1]] #make them match
                        for y,e, main_o_v in zip(year_act, ev, ef_io): #looping thorugh years and values
                            aux = expand_grid_name(['node', 'year_act', 'mode','emission',  'value'], #df columns 
                                                           [nod], [y], [m], [em], [e*main_o_v] )# df values
                            aux = aux[aux_B.columns]#oredring columns
                            aux_B = pd.concat([aux_B, aux], ignore_index=True) #concatenating into aux df
                
                aux_B = aux_B.merge(vtg_year, how='left', on='year_act')     #merging with aux time df
                aux_B = aux_B[emission_factor_df.columns] #oredering columsn
                emission_factor_df = pd.concat([emission_factor_df, aux_B], ignore_index=True) #concatenating into maind df    
    ##########################################################################################################################  
    # =============================================================================
    #         #INVESTMENT COST
    # =============================================================================      
            inv_cost_aux = pd.DataFrame(columns=['node',  'vintage',  'value']) #aux df

            inv_cost_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('investment_cost')].values.tolist())] #getting inv cost values

            if len(inv_cost_b) == 1: # if it is TRUE it means that the values are the same for all modes
                inv_cost_a = [float(v.strip().strip("'")) for v in inv_cost_b[0][0].strip('[]').split(',')] # get list of value as float
                    
            else: # the values are different for each mode           
                print("CHECK IT bc for inv_cost values should be the same for all modes")
                warning_list.append("CHECK IT bc for inv_cost values should be the same for all modes:  "+str(tn)+' - '+str(nod)) # it should not happen -->report into warning list                    
            
            if tn in ldr_original: #it is here to guarantee that ldr tecs will have capacity parameter associated to them
                if inv_cost_a == [0]:
                    inv_cost_a = [0.1]
                   
            if str(inv_cost_a[0]) != 'nan':
                inv_v = vtgs[:]
                    
                if len(inv_cost_a) == 1: #it means the values are the same for the entire period
                    
                    aux = expand_grid_name(['node', 'vintage',  'value'],  #df columns
                                                       [nod], inv_v,  inv_cost_a ) # df values
                    aux = aux[inv_cost_aux.columns]#ordering columns
                    inv_cost_aux = pd.concat([inv_cost_aux, aux], ignore_index=True) #concatenating into aux df
                    
                else:           # it means the values vary according to the years 
                    inv_cost_a = [float(sv) for sv in inv_cost_a if str(sv) !='nan']#removing nan values
    
                    while len(inv_cost_a) < len(inv_v):#check if length vintages and len values amtch
                        inv_cost_a = inv_cost_a + [inv_cost_a[-1]] #make them match
                    
                    for vin,v in zip(inv_v, inv_cost_a): #looping through vitnages and values
                        aux = expand_grid_name(['node', 'vintage',  'value'],  #df columns
                                                       [nod], [vin], [v] ) #df values
                        aux = aux[inv_cost_aux.columns]#ordering columns
                        inv_cost_aux = pd.concat([inv_cost_aux, aux], ignore_index=True) #concatenting into aux df   


            inv_cost_aux = inv_cost_aux[inv_cost_df.columns]    #ordering columns
            inv_cost_df = pd.concat([inv_cost_df, inv_cost_aux], ignore_index = True) #concatenating into final df
    
    
    ##########################################################################################################################         
    # =============================================================================
    #         #VARIABLE COST
    # =============================================================================      
            var_cost_aux = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'mode', 'time', 'value']) #aux var df
            var_cost_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('variable_cost')].values.tolist())] #getting var values in list for each mode
            if len(var_cost_b) == 1: # if it is TRUE it means that the values are the same for all modes
                if str(var_cost_b[0][0]) == 'nan': # if its not nan
                    #var_cost_a = [0.1]
                    var_cost_a = 'nan'
                else:
                    var_cost_a = [float(v.strip().strip("'")) for v in var_cost_b[0][0].strip('[]').split(',')] #get values in list as float
                if var_cost_a != 'nan' :
                    
                    if len(var_cost_a) == 1: #it means the values are the same for the entire period
                        var_cost_a = var_cost_a * len(year_act)
    
                    else:           # it means the values vary according to the years 
                        var_cost_a = [float(sv) for sv in var_cost_a if str(sv) !='nan'] #values in a list as float
    
                    if len(var_cost_a) > count_y:
                        var_cost_a = var_cost_a[count_y:] #removing value according to the first year_act
                    while len(var_cost_a) < len(year_act): #check if length of var values and year_act match
                        var_cost_a = var_cost_a + [var_cost_a[-1]] #make them match
                    
# =============================================================================
                    #INCLUDING EFFICIENCY TO COMPUTE VAR COST  
#                    if len(main_in_val) > 0:    
#                        mode_out_val = main_val[0]
#                        mode_inp_val = main_in_val[0]
#                        vef_io = list(np.array(mode_out_val)/np.array(mode_inp_val) )
#                    else:
#                        vef_io = main_val[m-1]                        
# =============================================================================
                    
                    #for y,var,m_out_val in zip(year_act, var_cost_a, vef_io): #looping through year a var values
                    for y,var in zip(year_act, var_cost_a): #looping through year a var values
                        aux = expand_grid_name(['node', 'year_act', 'mode', 'value'],  #df columns
                                                       [nod], [y], mod, 
                                                       #[var*m_out_val] ).merge( #df values
                                                       [var] ).merge( #df values
                                                      vtg_year_time, how='left', on='year_act') #merging into aux df time
                        aux = aux[var_cost_aux.columns] #oredering columns
                        var_cost_aux = pd.concat([var_cost_aux, aux], ignore_index=True) #concatenating into var aux
            
            else: # the values are different for each mode
                var_cost_a = [0] #maybe add it
                for va,m in zip(var_cost_b, mod): # looping thorugh values and node
                    var_val = va[0].strip('[]').split(',') #removing brackets and spliting values
                    var_val = [float(vv.strip("'")) for vv in var_val] #get values as float
                    
                    if len(var_val) == 1: # it means values vary over time
                        var_val1 = var_val * len(year_act)
                    else:
                        var_val1 = var_val[:]
                        
                    if len(var_val1) > count_y:
                        var_val1 =  var_val1[count_y:]  #removing value according to the first year_act
                    while len(var_val1) < len(year_act): #check if length of var values and year_act match
                        var_val1 =  var_val1 + [var_val1[-1]]  #make them match
                    
# =============================================================================
                    #INCLUDING EFFICIENCY TO COMPUTE VAR COST  
#                    if len(main_in_val) > 0:                    
#                        mode_out_val = main_val[m-1]
#                        mode_inp_val = main_in_val[m-1]
#                        vef_io = list(np.array(mode_out_val)/np.array(mode_inp_val) )
#                    else:
#                        vef_io = main_val[m-1]
# =============================================================================
                    
                    #for y,var,m_out_val in zip(year_act, var_val1, vef_io): #looping through year a var values
                    for y,var in zip(year_act, var_val1): #looping through year a var values
                        aux = expand_grid_name(['node', 'year_act', 'mode', 'value'], #df columns
                                                       [nod], [y], [m], 
                                                       #[var*m_out_val] ).merge( #df values
                                                       [var] ).merge( #df values
                                                      vtg_year_time, how='left', on='year_act')  #merging into aux df time
                        
                        aux = aux[var_cost_aux.columns] #oredering columns
                        var_cost_aux = pd.concat([var_cost_aux, aux], ignore_index=True) #concatenating into var aux

            var_cost_aux = var_cost_aux[var_cost_df.columns]    #oredering columns
            var_cost_df = pd.concat([var_cost_df, var_cost_aux], ignore_index = True) #concatenating into main df
    
    ##########################################################################################################################         
    # =============================================================================
    #         #FIXED COST
    # =============================================================================      
            fix_cost_aux = pd.DataFrame(columns=['node',  'vintage', 'year_act', 'value']) #aux df
            fix_cost_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('fixed_cost')].values.tolist())] #aux df
            if len(fix_cost_b) == 1: # if it is TRUE it means that the values are the same for all modes
                if str(fix_cost_b[0][0]) != 'nan': # if it is not nan
                    fix_cost_a = [float(v.strip().strip("'")) for v in fix_cost_b[0][0].strip('[]').split(',')] #get value as list of float
                    if len(fix_cost_a) == 1: #it means the values are the same for the entire period
                        aux = expand_grid_name(['node', 'vintage',  'value'], #df columns
                                                           [nod], vtgs,  fix_cost_a ).merge( #df values
                                                              vtg_year, how='left', on='vintage') #merging into aus time df
                        aux = aux[fix_cost_aux.columns] #ordering columns
                        fix_cost_aux = pd.concat([fix_cost_aux, aux], ignore_index=True) # concatenating into fixed cost aux
                    else:           # it means the values vary according to the years 
                        fix_cost_a = [float(sv) for sv in fix_cost_a if str(sv) !='nan'] # removing nan values
    
                        if len(fix_cost_a) > count_y:
                            fix_cost_a = fix_cost_a[count_y:] #removing value according to the first year_act
                        while len(fix_cost_a) < len(year_act): #check if length of var values and year_act match
                            fix_cost_a = fix_cost_a + [fix_cost_a[-1]] #make them match
                        
                        for y,v in zip(year_act, fix_cost_a):#looping through year all and values
                            aux = expand_grid_name(['node', 'year_act',  'value'],  #df columns
                                                           [nod], [y], [v] ).merge( #df values
                                                              vtg_year, how='left', on='year_act') #merging into aux time df
                            aux = aux[fix_cost_aux.columns] #oredering columns
                            fix_cost_aux = pd.concat([fix_cost_aux, aux], ignore_index=True) #concatenating into aus df
            
            else: # the values are different for each mode           
                #print("CHECK IT bc for fix_cost values should be the same for all modes")
                warning_list.append("CHECK IT bc for fix_cost values should be the same for all modes:  "+str(tn)+' - '+str(nod)) #it shoudl not happen add info into warning list
            fix_cost_aux = fix_cost_aux[fix_cost_df.columns]    #ordering columns
            fix_cost_df = pd.concat([fix_cost_df, fix_cost_aux], ignore_index = True) #concatenating into maind final df
    ##########################################################################################################################   
    # =============================================================================
    #         #HISTORICAL NEW CAPACITY HISC 
    # =============================================================================
            aux = pd.DataFrame(columns=['node',  'vintage',  'value']) #aux df
            #if tn1 not in tn_old_out: #because if tn1 equals tn_old_out the historical new capacity has been considered in the old tec
            #if tn1 == tn:

            hisc_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_new_capacity')].values.tolist())] #getting historical values
            hisc_y = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_new_capacity_years')].values.tolist())] #getting historical years
                
            if len(hisc_b)  == 1:
                if str(hisc_b[0][0]) != 'nan':
                    if vtgs[0] == years[0]: # it means the tec can be built in the historical period


                        hisc_y = [int(ttv) for ttv in hisc_y[0][0].strip('[]').split(',') if ttv != 'nan'] #removing nan and get values as int
                        hisc_b = [float(ttv.strip().strip("'")) for ttv in hisc_b[0][0].strip('[]').split(',') if ttv != 'nan'] #removing nan and get values as float

                        hisc_b2 = [] #aux list for values
                        #hisc_y2 = [yy for yy in hisc_y if yy>first_hist_new_cap] #values later greater than 2005 # not used without old tecs
                        for hh in hisc_y:   #looping through values defined                 
                            hisc_b2.append(hisc_b[hisc_y.index(hh)]) ##aggreagating all values related to historical years

                        hisc_v = float_round((sum(hisc_b2)/first_d_period), 0, ceil) #get the sum and divide by period in case lft is greater than last d period

                        aux = expand_grid_name(['node', 'vintage',  'value'],  #df columns
                                                               [nod], [vtgs[0]],  [hisc_v] ) #df values                

                    else:
                        warning_list.append("CHECK IT bc tec vtgs[0] is different from years[0]:  "+str(tn1)+' - '+str(nod)) #add into warning list
                    

            else:
                warning_list.append("CHECK IT bc for hisc values should be the same for all modes:  "+str(tn1)+' - '+str(nod)) #add into warning list
        

            aux = aux[historical_new_capacity_df.columns]    #ordering columns
            historical_new_capacity_df = pd.concat([historical_new_capacity_df, aux], ignore_index = True) #concat into final df    
            
    ##########################################################################################################################         
    # =============================================================================
    #         #BOUND INSTALLED BDI LO bound on total installed capacity LO
    # =============================================================================    
            bdi_lo_aux = pd.DataFrame(columns=['node',  'year_act', 'value']) #aux df    
            if tn1 == tn: #it means it is not an old tec --> old tecs do not have bdi

                if BLUES_tec_df.loc[idx[tn, [nod], :], ('land_constraint')].values.tolist()[0][0] != 'nan': #if land constraint is not nan
                    bdi_v = BLUES_tec_df.loc[idx[tn, [nod], :], ('land_constraint')].values.tolist()[0][0].strip('[]').strip("'") #get land constraint values
                    if bdi_v in [ll for ll in land_rel_df.relation.tolist() if ll.endswith('T')]: # land constraint value is in land relation df list of values that endswith T (it means that this cons is for the total balance of the related tpe of land)
                        bdi_lo_b = [[str([a.strip("'").strip() for a in land_rel_df[(land_rel_df['node']==nod) & (land_rel_df['relation'] == bdi_v)]['lower'].values[0].strip("[]").split(',')])]] #get the lower value from land relation df for this constraint


                ## ?????? IT IS NOT READY IF WE SET 2010 TO BE HISTORICAL add historical related with 2010 for the case we consider this year as historical
        # =============================================================================
        #               DEFINING HISTORICAL LAND AS BDI LO/UP IN THE FYEAR
        # =============================================================================
                        if BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_land_constraint')].values.tolist()[0][0] != 'nan': #if land constraint is not nan
                            hisc_v = BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_land_constraint')].values.tolist()[0][0].strip('[]').strip("'") #get land constraint values
                            if hisc_v in [ll for ll in land_rel_df.relation.tolist() if ll.endswith('a')]:  # land constraint value is in land relation df list of values that endswith T (it means that this cons is for the total balance of the related tpe of land)
                                hisc_bup = [a.strip("'").strip() for a in land_rel_df[(land_rel_df['node']==nod) & (land_rel_df['relation'] == hisc_v)]['upper'].values[0].strip("[]").split(',')] #get the upper value from land relation df for this constraint
                                hisc_blo = [a.strip("'").strip() for a in land_rel_df[(land_rel_df['node']==nod) & (land_rel_df['relation'] == hisc_v)]['lower'].values[0].strip("[]").split(',')] #get the lower value from land relation df for this constraint
                                if hisc_bup == hisc_blo:
                                    hisc_b = [[hisc_blo[0]]]
                                if len(hisc_b) >1: #values should be equal for historical capacity
                                    warning_list.append("CHECK IT bc for hisc values should be the same for all modes:  "+str(tn1)+' - '+str(nod)) #add into warning list
                                else:
                                    hisc_b = [float(ttv.strip().strip("'")) for ttv in hisc_b[0][0].strip('[]').split(',') if ttv != 'nan'] #removing nan and get values as float
                                            
                            else:
                                warning_list.append("CHECK IT bc for historical land_contraint upper and lower values are not equal:  "+str(tn)+' - '+str(nod)) #warn into warning list     
                        else: #the value is not in the land relation list
                            #print("CHECK IT bc for land_contraint values should be in the list")
                            warning_list.append("CHECK IT bc for hisc land_contraint values should be in the list:  "+str(tn)+' - '+str(nod)) #warn into warning list

                    else: #the value is not in the land relation list
                        #print("CHECK IT bc for land_contraint values should be in the list")
                        warning_list.append("CHECK IT bc for land_contraint values should be in the list:  "+str(tn)+' - '+str(nod)) #warn into warning list

                else: #there is no land constraint
                    bdi_lo_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('bound_total_capacity_lo')].values.tolist())] #get bound total cap lo values

                if len(bdi_lo_b) == 1: # if it is TRUE it means that the values are the same for all modes
                    if str(bdi_lo_b[0][0]) != 'nan': #if not nan
                        bdi_lo_b = [re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", bdi_lo_b[0][0])]
                        bdi_lo_a = [float(v.strip("'").strip().strip("'")) for v in bdi_lo_b[0]] # list of values as float
                        
                        if BLUES_tec_df.loc[idx[tn, [nod], :], ('land_constraint')].values.tolist()[0][0] != 'nan': #check if it is land
                            bdi_lo_a = hisc_b + bdi_lo_a #adding historical as bdi lo
                            
                        if len(bdi_lo_a) == 1: #it means the values are the same for the entire period
                            aux = expand_grid_name(['node', 'year_act',  'value'],  #f columns
                                                               [nod], year_act,  bdi_lo_a ) #df values
                            aux = aux[bdi_lo_aux.columns] #ordering columns
                            bdi_lo_aux = pd.concat([bdi_lo_aux, aux], ignore_index=True) #concatenating df
                        else:           # it means the values vary according to the years 
                            bdi_lo_a = [float(sv) for sv in bdi_lo_a if str(sv) !='nan'] #removing nan
                            if len(bdi_lo_a) > count_y:
                                bdi_lo_a = bdi_lo_a[count_y:] ##removing values according to the first year all value
        
                            while len(bdi_lo_a) < len(year_act): #check if length of muf and year all are matching
                                bdi_lo_a = bdi_lo_a + [bdi_lo_a[-1]]# make them amtch
                            
                            for yr,v in zip(year_act, bdi_lo_a): #looping through
                                aux = expand_grid_name(['node', 'year_act',  'value'], #df columns
                                                               [nod], [yr], [v] ) #df values
                                aux = aux[bdi_lo_aux.columns] #ordering columns
                                bdi_lo_aux = pd.concat([bdi_lo_aux, aux], ignore_index=True) # concat into aux df
        
                else: # the values are different for each mode           
                    #print("CHECK IT bc for bdi_lo values should be the same for all modes")
                    warning_list.append("CHECK IT bc for bdi_lo values should be the same for all modes:  "+str(tn)+' - '+str(nod)) #it shoudl not happen add into warning list
# =============================================================================
#             else: # to guarantee that hisc of old tecs will be kept
#                 if "hisc_b_aux" in locals():
#                     if len(hisc_b_aux) > 0: # adding hisc as bid lo for old tecs with las year
#                         if BLUES_tec_df.loc[idx[tn, [nod], :], ('last_year')].values.tolist()[0][0] != 'nan':
#                             aux = expand_grid_name(['node', 'year_act',  'value'],  #f columns
#                                                                [nod], year_act,  hisc_b_aux ) #df values
#                             aux = aux[bdi_lo_aux.columns] #ordering columns
#                             bdi_lo_aux = pd.concat([bdi_lo_aux, aux], ignore_index=True) #concatenating df                
# =============================================================================
                
            bdi_lo_aux = bdi_lo_aux[bdi_lo_df.columns]    #ordering columns
            bdi_lo_df = pd.concat([bdi_lo_df, bdi_lo_aux], ignore_index = True) #concat into final df

    ##########################################################################################################################         
    # =============================================================================
    #         #BOUND INSTALLED BDI UP bound on total installed capacity UP
    # =============================================================================
            bdi_up_aux = pd.DataFrame(columns=['node',  'year_act', 'value']) #aux df
            if tn1 == tn:     #if tn1 is different of tn it means it is an old technology that shoudl not consider variations on total isntalled capacity       
                if BLUES_tec_df.loc[idx[tn, [nod], :], ('land_constraint')].values.tolist()[0][0] != 'nan': #if land constraint is not nan
                    bdi_v = BLUES_tec_df.loc[idx[tn, [nod], :], ('land_constraint')].values.tolist()[0][0].strip('[]').strip("'") #get land constraint values
                    if bdi_v in [ll for ll in land_rel_df.relation.tolist() if ll.endswith('T')]:  # land constraint value is in land relation df list of values that endswith T (it means that this cons is for the total balance of the related tpe of land)
                        bdi_up_b = [[str([a.strip("'").strip() for a in land_rel_df[(land_rel_df['node']==nod) & (land_rel_df['relation'] == bdi_v)]['upper'].values[0].strip("[]").split(',')] )]]#get the upper value from land relation df for this constraint

                ##IT IS NOT READY IF WE SET 2010 TO BE HISTORICAL add historical related with 2010 for the case we consider this year as historical

                        if BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_land_constraint')].values.tolist()[0][0] != 'nan': #if land constraint is not nan
                            hisc_v = BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_land_constraint')].values.tolist()[0][0].strip('[]').strip("'") #get land constraint values
                            if hisc_v in [ll for ll in land_rel_df.relation.tolist() if ll.endswith('a')]:  # land constraint value is in land relation df list of values that endswith T (it means that this cons is for the total balance of the related tpe of land)
                                hisc_bup = [a.strip("'").strip() for a in land_rel_df[(land_rel_df['node']==nod) & (land_rel_df['relation'] == hisc_v)]['upper'].values[0].strip("[]").split(',')] #get the upper value from land relation df for this constraint
                                hisc_blo = [a.strip("'").strip() for a in land_rel_df[(land_rel_df['node']==nod) & (land_rel_df['relation'] == hisc_v)]['lower'].values[0].strip("[]").split(',')] #get the lower value from land relation df for this constraint
                                if hisc_bup == hisc_blo:
                                    hisc_b = [[hisc_blo[0]]]
                                if len(hisc_b) >1: #values should be equal for historical capacity
                                    warning_list.append("CHECK IT bc for hisc values should be the same for all modes:  "+str(tn1)+' - '+str(nod)) #add into warning list
                                else:
                                    hisc_b = [float(ttv.strip().strip("'")) for ttv in hisc_b[0][0].strip('[]').split(',') if ttv != 'nan'] #removing nan and get values as float
                                            
                            else:
                                warning_list.append("CHECK IT bc for historical land_contraint upper and lower values are not equal:  "+str(tn)+' - '+str(nod)) #warn into warning list     
                        else: #the value is not in the land relation list
                            #print("CHECK IT bc for land_contraint values should be in the list")
                            warning_list.append("CHECK IT bc for hisc_land_contraint values should be in the list:  "+str(tn)+' - '+str(nod)) #warn into warning list


                    else: #the value is not in the land relation list
                        #print("CHECK IT bc for land_contraint values should be in the list")
                        warning_list.append("CHECK IT bc for land_contraint values should be in the list:  "+str(tn)+' - '+str(nod)) #warn into warning list
                else: #there is no land constraint
                    bdi_up_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('bound_total_capacity_up')].values.tolist())] #get bound total cap lo values
                if len(bdi_up_b) == 1: # if it is TRUE it means that the values are the same for all modes
                    if str(bdi_up_b[0][0]) != 'nan': #if not nan
                        bdi_up_b = [re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", bdi_up_b[0][0])]
                        bdi_up_a = [float(v.strip().strip("'")) for v in bdi_up_b[0] ] # list of values as float
                        
                        if BLUES_tec_df.loc[idx[tn, [nod], :], ('land_constraint')].values.tolist()[0][0] != 'nan': #check if it is land
                            bdi_up_a = hisc_b + bdi_up_a #adding historical as bdi up                  
                            
                        if len(bdi_up_a) == 1: #it means the values are the same for the entire period
                            aux = expand_grid_name(['node', 'year_act',  'value'],  #f columns
                                                               [nod], year_act,  bdi_up_a ) #df values
                            aux = aux[bdi_up_aux.columns] #order columns
                            bdi_up_aux = pd.concat([bdi_up_aux, aux], ignore_index=True) #concat into aux df
                        else:           # it means the values vary according to the years 
                            bdi_up_a = [float(sv) for sv in bdi_up_a if str(sv) !='nan'] ##removing nan
                            if len(bdi_up_a) > count_y:
                                bdi_up_a = bdi_up_a[count_y:] ##removing values according to the first year all value
        
                            while len(bdi_up_a) < len(year_act): #check if length of muf and year all are matching
                                bdi_up_a = bdi_up_a + [bdi_up_a[-1]] #make them match
                            
                            for yr,v in zip(year_act, bdi_up_a): #looping through
                                aux = expand_grid_name(['node', 'year_act',  'value'],  #df columns
                                                               [nod], [yr], [v] ) #df values
                                aux = aux[bdi_up_aux.columns] #ordering columns
                                bdi_up_aux = pd.concat([bdi_up_aux, aux], ignore_index=True)  #concat into aux df
        
                else: # the values are different for each mode           
                    #print("CHECK IT bc for bdi_up values should be the same for all modes")
                    warning_list.append("CHECK IT bc for bdi_up values should be the same for all modes:  "+str(tn)+' - '+str(nod)) #it should not happen add into warning list
            bdi_up_aux = bdi_up_aux[bdi_up_df.columns]    #order columns
            bdi_up_df = pd.concat([bdi_up_df, bdi_up_aux], ignore_index = True) #concat into main df
    
    
    ##########################################################################################################################        
    # =============================================================================
    #         #BOUND activity BDA LO bound on ACTIVITY LO
    # =============================================================================
            bda_lo_aux = pd.DataFrame(columns=['node',  'year_act', 'mode', 'time', 'value']) #aux df
            if tn1 == tn:            

                bda_lo_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('bound_activity_lo')].values.tolist())] #getting bound act values by mode
                if len(bda_lo_b) == 1: # if it is TRUE it means that the values are the same for all modes
                    if str(bda_lo_b[0][0]) != 'nan': #if its not nan
                        bda_lo_a = [float(v.strip().strip("'")) for v in bda_lo_b[0][0].strip('[]').split(',')] #values in a list
                        if len(bda_lo_a) == 1: #it means the values are the same for the entire period
                            bda_lo_a = bda_lo_a * len(year_act)
    
                        else:           # it means the values vary according to the years 
                            bda_lo_a = [float(sv) for sv in bda_lo_a if str(sv) !='nan'] #removing nan values
                            
                        if len(bda_lo_a) > count_y:
                            bda_lo_a = bda_lo_a[count_y:] #adjusting value according to first year all
    
                        while len(bda_lo_a) < len(year_act): #check if length of year and values match
                            bda_lo_a = bda_lo_a + [bda_lo_a[-1]]#make them match

# =============================================================================
#                         if len(main_in_val) > 0:
#                             mode_out_val = main_val[0]
#                             mode_inp_val = main_in_val[0]
#                             bda_ef_io = list(np.array(mode_out_val)/np.array(mode_inp_val) )
#                         else:
#                             bda_ef_io = main_val[0]
# =============================================================================
                        
                        #for yr,bdv, bd_ef in zip(year_act, bda_lo_a, bda_ef_io): #looping through years and bda values
                        for yr,bdv in zip(year_act, bda_lo_a): #looping through years and bda values
                            aux = expand_grid_name(['node', 'year_act', 'mode', 'value'],  #df columns
                                                           [nod], [yr], mod, 
                                                           #[bdv/bd_ef] ).merge( #df values
                                                           [bdv] ).merge( #df values
                                                          year_time, how='left', on='year_act') #merging into aux time df
                            aux = aux[bda_lo_aux.columns] #adjusting columns
                            bda_lo_aux = pd.concat([bda_lo_aux, aux], ignore_index=True) #concat into aux df

                else: # the values are different for each mode
                    #warning_list.append(["check bda and var for this tec: " + tn1 + " - " + nod])
                    bda_up_b = BLUES_tec_df.loc[idx[tn, [nod], :], ('bound_activity_lo')].values.tolist() #getting values by mode
                    for bdva,m in zip(bda_lo_b, mod): # looping thorugh values and node
                        if str(bdva[0]) != 'nan':
                            bda_lo_val = bdva[0].strip('[]').split(',') #removing brackets and spliting values
                            bda_lo_val = [float(vv.strip().strip("'")) for vv in bda_lo_val] #get values as float
                            
                            if len(bda_lo_val) == 1: # it means values vary over time
                                bda_lo_val1 = bda_lo_val * len(year_act)
                            else:
                                bda_lo_val1 = bda_lo_val[:]
                                
                            if len(bda_lo_val1) > count_y:
                                bda_lo_val1 =  bda_lo_val1[count_y:]  #removing value according to the first year_act
                                
                            while len(bda_lo_val1) < len(year_act): #check if length of var values and year_act match
                                bda_lo_val1 =  bda_lo_val1 + [bda_lo_val1[-1]]  #make them match
    
# =============================================================================
#                             if len(main_in_val) > 0:                        
#                                 mode_out_val = main_val[m-1]
#                                 mode_inp_val = main_in_val[m-1]
#                                 vef_io = list(np.array(mode_out_val)/np.array(mode_inp_val) )
#                             else:
# =============================================================================
                                #vef_io = main_val[m-1]

                            #for y,var,m_out_val in zip(year_act, bda_lo_val1, vef_io): #looping through year a var values                            
                            for y,var in zip(year_act, bda_lo_val1): #looping through year a var values
                                aux = expand_grid_name(['node', 'year_act', 'mode', 'value'], #df columns
                                                               [nod], [y], [m], 
                                                               #[var/m_out_val] ).merge( #df values
                                                               [var] ).merge( #df values                                                                       
                                                              year_time, how='left', on='year_act')  #merging into aux df time
                                
                                aux = aux[bda_lo_aux.columns] #adjusting columns
                                bda_lo_aux = pd.concat([bda_lo_aux, aux], ignore_index=True) #concat into aux df

                            
            bda_lo_aux = bda_lo_aux[bda_lo_df.columns]    #ordering columns
            bda_lo_df = pd.concat([bda_lo_df, bda_lo_aux], ignore_index = True) #concat into final df
    
    ##########################################################################################################################     
    # =============================================================================
    #         #BOUND activity BDA UP bound on ACTIVITY UP
    # =============================================================================
            bda_up_aux = pd.DataFrame(columns=['node',  'year_act', 'mode', 'time', 'value']) #aux df
            if tn1 == tn:            

                bda_up_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('bound_activity_up')].values.tolist())] #getting bound act values by mode
                if len(bda_up_b) == 1: # if it is TRUE it means that the values are the same for all modes
# =============================================================================
                    if str(bda_up_b[0][0]) == 'nan':#if its not nan
                        if str(var_cost_b[0][0]) != 'nan': #ADDING BDA UP FOR TECS WITH NEGATIVE VOM
                            #if 'var_cost_a' in globals():
                            if var_cost_a[0] < 0:
                                bda_up_b  = [['[999999]']]
# =============================================================================
                            
                    if str(bda_up_b[0][0]) != 'nan':
                        bda_up_a = [float(v.strip().strip("'")) for v in bda_up_b[0][0].strip('[]').split(',')] #values in a list
                        if len(bda_up_a) == 1: #it means the values are the same for the entire period
                            bda_up_a = bda_up_a * len(year_act)
    
                        else:           # it means the values vary according to the years 
                            bda_up_a = [float(sv) for sv in bda_up_a if str(sv) !='nan'] #removing nan values
                            
                        if len(bda_up_a) > count_y:
                            bda_up_a = bda_up_a[count_y:] #adjusting value according to first year all
    
                        while len(bda_up_a) < len(year_act): #check if length of year and values match
                            bda_up_a = bda_up_a + [bda_up_a[-1]]#make them match

# =============================================================================
#                         if len(main_in_val) > 0:   
#                             mode_out_val = main_val[0]
#                             mode_inp_val = main_in_val[0]
#                             bda_ef_io = list(np.array(mode_out_val)/np.array(mode_inp_val) )
#                         else:
#                             bda_ef_io = main_val[0]
# =============================================================================
                        
                        for yr, bdv in zip(year_act, bda_up_a): #upoping through years and bda values
                            aux = expand_grid_name(['node', 'year_act', 'mode', 'value'],  #df columns
                                                           [nod], [yr], mod, 
                                                           #[bdv/bd_ef] ).merge( #df values
                                                           [bdv] ).merge( #df values                                                                   
                                                          year_time, how='left', on='year_act') #merging into aux time df
                            aux = aux[bda_up_aux.columns] #adjusting columns
                            bda_up_aux = pd.concat([bda_up_aux, aux], ignore_index=True) #concat into aux df
        
                else: # the values are different for each mode
                    bda_up_b = BLUES_tec_df.loc[idx[tn, [nod], :], ('bound_activity_up')].values.tolist() #getting values by mode

                    for bdva, m in zip(bda_up_b, mod): # upoping thorugh values and node
                        if str(bdva[0]) != 'nan':
                            bda_up_val = bdva[0].strip('[]').split(',') #removing brackets and spliting values
                            bda_up_val = [float(vv.strip().strip("'")) for vv in bda_up_val] #get values as fupat
                            
                            if len(bda_up_val) == 1: # it means values vary over time
                                bda_up_val1 = bda_up_val * len(year_act)
                            else:
                                bda_up_val1 = bda_up_val[:]     
                                
                            if len(bda_up_val1) > count_y:
                                bda_up_val1 =  bda_up_val1[count_y:]  #removing value according to the first year_act
                                
                            while len(bda_up_val1) < len(year_act): #check if length of var values and year_act match
                                bda_up_val1 =  bda_up_val1 + [bda_up_val1[-1]]  #make them match
    
# =============================================================================
#                             if len(main_in_val) > 0:                           
#                                 mode_out_val = main_val[m-1]
#                                 mode_inp_val = main_in_val[m-1]
#                                 vef_io = list(np.array(mode_out_val)/np.array(mode_inp_val) )
#                             else:
# =============================================================================
                                #vef_io = main_val[m-1]
                            
                            #for y,bdvar,ef_out in zip(year_act, bda_up_val1, vef_io): #upoping through year a var values
                            for y,bdvar in zip(year_act, bda_up_val1): #upoping through year a var values                            
                                aux = expand_grid_name(['node', 'year_act', 'mode', 'value'], #df columns
                                                               [nod], [y], [m], 
                                                               #[bdvar/ef_out] ).merge( #df values
                                                               [bdvar] ).merge( #df values
                                                              year_time, how='left', on='year_act')  #merging into aux df time
                                
                                aux = aux[bda_up_aux.columns] #adjusting columns
                                bda_up_aux = pd.concat([bda_up_aux, aux], ignore_index=True) #concat into aux df
                            
            bda_up_aux = bda_up_aux[bda_up_df.columns]     #ordering columns
            bda_up_df = pd.concat([bda_up_df, bda_up_aux], ignore_index = True)#concat into final df
     
    ##########################################################################################################################         
    # =============================================================================
    #         #BOUND NEW CAPACITY BDC LO bound on NEW installed capacity LO
    # =============================================================================
            bdc_lo_aux = pd.DataFrame(columns=['node',  'vintage', 'value']) #creating aux df
            #if tn1 == tn:           #nto considering  bdc for old tecs

            bdc_lo_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('bound_new_capacity_lo')].values.tolist())] #getting bound new cap lo for each mode
            if len(bdc_lo_b) == 1: # if it is TRUE it means that the values are the same for all modes
                if str(bdc_lo_b[0][0]) != 'nan': #if its not nan
                    bdc_lo_a = [float(v.strip().strip("'")) for v in bdc_lo_b[0][0].strip('[]').split(',')] #get values in a list of float
                    bdcv = [vv for vv in vtgs if vv>=year_act[0]]
                    if len(bdc_lo_a) == 1: #it means the values are the same for the entire period
                        aux = expand_grid_name(['node', 'vintage',  'value'], #df columns
                                                           [nod], bdcv,  bdc_lo_a )#df values
                        aux = aux[bdc_lo_aux.columns]#oredering columns
                        bdc_lo_aux = pd.concat([bdc_lo_aux, aux], ignore_index=True) #concat into aux df
                    else:           # it means the values vary according to the years 
                        bdc_lo_a = [float(sv) for sv in bdc_lo_a if str(sv) !='nan'] #removing nan
                        # check if I need to repeat the values for the 2 last years or if it is not considering the 2 first years
                        
                        #bdc_lo_a =  [bdc_lo_a[0]] +  bdc_lo_a #adjusting due to the different approach from old msg to IX -->it being done earlier by setting bdc =0 in 2010
                        
                        if len(bdc_lo_a) > count_y:
                            bdc_lo_a =  bdc_lo_a[count_y:]
                        while len(bdc_lo_a) < len(bdcv): #check if length of values and vintages match
                            bdc_lo_a = bdc_lo_a + [bdc_lo_a[-1]] #make them match
                        
                        for vin,v in zip(bdcv, bdc_lo_a): #looping through 
                            aux = expand_grid_name(['node', 'vintage',  'value'], #df columns
                                                           [nod], [vin], [v] ) #df values
                            aux = aux[bdc_lo_aux.columns]#oredering columns
                            bdc_lo_aux = pd.concat([bdc_lo_aux, aux], ignore_index=True) #concat into aux df
    
            else: # the values are different for each mode           
                #print("CHECK IT bc for bdc_lo values should be the same for all modes")
                warning_list.append("CHECK IT bc for bdc_lo values should be the same for all modes:  "+str(tn)+' - '+str(nod)) #it should not happen add into warning list
        
            bdc_lo_aux = bdc_lo_aux[bdc_lo_df.columns]     #adjust columns
            bdc_lo_df = pd.concat([bdc_lo_df, bdc_lo_aux], ignore_index = True) #concat into final df
    
    ##########################################################################################################################             
    # =============================================================================
    #         #BOUND NEW CAPACITY BDC UP bound on NEW installed capacity UP
    # =============================================================================
            if BLUES_tec_df.loc[idx[tn, [nod], :], ('add_bdc')].values.tolist()[0][0] == True:
                bdc_up_aux = pd.DataFrame(np.array([ [nod, 2010, 0 ]]), columns=['node',  'vintage', 'value'] ) #creating aux df
            else:
                bdc_up_aux = pd.DataFrame( columns=['node',  'vintage', 'value'] ) #creating aux df
            #if tn1 == tn:            
            bdc_up_b = [list(item) for item in unique(tuple(row) for row in BLUES_tec_df.loc[idx[tn, [nod], :], ('bound_new_capacity_up')].values.tolist())] #getting bound new cap lo for each mode
            if len(bdc_up_b) == 1: # if it is TRUE it means that the values are the same for all modes
                if str(bdc_up_b[0][0]) != 'nan': #if its not nan
                    bdc_up_a = [float(v.strip().strip("'")) for v in bdc_up_b[0][0].strip('[]').split(',')] #get values in a list of floa
                    bdcv = [vv for vv in vtgs if vv>=year_act[0]] #greater than year ac bc first year act cannot build
                    if len(bdc_up_a) == 1: #it means the values are the same for the entire period
                        aux = expand_grid_name(['node', 'vintage',  'value'],  #df columns
                                                           [nod], bdcv,  bdc_up_a ) #df values
                        aux = aux[bdc_up_aux.columns] #oredering columns
                        bdc_up_aux = pd.concat([bdc_up_aux, aux], ignore_index=True) #concat into aux df
                    else:           # it means the values vary according to the years 
                        bdc_up_a = [float(sv) for sv in bdc_up_a if str(sv) !='nan'] #removing nan
                        # check if I need to repeat the values for the 2 last years or if it is not considering the 2 first years
                        
                        #bdc_up_a =  [bdc_up_a[0]] +  bdc_up_a #adjusting due to the different approach from old msg to IX  -->it is being done earlier by setting bdc =0 in 2010
                        
                        if len(bdc_up_a) > count_y:
                            bdc_up_a =  bdc_up_a[count_y:]
                        while len(bdc_up_a) < len(bdcv): #check if length of values and vintages match
                            bdc_up_a = bdc_up_a + [bdc_up_a[-1]] #make them match
                        
                        for vin,v in zip(bdcv, bdc_up_a): #looping through 
                            aux = expand_grid_name(['node', 'vintage',  'value'],  #df columns
                                                           [nod], [vin], [v] )#df values
                            aux = aux[bdc_up_aux.columns] #oredering columns
                            bdc_up_aux = pd.concat([bdc_up_aux, aux], ignore_index=True)  #concat into aux df
                else:
                    pass
            else: # the values are different for each mode           
                #print("CHECK IT bc for bdc_up values should be the same for all modes")
                warning_list.append("CHECK IT bc for bdc_up values should be the same for all modes:  "+str(tn)+' - '+str(nod)) #it should not happen add into warning list
        
            bdc_up_aux = bdc_up_aux[bdc_up_df.columns]     #adjust columns
            bdc_up_df = pd.concat([bdc_up_df, bdc_up_aux], ignore_index = True) #concat into final df

    ####################################################################################################################
    # =============================================================================
    #         RELATION ACTIVITY
    # =============================================================================
            aux_rel = pd.DataFrame(columns=['relation','node_rel','year_rel','node_loc','year_act','mode', 'value'])
            
            
            if 'nan' not in list(set((list(itertools.chain(*BLUES_tec_df.loc[idx[tn, [nod], :], ('act_constraints')].values.tolist()))))): #if there is no NAN
                
                cc_aux = BLUES_tec_df.loc[idx[tn, [nod], :], ('act_constraints')].values.tolist() #get constraints
                #cc1 = unique(list(itertools.chain(*cc1))) #flattening the list 
                
                cv_aux = BLUES_tec_df.loc[idx[tn, [nod], :], ('act_constraint_values')].values.tolist() #get constraint values
                
                #if len(cc1) > 0: #check if len cc is bigger than 1 - it means there is more than 1 mode
                mod_aux = 0
                for cc_, cv_ in zip(cc_aux, cv_aux):

                    mod_aux+=1
                    cc = [c.strip('[]').split(',') for c in cc_] #removing brackects and and split by comma
                    cc = [[c.strip("'") if c[0] == "'" else c.strip(" ").strip("'") for c in cc1] for cc1 in cc][0] #removing empty characters and extra '
                    cc1 = cc[:] #copying it to not harm the original list of relations
                               
                    cv = [[c.strip('[]')] for c in cv_] #removing firts and last brackets
                    cv = [c[0].split('], [') if ']' in c[0] else c for c in cv] #removing mid brackets that split the list of values for each relation       
                    cv = [[re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", c2) for c2 in c] for c in cv][0]   
                    cv = [[float(c2.strip())  for c2 in c1] for c1 in cv]          #converting to value
                    
                    if len(cv) != len(cc): #check if len cc is not equal len cv
                        warning_list.append('check - len cv != len cc - '+tn+"    "+nod)
                    
                    else:       #if len cc ==  len cv 
                        for kv,kc in zip(cv,cc): #looping through each value and relation for each mode
                            #for kc, kv in zip(c1, v1): #looping through values and relations

                            if '-' not in kc:    #check if - is in kc
                                aa = rel_df[(rel_df['relation'] == kc) & (rel_df['node'] == nod)]['unit_group'].values.tolist()
                            else:  #if '-' is in kc
                                aa = rel_df[(rel_df['relation'] == kc.split('-')[0]) & (rel_df['node'] == kc.split('-')[-1])]['unit_group'].values.tolist()
                           
                            
                            if [] in aa or len(aa) == 0 : # if aa is empty
                                not_def_rel.append(list((list((set(itertools.chain(*cc1)))), tn))) #appending relations and tec name that have empty unit group
                            else: # if aa is not empty
                                #if 'activity' not in list(set(aa))[0]: #list(set([a for a in aa[0]])) : # check if unit group does not contain activity
          
                                #adjusting value based on the main input
# =============================================================================
#                                 if len(main_in_val) > 0:
#                                     iv = main_in_val[mod_aux-1] 
#                                     #ov = main_val[mod_aux-1]  
#                                     ef_io = list(np.array([1]) / np.array(iv) )
#                                 else:
#                                     #ef_io = main_val[mod_aux-1]  
# =============================================================================
                                ef_io = list(np.array([1]*len(year_act)) )
                                
                                if '-' in kc: # if - is in relation it means that the relation is related to Brazil region
                                    rc1 = kc.split('-')[0] #getting rleation
                                    n_rel = kc.split('-')[-1] #getting relation location
                                    cond = rel_df[(rel_df['relation'] == kc.split('-')[0]) & (rel_df['node'] == kc.split('-')[-1])]['condition'].values.tolist()[0]
                                    
                                else:
                                    rc1 =  kc
                                    n_rel = nod
                                    cond = rel_df[(rel_df['relation'] == kc) & (rel_df['node'] == nod)]['condition'].values.tolist()[0]
                                
                                if cond == 'i':
                                    ef_io2 = main_in_val[mod_aux-1] 
                                else:
                                    ef_io2 = main_val[mod_aux-1]
                                
                                if len(kv)>count_y: #+1: #if value varies through years
                                    kv = kv[count_y:] # do not get first year if it is after 2010
                                while len(kv) <len(year_act):
                                    kv = kv + [kv[-1]] #adding values to match with years length                                 
                                for y,kvv,ivv,oivv in zip(year_act,kv,ef_io,ef_io2): #looping through years and values
                                    aux = expand_grid_name(['relation','node_rel','year_rel','node_loc','year_act','mode',   'value'], #df columns
                                                           [rc1],       [n_rel],    [y],        [nod],     [y],   [mod_aux], [kvv*ivv*oivv]) #creating dataframe from these values
                               
                                    aux = aux[aux_rel.columns] 
                                    aux_rel = pd.concat([aux_rel, aux], ignore_index=True) #appending dataframe into auxiliary dataframe
    
                aux_rel = aux_rel[relation_activity_df.columns] #ordering columns
                relation_activity_df = pd.concat([relation_activity_df, aux_rel], ignore_index=True) #appending into the main dataframe

    ####################################################################################################################
    # =============================================================================
    #         RELATION CAP
    # =============================================================================
            aux_rel = pd.DataFrame(columns=['relation','node_rel','year_rel','value'])
            
            
            if 'nan' not in list(set((list(itertools.chain(*BLUES_tec_df.loc[idx[tn, [nod], :], ('cap_constraints')].values.tolist()))))): #if there is no NAN
                
                cc_aux = BLUES_tec_df.loc[idx[tn, [nod], :], ('cap_constraints')].values.tolist() #get constraints
                #cc1 = unique(list(itertools.chain(*cc1))) #flattening the list 
                
                cv_aux = BLUES_tec_df.loc[idx[tn, [nod], :], ('cap_constraint_values')].values.tolist() #get constraint values
                
                #if len(cc1) > 0: #check if len cc is bigger than 1 - it means there is more than 1 mode
                mod_aux = 0
                for cc_, cv_ in zip(cc_aux, cv_aux):

                    mod_aux+=1
                    cc = [c.strip('[]').split(',') for c in cc_] #removing brackects and and split by comma
                    cc = [[c.strip("'") if c[0] == "'" else c.strip(" ").strip("'") for c in cc1] for cc1 in cc][0] #removing empty characters and extra '
                    cc1 = cc[:] #copying it to not harm the original list of relations
                               
                    cv = [[c.strip('[]')] for c in cv_] #removing firts and last brackets
                    cv = [c[0].split('], [') if ']' in c[0] else c for c in cv] #removing mid brackets that split the list of values for each relation       
                    cv = [[re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", c2) for c2 in c] for c in cv][0]   
                    cv = [[float(c2.strip())  for c2 in c1] for c1 in cv]          #converting to value
                    
                    if len(cv) != len(cc): #check if len cc is not equal len cv
                        warning_list.append('check - len cv != len cc - '+tn+"    "+nod)
                    
                    else:       #if len cc ==  len cv 
                        for kv,kc in zip(cv,cc): #looping through each value and relation for each mode
                            #for kc, kv in zip(c1, v1): #looping through values and relations

                            if '-' not in kc:    #check if - is in kc
                                aa = rel_df[(rel_df['relation'] == kc) & (rel_df['node'] == nod)]['unit_group'].values.tolist()
                            else:  #if '-' is in kc
                                aa = rel_df[(rel_df['relation'] == kc.split('-')[0]) & (rel_df['node'] == kc.split('-')[-1])]['unit_group'].values.tolist()
                           
                            
                            if [] in aa or len(aa) == 0 : # if aa is empty
                                not_def_rel.append(list((list((set(itertools.chain(*cc1)))), tn))) #appending relations and tec name that have empty unit group
                            else: # if aa is not empty
                                #if 'activity' not in list(set(aa))[0]: #list(set([a for a in aa[0]])) : # check if unit group does not contain activity
          

                                
                                if '-' in kc: # if - is in relation it means that the relation is related to Brazil region
                                    rc1 = kc.split('-')[0] #getting rleation
                                    n_rel = kc.split('-')[-1] #getting relation location
                                    #cond = rel_df[(rel_df['relation'] == kc.split('-')[0]) & (rel_df['node'] == kc.split('-')[-1])]['condition'].values.tolist()[0]
                                    
                                else:
                                    rc1 =  kc
                                    n_rel = nod
                                    #cond = rel_df[(rel_df['relation'] == kc) & (rel_df['node'] == nod)]['condition'].values.tolist()[0]
                                

                                
                                if len(kv)>count_y: #+1: #if value varies through years
                                    kv = kv[count_y:] # do not get first year if it is after 2010
                                while len(kv) <len(year_act):
                                    kv = kv + [kv[-1]] #adding values to match with years length                                 
                                for y,kvv in zip(year_act,kv): #looping through years and values
                                    aux = expand_grid_name(['relation','node_rel','year_rel',   'value'], #df columns
                                                           [rc1],       [n_rel],    [y],        [kvv]) #creating dataframe from these values
                               
                                    aux = aux[aux_rel.columns] 
                                    aux_rel = pd.concat([aux_rel, aux], ignore_index=True) #appending dataframe into auxiliary dataframe
    
                aux_rel = aux_rel[relation_capacity_df.columns] #ordering columns
                relation_capacity_df = pd.concat([relation_capacity_df, aux_rel], ignore_index=True) #appending into the main dataframe
            
            
    ####################################################################################################################    
        if tn1[0].islower():
            tn1 = tn1[0].upper()+tn1[1:]
        tec_list = [        
                             vtg_year_df,
                             vtg_year_time_df,
                             nodes, 
                             year_act, 
                             list(set( vtg_year_time['time'].values.tolist() ) ),
                             vtgs,
                     #types = power_tec_df.loc[tn, 'type'],
                             list(set(mod)),
                             lft_df,
                             main_input_df,
                             inpt_df,
                             output_df,#
                             main_output_df,
                             emission_factor_df,
                             capacity_factor_df,#,
                             ctime_df,
                             inv_cost_df,
                             fix_cost_df,
                             var_cost_df,
                             muf_df,
                             historical_new_capacity_df,
                             #growth_new_capacity_up = growth_new_capacity_up #THERE ARE NO VALUES FOR THIS set(list(itertools.chain(*BLUES_tec_df['growth_new_capacity_up'].values.tolist())))
                             bdc_up_df,
                             bdc_lo_df,
                             bdi_lo_df,
                             bdi_up_df,
                             bda_up_df,
                             bda_lo_df,
                             relation_activity_df,
                             relation_capacity_df
                             ]
        
        return tec_list

    dict_names = ['vtg_year', 'vtg_year_time', 'nodes', 'year_act', 'times', 'vintages', #types = power_tec_df.loc[tn, 'type'],
                 'modes',  'technical_lifetime' ,'main_input', 'input' ,'output' , 'main_output', 'emission_factor','capacity_factor', 'construction_time' ,
                'inv_cost' ,'fix_cost' , 'var_cost' ,'min_utilization_factor' ,'historical_new_capacity' ,
                    #growth_new_capacity_up = growth_new_capacity_up #THERE ARE NO VALUES FOR THIS set(list(itertools.chain(*BLUES_tec_df['growth_new_capacity_up'].values.tolist())))
                'bound_new_capacity_up' ,'bound_new_capacity_lo' ,'bound_total_capacity_lo' ,'bound_total_capacity_up' ,
                'bound_activity_up' ,'bound_activity_lo' ,'relation_activity','relation_capacity'] 
    
    #tt = tecs_f[250:280]+tecs_f[1500:1530]
    #tec_list_test = list(map( aux__build_tec, tt)) #only for testing avoid to change the original value
    #tec_test = {key:{d:td for d,td in zip(dict_names, td1)} for key,td1 in zip(tt,  tec_list_test)  }  #for tests 
      
# =============================================================================
#     #theone below is the right one the above is for tests
# =============================================================================
    
    tec_list_df = list(map( aux__build_tec, tecs_)) #map is applying the building function to each of tecs in the list  
    
    if tec_dict == {}   :
        tec_dict = {key:{d:td for d,td in zip(dict_names, td1)} for key,td1 in zip(tecs_,  tec_list_df)  } #building dictionary
    else:
         tec_dict_aux = {key:{d:td for d,td in zip(dict_names, td1)} for key,td1 in zip(tecs_,  tec_list_df)  } #building dictionary
         tec_dict = {**tec_dict, **tec_dict_aux}       
    
    warning_list = list(set(warning_list))   
    
    dt_final =  datetime.now() - startTime
    print('Done.')
    print('Runtime of Building dictionary of technologies parameters:    '+ str(dt_final))
    return warning_list, tec_dict, not_def_rel

#%%

# =============================================================================

# Check constraint that might be Resource (not now)
# Talk to Roberto --> future, buy gams, decisions, Alex's ideas

#______________________________________________________________________________________
#Future Tasks:
#Time to 288 --> JAN 1 JAN 2... JAN 24 FEV 1...FEV24 MAR 1...
#Storage -->Like Newave
#Battery??
#Pumped Hydro??
#Spatial Resolution for power system??
#Industry??
# =============================================================================