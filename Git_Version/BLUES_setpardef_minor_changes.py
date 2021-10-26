# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:44:17 2018

@author: amendolf
"""


ldr_original = ldr[-1]
#%%

#aDDING NEW DEFICIT FOR ind_electricity Industrial Metallurgy demand
def_tecs = [t for t in tecs_f2 if'Deficit' in t]

#%%
if BLUES:
    def_tecs = [vd for vd in def_tecs if vd not in ['Deficit_Gq', 'Deficit_Pq', 'Deficit_Cq', 'Deficit_Wq']]
    tecs_f2_aux = tecs_f2[:]
    tecs_f2 = [v for v in tecs_f2 if v not in def_tecs]
for vt in def_tecs:
    del tec_dict[vt] #REMOVE
#%%
#ADDING muf = 1  to secondary land tecs
for st in sec_land_tec:
    aux1 = pd.DataFrame()
    for nnod in tec_dict[st]['nodes']:
        #print(nnod)
        aux = expand_grid_name(['node', 'vintage', 'year_act', 'value'], [nnod], tec_dict[st]['vtg_year'][tec_dict[st]['vtg_year']['node'] == nnod][['vintage', 'year_act']].values, [1] )
        aux1 = pd.concat([aux1, aux], ignore_index=True)
    tec_dict[st]['min_utilization_factor'] = aux1
#%%
##	node	world - regions - countries - grid cells
node_set = adb_files


    
#%%
## map_node

map_node_set = pd.DataFrame({'node1':adb_files+[main_region]*len([a for a in adb_files if a != main_region]), 
                                'node2': adb_files+[a for a in adb_files if a != main_region]})

#%%
## node categoriea cat_node
cat_node_set = pd.DataFrame( columns=['type_node', 'node'])
## type_node
type_node_set = list(set( cat_node_set.type_node ))
#%%
##	year_act		years (over entire model horizon)  
year_set = year_all[:]
##	year(year_act)		years included in a model instance (for myopic or rolling-horizon optimization)	
#year_set = [y for y in year_act if y> baseyear]


#%%
##	type_year		types of year aggregations 
type_year_set = year_act + ['firstmodelyear']  + ['cumulative'] #++ ['lastmodelyear']


#%%
##  cat_year(type_year,year_act)  	mapping of years to respective categories					
cat_year_set = pd.DataFrame({'type_year': type_year_set + [type_year_set[-1]]* (len(year_act)-1) , 
                             'year_act': year_act + [year_act[0]] + year_act #+ [year_act[-1]]
                             })
#%%

    
if has_ldr:
##	time		subannual time periods (seasons - days - hours) 	
    parent_time = 	'year'	
    time_set =  [parent_time] + time # Need to add year to allow for output on annual basis
    
    ## map_time
    map_time_set = pd.DataFrame({'time': ['year'] + ['year']*len(time) + time , 
                                 'time2': ['year'] + time + time})
else:
    time_set = time.copy()
    map_time_set = pd.DataFrame({'time': ['year'],
                                 'time2': ['year'] })
#%%							
##	lvl_spatial hierarchical levels of spatial resolution 
#lvl_spatial_set = ['bcu', 'national'] #keep it it might be helpful

lvl_spatial_set = ['national']
## map_spatial_hierarchy
map_spatial_hierarchy_map =  pd.DataFrame({'lvl_spatial': lvl_spatial_set[0],
                                           'node': [a for a in adb_files if a != main_region],
                                           'node_parent': main_region }) 


#%%	
if has_ldr:
    ##	lvl_temporal	hierarchical levels of temporal resolution 
    lvl_temporal_set = ['year', 'sub_anual']
    
    ## map_temporal_hierarchy
    
    map_temporal_hierarchy_map = pd.DataFrame({'lvl_temporal':[lvl_temporal_set[0]]+[lvl_temporal_set[1]]*len(time),
                                               'time': ['year']+time,
                                               'time2': 'year'})
				
#%%
##	rating		identifies the 'quality' of the renewable energy potential (bins acc. to Sullivan)     
rating_set = [ 1 ]

#%%
# =============================================================================
time_aux = { v:list(range(1+24*(v-1), 1+24*(v-1)+24)) for v in list(range(1,7))}
new_keys = list(reversed([1,3,5,7,9,11]))
# 
old_k = list(reversed(list(time_aux.keys())))
# 
# 
for key,n_key in zip(old_k, new_keys):
    time_aux[n_key] = time_aux.pop(key)
# =============================================================================
    
new_time_duration = [round(t/365, 6) for t in [31,28,31,30,31,30,31,31,30,31,30,31]]


#%%   
## 	duration_period(year_act)		duration of one multi-year period (in years)
year_reverse = list(reversed(year_set))
durat_par = [year_reverse[y] - year_reverse[y+1] for y in range(len(year_reverse)-1)]
if period_interp == 'recursive':
    durat_par = durat_par + [durat_par[-1]] #original recursive tme period interpretation
elif period_interp == 'progressive':
    durat_par = [durat_par[0]] + durat_par # for forwarding year instead of recursive
duration_period_par = pd.DataFrame({'year': year_set,
                                    'value': list(reversed(durat_par)) })

if BLUES:
    emission_duration_period_par = duration_period_par.copy()
    emission_duration_period_par['value'] = np.where(emission_duration_period_par['year'] == year_set[1], 5, emission_duration_period_par['value'])

if has_ldr:
    #for lts in lvl_temporal_set: make it flexible for more than sub-annual times
    ## 	duration_time(time)		duration of one time slice (relative to 1)
    #time_duration = np.array(time_duration) * (1/sum(time_duration))
    if sum(time_duration) != 1:
        print("time duration sum is different from 1")
        if round(sum(time_duration),3) == 1:
            print("... but it is being adjusted auttomatically because it is very close to 1")
            #time_duration[-1] = round(1 - sum(time_duration[:-1]),6)
            #time_duration2[-1] = round(1 - sum(time_duration2[:-1]),4)
            time_duration[-1] = float_round(1 - sum(time_duration[:-1]), 7, ceil) 
            time_duration = [round(v,7) for v in time_duration]
            if round(sum(time_duration),5) == 1:
                print("--> adjusted\n")
            else:
                print("\ncheck why duration time sum is not 1")
                sys.exit()        
        else:
            print("\ncheck why duration time sum is not 1")
            sys.exit()
        
    #time_duration = time_duration.tolist()  
    duration_time_par = pd.DataFrame({'time': ['year'] + time,
                                     'value':[1.0] + time_duration})
    
    m = len(time) 
    drate_time = ((1 + drate)**(1/m)) - 1
else:
        duration_time_par = pd.DataFrame({'time': ['year'],
                                     'value':[1.0]})
    
## 	interestrate(year_act)		interest rate (to compute discount factor)
interestrate_par = pd.DataFrame({'year': year_act,
                                 'value':[drate] * len(year_act)}) #drate #check if it is a annual rate


#%%
## 	demand(node,commodity,level,year_act,time) 		exogenous demand levels
#make script get demand
demand_par = demand_df[['node', 'commodity', 'level', 'year_act', 'time', 'value']].copy() 





#%%
##	tec		technologies
# Create list of technology parameters				
params_list = tec_dict.copy()
tec_set = sorted(list(tec_dict.keys()))

inv_tec_set = tec_set #


#%%
# get parameter names and remove redundancies
variables = list(set(list(itertools.chain(*[list(tec_dict[t].keys()) for t in tec_set]))))
variables = [v for v in variables if v not in ['nodes', 'modes', 'year_act', 'times', 'vintages', 'vtg_year', 'vtg_year_time', 'types']]

#%%
#len_tec = pd.DataFrame(columns = sorted(variables), index=tec_set)  #length of each variable by tech
def aux_f1(tt):  
    #print(tt)
    if v in params_list[tt].keys():
        df_aux = params_list[tt][v]
        #len_tec[v].loc[[tt]] = len(df_aux)
        if len(df_aux) > 0:
            df_aux['tec'] = tt
            df_aux = df_aux[[df_aux.columns[0]]+['tec']+list(df_aux.columns[1:-1])]
        else:
            df_aux = pd.DataFrame( columns=df_aux.columns.tolist() + ['tec'] )

    return df_aux

# =============================================================================
# v ='emission_factor'
# df1[(df1['value'] < 0.000001) & (df1['value'] > 0.0)]
# AD=df1[df1['level'].isin(['Water_block', 'Water_Treated', 'Raw_Water'])]
# =============================================================================
gb = {}
var_dict = {}  #variable dictionary
for v in variables: #looping through defined variables
    print('VARIABLE: '+v)
    df = list(map( aux_f1, tec_set)) #map is applying the building function to each of tecs in the list   
    df1 = pd.concat(df)
    df1.reset_index(drop=True, inplace=True)
    #think of a way to adjust columns
    #print(df1[df1['value']>0]['value'].min())


    df1['value'] = df1['value'].astype(float)
    ycol = [ll for ll in df1.columns if 'vintage' in ll or 'year' in ll]
    for ycol1 in ycol:
        df1[ycol1] = df1[ycol1].astype(int)
    #df1['value'] = np.where(np.logical_and( df1['value']<0.000001 , df1['value']>0), 0.0, df1['value']) #TRY WITH 0.0001
    

    
    var_dict[v+'_par'] = df1    
    
[ v for v in df1['value'] if type(v) == str ]
print('finished')
locals().update(var_dict)

#Compare to the R version
#a = [te for te in tec_set if te.startswith('irr_')]
#list(set(relation_activity_par['relation'].values.tolist() ) )
upper = [vu for vu in list(set(relation_activity_par['relation'].values.tolist() ) ) if vu.endswith('U')]
lower = [vl for vl in list(set(relation_activity_par['relation'].values.tolist() ) ) if vl.endswith('L')]

# =============================================================================
#     attempts to reduce infeasibility
# =============================================================================

relation_activity_par['value'] = np.where( (relation_activity_par['relation'].isin(upper) ) & (relation_activity_par['year_rel']==2010 ),  relation_activity_par['value'].apply(lambda x: float_round(x, 8, ceil)), relation_activity_par['value'])
relation_activity_par['value'] = np.where( (relation_activity_par['relation'].isin(lower) ) & (relation_activity_par['year_rel']==2010 ),  relation_activity_par['value'].apply(lambda x: float_round(x, 8, floor)), relation_activity_par['value'])


relation_activity_par['value'] = np.where( (relation_activity_par['relation'].isin(upper) ) & (relation_activity_par['year_rel']!=2010 ),  relation_activity_par['value'].apply(lambda x: float_round(x, 6, ceil)), relation_activity_par['value'])
relation_activity_par['value'] = np.where( (relation_activity_par['relation'].isin(lower) ) & (relation_activity_par['year_rel']!=2010 ),  relation_activity_par['value'].apply(lambda x: float_round(x, 6, floor)), relation_activity_par['value'])






input_par.drop_duplicates(keep='first', inplace=True, ignore_index=True)




#%%
relation_total_capacity_par = relation_capacity_par.copy()


#%% ##	emission		greenhouse gases - pollutants - etc.	

emission_set = list(set(list(itertools.chain(*[list(params_list[tt]['emission_factor'].emission) for tt in tec_set if 'emission_factor' in params_list[tt].keys()]))))

#LOOKING FOR emissions WITH SAME LETTER
emission_set_aux = [rr.lower() for rr in emission_set]

EE_set = []
for ee in emission_set:
    lee = len([raux for raux in emission_set_aux if raux==ee.lower()])
    if lee >1:
        EE_set.append(ee)
        
if len(EE_set) >0:
    print("adjust emission set __ there are duplicated values")
    sys.exit()

##	type_emission		types of emission aggregations
# =============================================================================
# #ADJUST TYPE EMISSION --> CO2 N2O CH4 GHG
# =============================================================================
type_emission_set = [ll[:3] if ll not in ['TCO2', 'cGHG'] else ll for ll in  emission_set]


##	cat_emission(type_emission,emission)		mapping of emissions to respective categories					
cat_emission_map = pd.DataFrame({'type_emission': type_emission_set,
                                 'emission': emission_set})					

tp_emi = dict(zip(emission_set, type_emission_set ))

emission_cost['type_emission'] = emission_cost['type_emission'].map(tp_emi)
emission_bound_par = emission_bound.copy()
emission_cost.drop_duplicates(keep='first', inplace=True, ignore_index=True)
emission_cost_par = emission_cost.copy()


emission_bound_par = emission_bound_par[emission_bound_par['value'] < 9999999999]
emission_bound_par.drop_duplicates(keep='first', inplace=True, ignore_index=True)

#%%

# =============================================================================
# #CONVERTING CUMULATIVE RELATION INTO RESOURCE
# =============================================================================
relation_upper_par = relation_upper.copy()
relation_lower_par = relation_lower.copy()
relation_cost_par = relation_cost.copy()

#relation_activity_par[relation_activity_par['year_rel'] == 'cumulative']['relation'].values.tolist()
relation_upper_par[relation_upper_par['year_rel'] == 'cumulative']['relation'].values.tolist()
CUM_RELATION = list(set( relation_upper_par[relation_upper_par['year_rel'] == 'cumulative']['relation'].values.tolist() ))

rel2res_upper = pd.DataFrame(columns=relation_upper_par.columns) 
rel2res_lower = pd.DataFrame(columns=relation_lower_par.columns)     
rel2res_cost = pd.DataFrame(columns=relation_cost_par.columns) 
rel2_res_act = pd.DataFrame(columns=relation_activity_par.columns)
    
for cr in CUM_RELATION:
    rel2res_upper = pd.concat([rel2res_upper, relation_upper_par[relation_upper_par['relation'] == cr] ], ignore_index=True )
    rel2res_lower = pd.concat([rel2res_lower, relation_lower_par[relation_lower_par['relation'] == cr]       ], ignore_index=True )
    rel2res_cost = pd.concat([rel2res_cost, relation_cost_par[relation_cost_par['relation'] == cr]   ], ignore_index=True )
    rel2_res_act = pd.concat([rel2_res_act, relation_activity_par[relation_activity_par['relation'] == cr] ], ignore_index=True )
     
    relation_upper_par = relation_upper_par[relation_upper_par['relation'] != cr]    
    relation_lower_par = relation_lower_par[relation_lower_par['relation'] != cr]    
    relation_cost_par = relation_cost_par[relation_cost_par['relation'] != cr]       
    relation_activity_par = relation_activity_par[relation_activity_par['relation'] != cr]
    

rel2res_upper.rename(columns={'relation': 'commodity', 'node_rel':'node', 'year_rel':'grade' }, inplace=True)
rel2res_upper['grade'] = 'a'
rel2res_upper = rel2res_upper[['node', 'commodity', 'grade', 'value']]
#resource_cost(node,commodity,grade,year_all) 
rel2res_cost.rename(columns={'relation': 'commodity', 'node_rel':'node', 'year_rel':'year_act' }, inplace=True)
rel2res_cost['grade'] = 'a'
rel2res_cost = rel2res_cost[['node','commodity','grade', 'year_act', 'value']]


for t in list(set( rel2_res_act['tec'].values.tolist() )):
    input_par_rel2res = input_par[input_par['tec'] == t] #input par of the current technology
    
    input_par_rel2res_aux = input_par_rel2res.copy()
    input_par_rel2res_aux['commodity'] = list(set( rel2_res_act[rel2_res_act['tec'] == t]['relation'].values.tolist() ))[0]
    input_par_rel2res_aux['level'] = 'Resources'
    input_par_rel2res_aux['value'] = list(set( rel2_res_act[rel2_res_act['tec'] == t]['value'].values.tolist() ))[0]
    input_par_rel2res_aux.drop_duplicates(keep='first', inplace=True, ignore_index=True)
    
    input_par_rel2res = pd.concat([input_par_rel2res, input_par_rel2res_aux], ignore_index=True)
    
    input_par = input_par[input_par['tec'] != t]
    input_par = pd.concat([input_par, input_par_rel2res], ignore_index=True)
    
#%%
#RELATIONS
# =============================================================================
# THERE ARE MANY RELATIONS THAT CAN BE REPRESENTED BY OTHER PARAMETERS
# =============================================================================
    
relation_set = sorted(list(set(relation_activity_par['relation'].values.tolist() ) ) + list(set(relation_upper_par['relation'].values.tolist() ) )  + list(set(relation_lower_par['relation'].values.tolist() ) ) + list(set(relation_cost_par['relation'].values.tolist() ) ))

relation_set = sorted(list(set(relation_set)) )

original_relation_set = relation_set[:]

# =============================================================================
# #LOOKING FOR RELATIONS WITH SAME LETTER
# =============================================================================
relation_set_aux = [rr.lower() for rr in relation_set]
dupe = [item for item, count in collections.Counter(relation_set_aux).items() if count > 1]
for RR in dupe:
    idx = [i for i,val in enumerate(relation_set_aux ) if val==RR]
    for RR2 in idx:
        if RR2 == min(idx):
            relation_set[RR2] = relation_set[RR2]+'_A'
        else:
            relation_set[RR2] = relation_set[RR2]+'_B'
            

relation_set = sorted(relation_set) 
        
new_rel_dict = dict(zip(original_relation_set, relation_set ))

relation_upper_par['relation'] = [new_rel_dict[rrr] for rrr in relation_upper_par['relation']]
relation_lower_par['relation'] = [new_rel_dict[rrr] for rrr in relation_lower_par['relation']]
relation_cost_par['relation'] = [new_rel_dict[rrr] for rrr in relation_cost_par['relation']]
relation_activity_par['relation'] = [new_rel_dict[rrr] for rrr in relation_activity_par['relation']]
#main_rel_par['relation'] = [new_rel_dict[rrr] for rrr in main_rel_par['relation']]


#%%
##	grade		grades of extraction of raw materials    				
grade_set =  ['a']
resource_volume_par = resource_df.copy()#(node,commodity,grade)               volume of resources in-situ at start of the model horizon


resource_cost = resource_volume_par.copy() #(node,commodity,grade,year_act)      
resource_cost['value'] = 0.0 # -----------------> make it as it is #Try without this line
resource_cost_par = expand_grid_name(resource_cost.columns.tolist()+['year_act'], resource_cost.values, year_act)


resource_volume_par = pd.concat([resource_volume_par, rel2res_upper], ignore_index=True)

resource_cost_par = resource_cost_par[['node', 'commodity', 'grade', 'year_act' ,'value']]
resource_cost_par = pd.concat([resource_cost_par, rel2res_cost], ignore_index=True)
 
level_options = list(set(input_par[input_par['commodity'].isin(resource_volume_par['commodity'].values.tolist())]['level'].values.tolist() )) #    level_resource (level)                  subset of 'level' to mark all levels related to make hfossil resources
level_resource_par = ['Resources']
#level_renewable(level)                  subset of 'level' to mark all levels related to renewable resources
    
#%%
lvl_com_cm_df = lvl_com_check(tpp, adb_files, ldb) #level commodity pairs --> original value taken from oriignal database
    
##	commodity		resources - electricity - water - land availability - etc.
commodity_set = list(set(input_par.commodity.values.tolist()+output_par.commodity.values.tolist()+demand_par.commodity.values.tolist()))

##	level		levels of the reference energy system or supply chain ( primary - secondary - ... - useful )			
level_set = list(set(input_par.level.values.tolist()+output_par.level.values.tolist()+demand_par.level.values.tolist()))

#getting all commodity level pairs
#list(set([(c,l) for c,l in zip(input_par.commodity.values.tolist(), input_par.level.values.tolist())]))

com_lvl_set = list(set([(c,l) for c,l in zip(input_par.commodity.values.tolist(), input_par.level.values.tolist())] 
                + list(set([(c,l) for c,l in zip(output_par.commodity.values.tolist(), output_par.level.values.tolist())])) 
                + list(set([(c,l) for c,l in zip(demand_par.commodity.values.tolist(), demand_par.level.values.tolist())])) ) )

com_lvl_df = pd.DataFrame(com_lvl_set, columns=['commodity', 'level'])   
com_lvl_df = com_lvl_df[['level', 'commodity']]   
com_lvl_df.sort_values(['level', 'commodity'], inplace=True)        
com_lvl_df.reset_index(inplace=True, drop=True)  
com_lvl_df.to_csv(os.path.join(path_case,project,ldb,project+'_'+task +'_commodity_level.csv'), index=False, sep=';')

# =============================================================================
print("Building lvl_com_dict")

lvl_com_list = [(l,c) for l,c in zip(com_lvl_df['level'].values.tolist(), com_lvl_df['commodity'].values.tolist())]
lvl_com_t_lst = list(lvl_com_t.keys())
new_lc = [lc for lc in lvl_com_list if lc not in lvl_com_t_lst]
for lvl,com in new_lc:
#def aux_lc_f(lst):
    #lvl, com = lst
    lvl_com_t[(lvl,com)] ={}   #['tec_input'] = {}
    lvl_com_t[(lvl,com)] ={}   #['tec_outut'] = {}
#     
    inp_t = list(set( input_par[(input_par['level'] == lvl) & (input_par['commodity'] == com)]['tec'].values.tolist() ))
    #inp_t = list(set(inp_t))
    lvl_com_t[(lvl,com)]['tec_input'] = inp_t
# 
    out_t = list(set( output_par[(output_par['level'] == lvl) & (output_par['commodity'] == com)]['tec'].values.tolist() ))
    #out_t = list(set(out_t))
    lvl_com_t[(lvl,com)]['tec_outut'] = out_t
    #return com_lvl_dic[lst]
# 
# aux_lst = list(map(aux_lc_f, lvl_com_list))
#com_lvl_dic = reduce(aux_lc_f, lvl_com_list[:3])
# del aux_lst
# =============================================================================
   
a = com_lvl_df[['level', 'commodity']].values.tolist()
b = lvl_com_cm_df[['level', 'commodity']].values.tolist()

#dd = tec_dict['Port_Impcoal_Reg_Co']

#comparing original com-lvl pairs with the ones generated in the conversion
lvl_com_cm_df['status2'] = ['Used' if bb in a else 'Not_Used' for bb in b]
lvl_com_cm_df['status3'] = np.where(lvl_com_cm_df['status2'] == lvl_com_cm_df['status'], True, False)
#lvl_com_cm_df[lvl_com_cm_df['status3'] == False]
#%%
#LOOKING FOR commodities WITH SAME LETTER
commodity_set_aux = [rr.lower() for rr in commodity_set]
CC_set = []
for cc in commodity_set:
    lcc = len([raux for raux in commodity_set_aux if raux==cc.lower()])
    if lcc >1:
        CC_set.append(cc)

if len(CC_set) >0:
    print("adjust commodity set __ there are duplicated values")
    sys.exit()
        
#LOOKING FOR commodities WITH SAME LETTER
level_set_aux = [rr.lower() for rr in level_set]
LL_set = []
for ll in level_set:
    lll = len([raux for raux in level_set_aux if raux==ll.lower()])
    if lll >1:
        LL_set.append(ll)
        
if len(LL_set) >0:
    print("adjust level set __ there are duplicated values")
    sys.exit()
#%%
#commodity level balance equality
if BLUES:
    land_primary = [l for l in land_tec if l not in sec_land_tec]
    land_bal_eq = land_generation + land_primary
    
    com_eq = []
    lvl_eq = []
    for ll in land_bal_eq:
        com_eq.append( list(set(tec_dict[ll]['output']['commodity'].values.tolist()) ) )
        lvl_eq.append( list(set(tec_dict[ll]['output']['level'].values.tolist()) ) )    
     
    com_eq = list(itertools.chain(*com_eq))   
    lvl_eq = list(itertools.chain(*lvl_eq))     
    land_bl_eq_df = pd.DataFrame({'level': lvl_eq,
                                  'commodity': com_eq 
                                  })
    
    bal_eq_df = pd.concat([balance_eq_df, land_bl_eq_df], ignore_index=True)
    bal_eq_df.drop_duplicates(keep='first', inplace=True, ignore_index=True) # in old pandas version ignore_index is not an option for drop duplicate

    bal_eq_cum_rel = pd.DataFrame({'level':  ['Resources']*len(CUM_RELATION ),
                                      'commodity': CUM_RELATION 
                                      })
    bal_eq_df = pd.concat([bal_eq_df, bal_eq_cum_rel], ignore_index=True)

    co2_com = ['captured_co2', 'transported_co2']
    bal_eq_co2 = pd.DataFrame({'level':  ['Auxiliary'] *len(co2_com),
                                      'commodity': co2_com 
                                      })
    
    bal_eq_df = pd.concat([bal_eq_df, bal_eq_co2], ignore_index=True)
    
else:
    bal_eq_df = balance_eq_df.copy()
    bal_eq_df.drop_duplicates(keep='first', inplace=True, ignore_index=True) # in old pandas version ignore_index is not an option for drop duplicate


bal_eq_df_map = bal_eq_df[bal_eq_df['commodity'].isin(commodity_set)]
#%%
##	type_tec		types of technologies			
type_tec_set = ['all']

##	cat_tec(type_tec,tec)		mapping of technologies to respective categories
cat_tec_set = pd.DataFrame({'type_tec': type_tec_set * len(tec_set),
                            'tec': tec_set})


##  inv_tec(tec)		technologies that have explicit investment and capacity decision variables
inv_tec_set = list(set(inv_cost_par['tec'].tolist()))
 
##	mode		modes of operation
mode_set = list(set(input_par['mode'].values.tolist()+output_par['mode'].values.tolist()))

#%%
# =============================================================================
# #SET PARAMETERS RELATED TO HISTORICAL TO GUARANTEE LAND BALANCE IN THE FIRST YEAR
# =============================================================================
if BLUES: 
    
    
    dct_tec = dict(zip(land_primary, sec_land_tec))
    if not all([True if dct_tec[k].startswith(k.split('Land_')[-1])  else False for k in list(dct_tec.keys())]):
        raise Exception("land primary tec is not matching land secondary tec")
        
 
    class_land_tec = ['class_'+clt for clt in class_land_tec]
    
    cat_land_tec_df['type_land_tec'] = ['class_'+clt for clt in cat_land_tec_df['type_land_tec'].values.tolist() ]
    cat_land_tec_df.rename(columns={'type_land_tec':'BLUES_type_land_tec', 'technology':'tec'}, inplace=True)
    
    map_land_rel_df['type_land_tec'] = ['class_'+clt for clt in map_land_rel_df['type_land_tec'].values.tolist() ]
    map_land_rel_df.rename(columns={'type_land_tec':'BLUES_type_land_tec', 'technology':'tec'}, inplace=True)    
    

#%%
# adjusting lifetime if tec has hisc bdc or bdi
bdc_t = bound_new_capacity_up_par['tec'].values.tolist() + bound_new_capacity_lo_par['tec'].values.tolist()
bdi_t = bound_total_capacity_up_par['tec'].values.tolist() + bound_total_capacity_lo_par['tec'].values.tolist()
hisc_t = historical_new_capacity_par['tec'].values.tolist()

cap_t = bdc_t + bdi_t+ hisc_t
cap_t = list(set(cap_t))
missing_pll = [t for t in cap_t if t not in technical_lifetime_par['tec'].values.tolist()]
if len(missing_pll) >0:
    print('there are tecs missing pll info')
    sys.exit()

#no_pll_tec = [t for t in tec_set if t not in technical_lifetime_par['tec'].values.tolist()]   

#%%
# =============================================================================
#     FINE ADJUSTMENTS - REVIEW
# =============================================================================
if BLUES:
    #removing unnecessary bounds --> it keeps bda for no_pll tecs
    #bound_activity_up_par1 = bound_activity_up_par[(~bound_activity_up_par['tec'].isin(no_pll_tec) )]
    #bound_activity_up_par2 = bound_activity_up_par[(bound_activity_up_par['tec'].isin(no_pll_tec) )]
    #bound_activity_up_par = bound_activity_up_par[(bound_activity_up_par['value'] < 990000) | (bound_activity_up_par['tec'].isin(no_pll_tec) )]
    
    #bound_activity_up_par['value'] = np.where(bound_activity_up_par['value'] >990000, 9999999, bound_activity_up_par['value'])
    #list(set(bound_activity_up_par[(bound_activity_up_par['value'] > 99999999) ]['tec'].tolist() ) ) 999999999
    #bound_activity_up_par = bound_activity_up_par[(bound_activity_up_par['value'] <= 99999999) ]	#99999999
    bound_activity_lo_par = bound_activity_lo_par[bound_activity_lo_par['value'] > 0]		
    
    #bound_new_capacity_up_par = bound_new_capacity_up_par[bound_new_capacity_up_par['value'] <= 99999999] #9999999
    #bound_new_capacity_up_par[bound_new_capacity_up_par['value'] >= 9999999]
    #bound_new_capacity_up_par['value'] = np.where(bound_new_capacity_up_par['value'] >990000, 990000, bound_new_capacity_up_par['value'])
    bound_new_capacity_lo_par = bound_new_capacity_lo_par[bound_new_capacity_lo_par['value'] > 0]
    
    #bound_total_capacity_up_par = bound_total_capacity_up_par[bound_total_capacity_up_par['value'] < 99999999] #9999999
    #bound_total_capacity_up_par[bound_total_capacity_up_par['value'] >= 9999999]
    #bound_total_capacity_up_par['value'] = np.where(bound_total_capacity_up_par['value'] >990000, 990000, bound_total_capacity_up_par['value'])
    bound_total_capacity_lo_par = bound_total_capacity_lo_par[bound_total_capacity_lo_par['value'] > 0]
    
    construction_time_par = construction_time_par[construction_time_par['value'] > 0]
    construction_time_par['value'] = np.where(np.logical_and(construction_time_par['value'] > 0, construction_time_par['value'] < 1), 1, construction_time_par['value'])
    
    historical_new_capacity_par = historical_new_capacity_par[historical_new_capacity_par['value'] >= 0]
    
    #relation_lower_par = relation_lower_par[relation_lower_par['value'] > -9999999]
    #relation_upper_par = relation_upper_par[relation_upper_par['value'] < 9999999] #9999999
    resource_cost_par = resource_cost_par[resource_cost_par['value'] != 0]
    emission_factor_par = emission_factor_par[emission_factor_par['value'] != 0]
    emission_cost_par = emission_cost_par[emission_cost_par['value'] != 0]
    #emission_bound_par = emission_bound_par[emission_bound_par['value'] < 9999999] #9999999
    

    #var_cost_par = var_cost_par[var_cost_par['value'] == 0 ]
    #inv_cost_par = inv_cost_par[inv_cost_par['value'] == 0 ]
    
    output_par = output_par[output_par['value'] != 0]
    input_par = input_par[input_par['value'] != 0]    

    var_cost_par = var_cost_par[var_cost_par['value'] != 0 ]
    
    fix_cost_par = fix_cost_par[fix_cost_par['value'] > 0 ]
    
    fix_cost_par['value'] = np.where(fix_cost_par['value'] >=1000, round(fix_cost_par['value'], 2), fix_cost_par['value'])
    var_cost_par['value'] = np.where(var_cost_par['value'] >=1000, round(var_cost_par['value'], 2), var_cost_par['value'])



    
    inv_cost_par['value'] = np.where(inv_cost_par['value'] >=1000, round(inv_cost_par['value'], 2), inv_cost_par['value'])
    #inv_cost_par['value'] = np.where(inv_cost_par['value'] == 0, 0.01, inv_cost_par['value'])
    #bound_total_capacity_up_par = bound_total_capacity_up_par[bound_total_capacity_up_par['value'] < 990000]

    capacity_factor_par['value'] =  [ float_round(vv, 2, ceil) for vv in capacity_factor_par['value']]
    min_utilization_factor_par['value'] =  [ float_round(vv, 2, floor) for vv in min_utilization_factor_par['value'] ]






        

#%%

if BLUES:

    #auxiliary relation parametres and sets
    map_relation_par = relation_activity_par[['relation','node_rel', 'year_act']].copy() # (relation, node_rel, year_act)        
    is_relation_upper_par = relation_upper_par[['relation','node_rel', 'year_rel']].copy()
    is_relation_lower_par = relation_lower_par[['relation','node_rel', 'year_rel']].copy()    
    
    def_tec = [t for t in tec_set if t.startswith('Deficit')]
    
    bound_activity_lo_par = bound_activity_lo_par[~bound_activity_lo_par['tec'].isin(def_tec)]

#%%
# =============================================================================
# INCREASING MUF VALUE FOR NUCLEAR TO REDUCE VARIATIONS
# =============================================================================

for tt in ldr_original:
    if tt.startswith('Nuclear_'):
        #print(tt)
        new_v = capacity_factor_par[capacity_factor_par['tec'] == tt]['value'].min() -0.1
        if new_v <0.7:
            new_v = 0.7
        #print(new_v)
                
        min_utilization_factor_par['value'] =  np.where( (min_utilization_factor_par['tec'].str.startswith('Nuclear')) &        
                                                      (min_utilization_factor_par['tec'].isin(ldr_original)),
                                                      new_v,
                                                      min_utilization_factor_par['value'] )
        
#%%
# =============================================================================
# #adding MIN UTIL TIME FACTOR FOR LDR and RELAXING ELECT FINAL - CF=MUF
# =============================================================================

#adding MIN UTIL TIME FACTOR FOR LDR only for Elect_Final
if has_ldr == False:
    min_utilization_time_factor_par = pd.DataFrame(columns=[ 'node', 'tec', 'vintage', 'year_act', 'time', 'value'])
    
    
    aux = capacity_factor_par[capacity_factor_par['tec'] == 'Elect_Final'].copy()
    aux['value'] = aux['value'].values - 0.01
    aux = aux[['node', 'tec', 'vintage', 'year_act', 'time', 'value']]
    min_utilization_time_factor_par = pd.concat([min_utilization_time_factor_par, aux], ignore_index=True)
    
    min_utilization_factor_par = min_utilization_factor_par[min_utilization_factor_par['tec'] != 'Elect_Final']  
    min_utilization_time_factor_par = min_utilization_time_factor_par[[ 'node', 'tec', 'vintage', 'year_act', 'time', 'value']]


    


 #%%       
#ADJUSTING COLUMNS
capacity_factor_par = capacity_factor_par[['node', 'tec', 'vintage', 'year_act', 'time', 'value']]
construction_time_par = construction_time_par[['node', 'tec', 'vintage', 'value']]
technical_lifetime_par = technical_lifetime_par[['node', 'tec', 'vintage', 'value']]
historical_new_capacity_par = historical_new_capacity_par[['node', 'tec',  'vintage', 'value']]
min_utilization_factor_par = min_utilization_factor_par[[ 'node', 'tec', 'vintage', 'year_act', 'value']]

input_par = input_par[['node', 'tec', 'vintage', 'year_act', 'mode', 'node_in', 'commodity', 'level', 'time', 'time_in' , 'value']]
output_par = output_par[['node', 'tec', 'vintage', 'year_act', 'mode', 'node_out', 'commodity', 'level', 'time', 'time_out' , 'value']]

var_cost_par = var_cost_par[['node', 'tec', 'vintage', 'year_act', 'mode', 'time', 'value']]        
inv_cost_par = inv_cost_par[['node', 'tec', 'vintage', 'value']]           
fix_cost_par = fix_cost_par[['node', 'tec', 'vintage', 'year_act', 'value']]   
           
bound_activity_lo_par = bound_activity_lo_par[['node', 'tec', 'year_act', 'mode', 'time', 'value']]
bound_activity_up_par = bound_activity_up_par[['node', 'tec', 'year_act', 'mode', 'time', 'value']]
bound_total_capacity_up_par = bound_total_capacity_up_par[['node', 'tec', 'year_act', 'value']]
bound_total_capacity_lo_par = bound_total_capacity_lo_par[['node', 'tec', 'year_act', 'value']]
bound_new_capacity_up_par = bound_new_capacity_up_par[['node', 'tec', 'vintage', 'value']]
bound_new_capacity_lo_par = bound_new_capacity_lo_par[['node', 'tec', 'vintage', 'value']]

resource_volume_par = resource_volume_par[['node', 'commodity', 'grade', 'value']]
resource_cost_par = resource_cost_par[['node', 'commodity', 'grade', 'year_act', 'value']]

emission_factor_par = emission_factor_par[['node', 'tec', 'vintage', 'year_act', 'mode', 'emission', 'value']]
emission_cost_par = emission_cost_par[['node', 'type_emission', 'type_tec', 'year_act', 'value']] #adjust this 
emission_bound_par = emission_bound_par[['node', 'type_emission', 'type_tec', 'type_year', 'value']]
relation_activity_par = relation_activity_par[['relation', 'node_rel', 'year_rel', 'node_loc', 'tec', 'year_act', 'mode', 'value']]
relation_total_capacity_par = relation_total_capacity_par[['relation', 'node_rel', 'year_rel', 'tec', 'value']]

#bound_activity_up_par[bound_activity_up_par['time'] != 'year']
#%%
emission_factor_par.to_csv(os.path.join(path_case,project,'IX', ldb,project+'_'+task+ldb_aux+'_emission_rel_par_'+VV+'.csv'), sep =';')

    

#%%
# Remove large redundant list of data	

#del var_dict, params_list, tec_dict, BLUES_tec_df, link_t
if Load_info:
    del original_BLUES_df

    #%%
# =============================================================================
# capacity_factor_par_old = capacity_factor_par.copy()
# print('working on capacity factor adjustment')
# #capacity_factor_par = capacity_factor_par_old.copy()
# for t in list( set( capacity_factor_par['tec'].values.tolist() ) ):
#     if all(capacity_factor_par[capacity_factor_par['tec'] == t]['value'] == 1 ): #is_unique
#         capacity_factor_par = capacity_factor_par[capacity_factor_par['tec'] != t]
#         
# del capacity_factor_par_old
# =============================================================================
#%%#

	