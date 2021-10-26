# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 15:08:11 2018

@author: amendolf
"""
#mp.close_db()
#%%

#NO BALANCE EQUALITY
#NO BDA LO
#NO BDC LO
#NO BDI LO

#NO - REL LO
#del relation_lower_par


#%%
from datetime import datetime
istime = datetime.now()


# In[3]:

import ixmp
import message_ix
mp = ixmp.Platform(name='local', jvmargs=['-Xmx8G']) #['default', 'local'] #4g 6G 8G
#mp.jdbc.start_jvm(jvmargs=['-Xmx4G'])
#mp.backend.jdbc.JDBCBackend(jvmargs=['-Xmx4G'])
# =============================================================================
# scenario = message_ix.Scenario(mp, model='Westeros Electrified', 
#                                scenario='baseline', version='new')
# =============================================================================                                     


scenario = message_ix.Scenario(mp, model = model_name, 
                                   scenario=sc_name, version='new')
#scen = message_ix$Scenario(mp,model, scen)

#%%
# list of sets, maps and parameter available
# set_check_list = objects()[grep('.set',objects())]    # single index sets, vectors in R
# map_check_list = objects()[grep('.map',objects())]    # which have more that one index, dataframes in R
# par_check_list = objects()[grep('.par',objects())]    # parameters, dataframes in R

# Create lists containing parameters and force to character string
set_list_dic = dict(   node = node_set.copy(),
        			   year = year_set.copy(),
					   lvl_spatial = lvl_spatial_set.copy(),
					   #lvl_temporal = lvl_temporal_set.copy(),
					   time = time_set.copy(),
					   #type_tec = type_tec_set.copy(),
					   type_year = type_year_set.copy(),
					   type_emission = type_emission_set.copy(),

					   grade = grade_set.copy(),
					   mode = mode_set.copy(),
                   commodity = commodity_set.copy(),
					   level = level_set.copy(),
					   technology = tec_set.copy(),
					   emission = emission_set.copy(),
					   rating = rating_set.copy(),
                       
					)


        
map_list_dic = dict(   
                   cat_emission = cat_emission_map.copy(),
					   #cat_year = cat_year_set.copy(),
					   #cat_tec = cat_tec_set.copy(),
					   map_spatial_hierarchy = map_spatial_hierarchy_map.copy(),
					   #map_temporal_hierarchy = map_temporal_hierarchy_map.copy(),
					   map_time = map_time_set.copy(),
					   map_node = map_node_set.copy(),
                   #bal_eq_map = bal_eq_df_map.copy()

					)


if has_ldr:
    set_list_dic['lvl_temporal'] = lvl_temporal_set    
    map_list_dic['map_temporal_hierarchy'] = map_temporal_hierarchy_map.copy()

for c in set_list_dic.keys():
    if all(isinstance(item, str) for item in set_list_dic[c]) != True:
        set_list_dic[c] = list(map(str, set_list_dic[c]))

for m in map_list_dic.keys():
    map_list_dic[m] = map_list_dic[m].applymap(str)

#%%        
par_list_dic = dict(    
                        input = input_par.copy(),
                        capacity_factor = capacity_factor_par.copy(),
					    output = output_par.copy(),
					    inv_cost = inv_cost_par.copy(),

					    var_cost = var_cost_par.copy(),
					    emission_factor = emission_factor_par.copy(),
					    technical_lifetime = technical_lifetime_par.copy(), 
					    construction_time = construction_time_par.copy(),  
					    demand = demand_par.copy(),
#					    duration_period = duration_period_par, 
					    duration_time = duration_time_par.copy(),

					    fix_cost = fix_cost_par.copy(),
					    historical_new_capacity = historical_new_capacity_par.copy(), 
					    interestrate = interestrate_par.copy(), 
					    
                        min_utilization_factor = min_utilization_factor_par.copy(),
                        bound_activity_lo = bound_activity_lo_par.copy(),
                        
                        bound_activity_up = bound_activity_up_par.copy(),
                        bound_total_capacity_up = bound_total_capacity_up_par.copy(),
                        
                        bound_total_capacity_lo = bound_total_capacity_lo_par.copy(),
                        bound_new_capacity_up = bound_new_capacity_up_par.copy(),
                        
                        bound_new_capacity_lo = bound_new_capacity_lo_par.copy()                       
					)

#RESOURCE

if 'resource_volume_par' in locals():#if EMISS== True or SOL_WIND==True:
    par_list_dic['resource_volume'] = resource_volume_par.copy()
if 'resource_cost_par' in locals():#if EMISS== True or SOL_WIND==True:
    par_list_dic['resource_cost'] = resource_cost_par.copy()

#EMISSION
if 'emission_bound_par' in locals():#if EMISS== True or SOL_WIND==True:
    par_list_dic['bound_emission'] = emission_bound_par.copy()
if 'emission_cost_par' in locals():#if EMISS== True or SOL_WIND==True:
    par_list_dic['tax_emission'] = emission_cost_par.copy()


model_size = 0		  		  
for sp in par_list_dic.keys():
    print(sp)
    par_list_dic[sp] = par_list_dic[sp].applymap(str)   
    model_size+= len(par_list_dic[sp]) #last print --> 8950910


#par_list_dic['min_utilization_factor'].drop_duplicates(keep='first', inplace=True)
# =============================================================================
# sp = 'output'    
# par_list_dic[sp][par_list_dic[sp]['commodity']=='water_with']
# =============================================================================
#%%
count_nan = 0
for sp in par_list_dic.keys():
    if sp not in ['demand', 'duration_time', 'interestrate','resource_volume', 'resource_cost', 'tax_emission', 'bound_emission'] :
        print(sp)
        #par_list_dic[sp][(par_list_dic[sp].isin(['nan'])==True).any(1)]
        if len(list(set( par_list_dic[sp][(par_list_dic[sp].isin(['nan'])==True).any(1)]['tec'].values.tolist() )) ) > 0:
            count_nan += 1
            print(sp)
            print('\n')
            print(list(set( par_list_dic[sp][(par_list_dic[sp].isin(['nan'])==True).any(1)]['tec'].values.tolist() )) )
            print('____________________________________________________________________')
            print('\n')


# =============================================================================
#a = par_list_dic[sp][(par_list_dic[sp].isin(['nan'])==True).any(1)]
#ab = a[a['tec'] == 'Hh_Light_Inc_Base']
# =============================================================================

#check those parameters
if len(par_list_dic['demand'][(par_list_dic['demand'].isin(['nan'])==True).any(1)]) > 0:
    print(par_list_dic['demand'][(par_list_dic['demand'].isin(['nan'])==True).any(1)])
    count_nan += 1
if len(par_list_dic['resource_volume'][(par_list_dic['resource_volume'].isin(['nan'])==True).any(1)]) > 0:
    print(par_list_dic['resource_volume'][(par_list_dic['resource_volume'].isin(['nan'])==True).any(1)])
    count_nan += 1
if len(par_list_dic['duration_time'][(par_list_dic['duration_time'].isin(['nan'])==True).any(1)]) > 0:
    print(par_list_dic['duration_time'][(par_list_dic['duration_time'].isin(['nan'])==True).any(1)])
    count_nan += 1
if len(par_list_dic['interestrate'][(par_list_dic['interestrate'].isin(['nan'])==True).any(1)]) > 0:
    print(par_list_dic['interestrate'][(par_list_dic['interestrate'].isin(['nan'])==True).any(1)])
    count_nan += 1
if len(par_list_dic['resource_cost'][(par_list_dic['resource_cost'].isin(['nan'])==True).any(1)]) > 0:
    print(par_list_dic['resource_cost'][(par_list_dic['resource_cost'].isin(['nan'])==True).any(1)])
    count_nan += 1
if len(par_list_dic['tax_emission'][(par_list_dic['tax_emission'].isin(['nan'])==True).any(1)]) > 0:
    print(par_list_dic['tax_emission'][(par_list_dic['tax_emission'].isin(['nan'])==True).any(1)])
    count_nan += 1
if len(par_list_dic['bound_emission'][(par_list_dic['bound_emission'].isin(['nan'])==True).any(1)]) > 0:
    print(par_list_dic['bound_emission'][(par_list_dic['bound_emission'].isin(['nan'])==True).any(1)])
    count_nan += 1
if count_nan >0:
    print("it has nan values")        
    sys.exit()
#%%
# =============================================================================
# LOOKING FOR DUPLICATED ROWS IN PARAMETERS BEFORE ADD THEM TO IXMP
# =============================================================================   

count_duplicate = 0
for sp in par_list_dic.keys():
    if sp not in ['demand', 'duration_time', 'interestrate','resource_volume', 'resource_cost', 'tax_emission', 'bound_emission'] :

        duplicateDFRow = par_list_dic[sp][par_list_dic[sp].duplicated()]
        #par_list_dic[sp][(par_list_dic[sp].isin(['nan'])==True).any(1)]
        if len(duplicateDFRow) > 0:
            print(sp + ' has duplicate rows')
            count_duplicate += 1
            print('____________________________________________________________________')
            print('\n')
            



# =============================================================================
#a = par_list_dic[sp][(par_list_dic[sp].isin(['nan'])==True).any(1)]
        
#ab = a[a['tec'] == 'Hh_Light_Inc_Base']
# =============================================================================

#check those parameters
duplicateDFRow = par_list_dic['demand'][par_list_dic['demand'].duplicated()]
if len(duplicateDFRow) >0:
    print('demand'+ ' has duplicate rows')
    count_duplicate += 1
    
duplicateDFRow = par_list_dic['resource_volume'][par_list_dic['resource_volume'].duplicated()]
if len(duplicateDFRow) >0:
    print('resource_volume' + ' has duplicate rows')
    count_duplicate += 1
    
duplicateDFRow = par_list_dic['duration_time'][par_list_dic['duration_time'].duplicated()]
if len(duplicateDFRow) >0:
    print('duration_time' + ' has duplicate rows')
    count_duplicate += 1
    
duplicateDFRow = par_list_dic['interestrate'][par_list_dic['interestrate'].duplicated()]
if len(duplicateDFRow) >0:
    print('interestrate' + ' has duplicate rows')
    count_duplicate += 1
    
duplicateDFRow = par_list_dic['resource_cost'][par_list_dic['resource_cost'].duplicated()]    
if len(duplicateDFRow) >0:
    print('resource_cost' + ' has duplicate rows')
    count_duplicate += 1
    
duplicateDFRow = par_list_dic['tax_emission'][par_list_dic['tax_emission'].duplicated()]    
if len(duplicateDFRow) >0:
    print('tax_emission' + ' has duplicate rows')
    count_duplicate += 1
    
duplicateDFRow = par_list_dic['bound_emission'][par_list_dic['bound_emission'].duplicated()]    
if len(duplicateDFRow) >0:
    print('bound_emission' + ' has duplicate rows')
    count_duplicate += 1
 
if count_duplicate>0:
    print("\n\nit has duplicated values\n\n")
    sys.exit()
#else:
#    del par_list_dic
#%%
    
del input_par,capacity_factor_par, #output_par
del    inv_cost_par,var_cost_par,emission_factor_par
del    technical_lifetime_par,construction_time_par,  demand_par
#					    duration_period = duration_period_par, duration_time_par,
del     fix_cost_par, historical_new_capacity_par,  interestrate_par
del     min_utilization_factor_par, bound_activity_lo_par
del     bound_activity_up_par, bound_total_capacity_up_par
del     bound_total_capacity_lo_par, bound_new_capacity_up_par
del     bound_new_capacity_lo_par        
        
#%%
scenario.add_horizon({'year': year_all,  
                      'firstmodelyear': year_act[0]})    
    

    ## Add SETs
for s in set_list_dic.keys():
    print(s)
    scenario.add_set(s, set_list_dic[s])
    

# =============================================================================
#for t in  set_list_dic[s]:
#    print(t)
#    scenario.add_set(s, t)
# =============================================================================
    
scenario.add_set("level_resource", level_resource_par)   
scenario.add_par('duration_period', duration_period_par.applymap(str)   ) 
scenario.add_par('duration_time', duration_time_par.applymap(str)    )

scenario.add_set("balance_equality",  bal_eq_df_map)    


## Add MAPPING SETs
for m in map_list_dic.keys():
    print(m)
    tmp_name = list(scenario.idx_names(m))
    tmp_map = map_list_dic[m]
    #print(tmp_name)
    #print(list(tmp_map.columns))
    tmp_map.columns = tmp_name
    scenario.add_set(m, tmp_map)

del tmp_map, map_list_dic, set_list_dic
#%%
##### INITIALIZE NEW SETS, VARIABLES(PARAMETERS) AND EQUATIONS for RELATIONS


scenario.init_set( "relation2" )
scenario.init_set( "map_relation2", ['relation2','node','year'], ['relation2','node_rel', 'year_act'])

scenario.init_set( "is_relation_upper2", ['relation2','node','year'], ['relation2','node_rel', 'year_rel'] )
scenario.init_set( "is_relation_lower2", ['relation2','node','year'], ['relation2','node_rel', 'year_rel'] )

scenario.init_par('relation_upper2', ['relation2','node','year'],
                                      ['relation2', 'node_rel','year_rel'])

scenario.init_par('relation_lower2', ['relation2','node','year'],
                                     ['relation2', 'node_rel','year_rel'])

scenario.init_par('relation_cost2', ['relation2','node','year'], ['relation2','node_rel','year_rel'])

scenario.init_par('relation_new_capacity2', ['relation2','node','year', 'technology'])
scenario.init_par('relation_total_capacity2', ['relation2','node','year', 'technology'])

scenario.init_par('relation_activity2', ['relation2','node','year', 'node', 'technology', 'year', 'mode'], 
                  ['relation2','node_rel','year_rel','node_loc', 'tec', 'year_act','mode'])

scenario.init_equ("RELATION_EQUIVALENCE2", [ 'relation2', 'node','year' ] )
scenario.init_equ("RELATION_CONSTRAINT_UP2", [ 'relation2', 'node','year' ] )
scenario.init_equ("RELATION_CONSTRAINT_LO2", [ 'relation2', 'node','year' ] )

#%%
# ADD VALUES for RELATIONS
if 'relation_set' in locals():#
    scenario.add_set('relation2', relation_set )
    par_list_dic['relation_activity2'] = relation_activity_par
    map_relation_par.rename(columns={'relation':'relation2'}, inplace=True)
    #map_relation_ix = map_relation_par.copy()
    #map_relation_ix['value'] = True
    scenario.add_set("map_relation2",  map_relation_par) #do it   
        
if 'relation_total_capacity_par' in locals():#
    par_list_dic['relation_total_capacity2'] = relation_total_capacity_par
    
if 'relation_new_capacity_par' in locals():#
    par_list_dic['relation_new_capacity2'] = relation_new_capacity_par    
    
#RELATION VALLUES
if 'relation_cost_par' in locals():#
    relation_cost_par.rename(columns={'relation':'relation2'}, inplace=True)
    par_list_dic['relation_cost2'] = relation_cost_par
    
if 'relation_upper_par' in locals():#
    relation_upper_par.rename(columns={'relation':'relation2'}, inplace=True)
    par_list_dic['relation_upper2'] = relation_upper_par

    is_relation_upper_par.rename(columns={'relation':'relation2'}, inplace=True)
    #is_relation_upper_ix['value'] = True
    scenario.add_set("is_relation_upper2", is_relation_upper_par)  
        
if 'relation_lower_par' in locals():#
    relation_lower_par.rename(columns={'relation':'relation2'}, inplace=True)
    par_list_dic['relation_lower2'] = relation_lower_par
    
    is_relation_lower_par.rename(columns={'relation':'relation2'}, inplace=True)
    #is_relation_lower_ix['value'] = True
    scenario.add_set("is_relation_lower2", is_relation_lower_par) #do it         

#%%
    
# =============================================================================
# CONVERT INJC --> relationsc --> into RESOURCE
    #CHECK IF TECS WITH NO MAP_TEC_LIFETIME ARE BEING CONSIDERED ON RELATION EQUATIONS -->APPARENTLY IT IS OK
# =============================================================================
    
#del tec_dict
startTime2 = datetime.now()
## Add PARAMETERs
for p in list(par_list_dic.keys()):#[4:]:
    print('\n...'+p+' ---> Loading...')
    print(list(scenario.idx_names(p)))
    print(par_list_dic[p].columns.tolist())
    tmp_name = list(scenario.idx_names(p)) + ['value']
    tmp_par = par_list_dic[p]
    #print('tmp_name:\n')
    #print(tmp_name)
    #print('cols:\n')
    #print(list(tmp_par.columns))
    tmp_par.columns = tmp_name
    #tmp_par['unit'] = '-'
    if len(tmp_par) >0:
        scenario.add_par(p, tmp_par)
        par_list_dic[p] = None

    print('\n...'+p+' ---> Done.')

Runtime_par = datetime.now() - startTime2 
print(Runtime_par)       
del tmp_par, par_list_dic

#%%

scenario.init_par('emission_duration_period', ['year'] )
scenario.add_par('emission_duration_period', emission_duration_period_par)

scenario.init_par('duration_period2', ['year'] )
scenario.add_par('duration_period2', duration_period_par.applymap(str) )
#%%
#adjust it
print("\nadding main output value")
scenario.init_par('main_output_val', ['node', "technology",    'year',     'year',   'mode', 'time'],
                                      ['node','tec',  'vintage', 'year_act', 'mode',  'time'])    

scenario.add_par('main_output_val', main_output_par)

print("\nadding main input value")
scenario.init_par('main_input_val', ['node', "technology",    'year',     'year',   'mode', 'time'],
                                      ['node','tec',  'vintage', 'year_act', 'mode',  'time'])    

scenario.add_par('main_input_val', main_input_par)
#%%
# =============================================================================
print("\nadding map tec vtg")
scenario.init_set( "map_tec_vtg", ['node','technology','year'], ['node','tec','vintage'] ) 
# 
map_tec_vtg_df = output_par[['node', 'tec', 'vintage']]
map_tec_vtg_df.drop_duplicates(keep='first', inplace=True)
# 
scenario.add_set("map_tec_vtg",  map_tec_vtg_df)
# =============================================================================
print("\nadding map tec vtg act")    
scenario.init_set( "map_tec_vtg_act", ['node','technology','year', 'year'], ['node','tec','vintage', 'year_act'] ) 

map_tec_vtg_act_df = output_par[['node', 'tec', 'vintage', 'year_act']]
map_tec_vtg_act_df.drop_duplicates(keep='first', inplace=True)

scenario.add_set("map_tec_vtg_act",  map_tec_vtg_act_df) 
#%%
##### INITIALIZE NEW SETS, VARIABLES(PARAMETERS) AND EQUATIONS for BLUES BALANCE

#New set, par and eq.
scenario.init_set('BLUES_land_tec', idx_sets='technology')
scenario.init_set('BLUES_sec_land_tec', idx_sets='technology')
scenario.init_set( "BLUES_type_land_tec" )
scenario.init_set( "BLUES_cat_land_tec", [ 'BLUES_type_land_tec', 'technology' ], [ 'BLUES_type_land_tec', 'tec' ])
scenario.init_set("map_BLUES_land_rel", [ 'node', "technology", 'BLUES_type_land_tec' ], [ 'node', "tec", 'BLUES_type_land_tec' ] )

scenario.init_equ("BLUES_LAND_BALANCE", [ 'node','technology', 'year', 'BLUES_type_land_tec' ] )
#scenario.init_equ("BLUES_LAND_HIST_BAL", [ 'node','technology', 'year' ] )
#%%
# ADD VALUES for LAND SETS
print("adding new BLUIES land set and parameters")
for lt in land_tec:
    scenario.add_set('BLUES_land_tec', [lt])

for slt in sec_land_tec:
    scenario.add_set('BLUES_sec_land_tec', [slt])

scenario.add_set("BLUES_type_land_tec", class_land_tec ) #['f2', 'p2', 'c2'])

scenario.add_set("BLUES_cat_land_tec", cat_land_tec_df)

scenario.add_set("map_BLUES_land_rel", map_land_rel_df)
#scenario.set("map_BLUES_land_rel")
#%%
##### INITIALIZE NEW VARIABLES(PARAMETERS) AND EQUATIONS for MIN UTILIZATION FACTOR
print("New MUF constraints")
scenario.init_par('min_utilization_time_factor', ['node','technology','year', 'year', 'time'],
                                                     ['node','tec','vintage', 'year_act', 'time'])


scenario.init_equ("MIN_UTILIZATION_CONSTRAINT_TIME", ['node','technology', 'year', 'year', 'time' ], 
                                                          ['node','tec','vintage', 'year_act', 'time'] )

if has_ldr:#if muf_time:
    scenario.add_par('min_utilization_time_factor', min_utilization_time_factor_par)


#%%
##### INITIALIZE NEW VARIABLES(PARAMETERS) AND EQUATIONS for MIN UTILIZATION FACTOR

scenario.init_par('var_cost2', ['node','technology','year', 'year', 'mode'],
                                                     ['node','tec','vintage', 'year_act', 'mode'])

if no_vom_time == True:#if muf_time:
    scenario.add_par('var_cost2', var_cost_par2)
    
# =============================================================================
# import inspect
# print(os.path.abspath(inspect.getfile(scenario.add_par)))
# print(os.path.abspath(inspect.getfile(scenario.solve)))
# =============================================================================
#%%
#### COMMIT and output to GDX
startTime3 = datetime.now()

from message_ix import log
print( 'Generating GDX' )

log.info('version number prior to commit: {}'.format(scenario.version))
print('log info --> Done')

scenario.commit(comment='BLUES_LB_TEST')
print('commit --> Done')

#log.info('version number prior committing to the database: {}'.format(scenario.version))

Runtime_par2 = datetime.now() - startTime3 
print(Runtime_par2)   

# An `ixmp` database can contain many scenarios, and possibly multiple versions of the same model and scenario name.
# These are distinguished by unique version numbers.
# 
# To make it easier to retrieve the "correct" version (e.g., the latest one), you can set a specific scenario as the default version to use if the "Westeros Electrified" model is loaded from the `ixmp` database.

# In[33]:


scenario.set_as_default()


#%%
#### SOLVE LOCALLY
startTime5 = datetime.now()
print( 'Solving...' )

#current_working_drive = os.getcwd()
#os.chdir(os.path.join( basin_ix_path,'model' ) ) # set based on local machine

#runm = 'MESSAGE'
#cmd =  "gams "+runm+"_run.gms --in=data\\MSGdata_"+scname+".gdx --out=output\\MSGoutput_"+scname+".gdx"
#res = system(cmd)
pathix = "\\".join(message_ix.__file__.split('\\')[:-1] + [ 'model', 'data'] ) 
datafname = 'MsgData_'+model_name+'_'+sc_name+'.gdx'
try:
    scenario.solve()
except:
    if datafname in os.listdir(pathix):
        print('\n'+datafname + ' is in '+ pathix )
        if datafname not in os.listdir(os.path.join(path_model, 'data')):
            shutil.copy(os.path.join(pathix, datafname), os.path.join(path_model, 'data', datafname) )    
        else:
            oldd = [ll.split('old_')[-1].split('.')[0] for ll in os.listdir(os.path.join(path_model, 'data')) if ll.startswith(datafname.split('.gdx')[0]+'_old')]
            if len(oldd) == 0:
                shutil.copy(os.path.join(path_model, 'data', datafname), os.path.join(path_model, 'data', datafname.split('.gdx')[0]+'_old_01.gdx') )
            else:
                olddv = max(oldd)
                if olddv[0] == '0':
                    olddv = str(int(olddv[-1])+1)
                else:
                    olddv = str(int(olddv)+1)
                    
                if len(olddv) == 1:
                    olddv = '0'+olddv
                    
                shutil.copy(os.path.join(path_model, 'data', datafname), os.path.join(path_model, 'data', datafname.split('.gdx')[0]+'_old_'+olddv+'.gdx')  )              
            
            shutil.copy(os.path.join(pathix, datafname), os.path.join(path_model, 'data', datafname) )
        print('\n'+datafname + ' was copied for '+ os.path.join(path_model, 'data') )
        
    else:
        print('solve did not generate gdx input data file')
        
Runtime_par5 = datetime.now() - startTime5
print(Runtime_par5)  

Runtime_par6 = datetime.now() - istime
print("\nTotal elapsed time:")
print(Runtime_par6)  
#os.chdir(current_working_drive)

	
   
#%%
mp.close_db()    

