# -*- coding: utf-8 -*-
"""
Created on Sat May  9 18:53:28 2020

@author: FD
"""
#%%
startTimeg = datetime.now()

#os.chdir(path) #setting the path of your working drive as the working drive directory
#sys.path.append('D:\\Fabio\\Diuana_Pdrive_IIASA\\Models\\NEST-BLUES\\Input_data_scripts') #set path into the system to avoid issues #REPLACE THIS ONE
sys.path.append(os.getcwd()) #set paths into the system to avoid issues
sys.path.append(path_script) #set paths into the system to avoid issues
sys.path.append(path_case) #set paths into the system to avoid issues
#from xfuncs import * #importing auxiliary functions
from BLUES_conversion_r169_functions import * #import from BLUES conversion functios
#from BLUES_ldb_building_technologies_FD_r13 import * #import from BLUES conversion functios    
from BLUES_building_technologies_FD_r163 import * #import from BLUES conversion functios   
    #%%
# =============================================================================
# GETTING BASIC INFO AND DEFINING TIME AND YEARS PARAMETERS
# =============================================================================
#reading line by line
l = open(os.path.join(path_case, project, main_region,"data",main_region+".adb"), "r") #open main region file
l1 = l.read() #reading adb
l2 = l1.split('\n') # split by line

#getting years
years = [t for t in l2[0:10] if t.startswith('timest')] 
years = years[0].split('timesteps: ')[-1] 
years = [int(y) for y in years.split(' ')]
year_all = years[:]

try:
    count_y =  years.index(first_model_year)-1 #define year index according to the first year model defined
except ValueError:
    raise Exception("first_model_year is not among years in this model - check it") 

year_all = year_all[count_y:-1] #defining years according to the first year model
vtgs = years[::] #set vintages
year_act = year_all[count_y+1:] #year_act



main_d_period = year_act[1] - year_act[0]
first_d_period = year_act[0] - vtgs[0] #last duration period
drate = [float(t.split(' ')[-1]) for t in l2[0:10] if t.startswith('drate')][0]/100 #discount rate
#first_hist_new_cap = first_model_year - year_hist_new_cap
start_year = {year_act[y]:years[y] +1   for y in range(len(year_act))}
#start_year[vtgs[-1]] = vtgs[-1] 
start_year[vtgs[0]] = vtgs[0] 

del l, l1, l2 #DELETING SOME VARIABLES

#time values
if has_ldr:
    if time1 not in ['months', 'month-hours', 'original']: #check option defined
        raise Exception("time should be defined as 'months', 'month-hours' or 'original'")
    if time1 == 'months': #if months it is from 1 to 12
        time_duration = [t/365 for t in [31,28,31,30,31,30,31,31,30,31,30,31]] #get time duration for months
        time = list(range(1,13))  #time for months from 1 to 12
    elif time1 == 'month-hours': #if it is month-hours it is 12 months * 24 hours == 288
        time = list(range(1,289)) 
        time_duration = [[tt/24]*24 for tt in [t/365 for t in [31,28,31,30,31,30,31,31,30,31,30,31]]] #get time duration for 288 teme steps
        time_duration = list(itertools.chain(*time_duration)) #flattening the list
    elif time1 == 'original': #if it is original it will take the values from adb/ldb database
        dt = os.path.join(path_case, project, main_region,"data", main_region+".adb")  #adb path
        if BLUES:
            dt = os.path.join(path_case, project, 'SE',"data", "SE.adb")  #adb path
        dt = open(dt) #reading abd
        
        time_block = "" #creating block
        found = False #starting
        for line in dt: #looping through lines
            if found: #if found is true add line into block
                time_block += line
                if line.strip() == 'energyforms:':  #break loop if line is equal this condition
                    break
            else: #if found is false
                if line.startswith('length'):   #set found as true if line attend this condition                 
                    found = True
                    time_block += line
    
        
        time_duration = [r_t_z("{:.11f}".format(float(t))) for t in re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", time_block)] #getting values
        time  = list(range(1, len(    time_duration)+1))        #time according to number of values
        #time_duration = [float_round(t, 4,floor) for t in time_duration]
        #time_duration = [round(t, 6) for t in time_duration]
        time_duration = [float(t) for t in time_duration]
else:        
    time = ['year']   
    time_duration = [1]
    
all_time_info = (year_all, count_y, time, first_d_period, main_d_period, start_year)#, first_hist_new_cap)
#%%
#to guarantee that ldb is empty in case task = 'ADB'
if task == 'ADB':
    ldb = ""
    ldb_aux = ldb
    ldb_from_adb_dict = False
elif task == 'LDB':
    if ldb == "":
        ldb = input("Set your ldb")
    ldb_aux = '_'+ldb
else:
    raise Exception("Task must be 'ADB' or 'LDB'")    
    
tpp = (task, path_case, project)    
#%%

Load_variables = [ 
                   str("original BLUES_tec_df") + ' at ' + os.path.join(path_case,project,'IX',file_version,ldb,project+'_'+task+ldb_aux+"_original_df.csv"),
                   str('BLUES_tec_df') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb, project+'_'+task+'_'+'df_'+file_version+'.csv'),
                   str('demand_df') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb,project+'_' + task+'_demand_df'+ldb_aux+'.csv'),
                   str('resource_df') + ' at '+os.path.join(path_case,project, 'IX',file_version,ldb, project+'_' + task+'_resource_'+'df_'+ldb+'.csv'),
                   str('land_rel_df') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_land_relation_df'+ldb_aux+'.csv'),
                   str('emission_df') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_emission_relation_df'+ldb_aux+'.csv'),
                   str('rel_df_no_land') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_relation_df'+ldb_aux+'.csv'),
                   str('link_t') + ' at '+ os.path.join(path_case,project,'IX',file_version,ldb,project+'_' + task+'_link_t_info'+ldb_aux+'.npy'),
                   str('lvl_com_t') + ' at '+ os.path.join(path_case,project,'IX',file_version,ldb,project+'_'+task+'_lvl_com_t_info'+ldb_aux+'.npy'),
                   str('ldr') + ' at '+os.path.join(path_case, project,'IX',file_version,ldb, project+"_"+task+"_ldr"+ldb_aux),
                   str('bal_eq.') + ' at '+os.path.join(path_case,project,'IX',file_version, ldb,project+'_'+task+ldb_aux+'bal_eq_df'+'.csv'),
                   str('ix data') + ' at '+os.path.join(path_case, project,'IX',file_version,ldb, project+'_' + task + '_land_ix_data' + ldb_aux),
                   str('emission_bound') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_emission_bound'+ldb_aux+'_df'+'.csv'),
                   str('emission_cost') + ' at '+ os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_emission_cost'+ldb_aux+'_df'+'.csv'),
                   str('relation_upper') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_relation_upper'+ldb_aux+'_df'+'.csv'),
                   str('relation_lower') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_relation_lower'+ldb_aux+'_df'+'.csv'),
                   str('relation_cost') + ' at '+os.path.join(path_case,project,'IX', file_version,ldb,project+'_' + task+'_relation_cost'+ldb_aux+'_df'+'.csv'),
                   str("old out tecs")  + ' at ' + os.path.join(path_case, project,'IX',file_version, ldb,project + '_' + task+"_tn_old_out"+ldb_aux),
                   str("tecs_f") + ' at ' + os.path.join(path_case, project,'IX',file_version, ldb,project + '_' + task+"_tecs_f"+ldb_aux),
                   str('warning list') + ' at '+os.path.join(path_case, project,'IX',file_version,ldb, project+'_'+task+"_warning_list"+ldb_aux+'_'+file_version),
                   str('tec_dict') + ' at '+os.path.join(path_case,project,'IX',file_version,ldb,project + '_'+task+ "_tec_dict_ix_format"+ldb_aux+'_'+file_version+".npy")
              ]

#%%
k = 0
original_BLUES_df = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb,project+'_'+task+ldb_aux+"_original_df.csv"), index_col=[0,1,2], header=[0,1],sep=';',dtype='object')

#land_data = pd.read_csv(os.path.join(path,project, project+'_'+task+'_'+'land_data_df'+'.csv')) land data

k+=1
#LOAD BLUES TEC DF DATA
BLUES_tec_df = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb, project+'_'+task+ldb_aux+'_df_'+file_version+'.csv'), index_col=[0,1,2], header=[0,1],sep=';', dtype='object') #exporting adb_df
BLUES_tec_df = BLUES_tec_df.fillna('nan')
#%%
#load them
#demand
k+=1
demand_df = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb,project+'_' + task+'_demand_df'+ldb_aux+'.csv'),sep=';')    #export to csv

k+=1
#get the resource info from .adb files --> OK
resource_df = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_resource_'+'df'+ldb_aux+'.csv'),sep=';') #exporting to csv    --> level Resources  

k+=1
if BLUES:
    #land relations

    land_rel_df = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_land_relation_df'+ldb_aux+'.csv'),sep=';')

k+=1
#emissions
emission_rel_df = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_emission_relation_df'+ldb_aux+'.csv'),sep=';')        
k+=1
#relations
rel_df_no_land = pd.read_csv(os.path.join(path_case,project, 'IX',file_version,ldb,project+'_' + task+'_relation_df'+ldb_aux+'.csv'),sep=';')

#%%
#LOAD LINK TEC
k+=1
link_t = np.load(os.path.join(path_case,project,'IX',file_version,ldb,project+'_' + task+'_link_t_info'+ldb_aux+'.npy'), allow_pickle=True).item()
k+=1
lvl_com_t = np.load(os.path.join(path_case,project,'IX',file_version,ldb,project+'_'+task+'_lvl_com_t_info'+ldb_aux+'.npy'), allow_pickle=True).item()
#%%

# =============================================================================
# #DO I NEED TO LOAD LDR??
# =============================================================================
k+=1
with open(os.path.join(path_case, project,'IX', file_version,ldb,project+"_"+task+"_ldr"+ldb_aux), 'rb') as fp: #loading list
    ldr = (pickle.load(fp)) #loading ldr
    ldr = ldr[0] #adjusting it
    
ldr_info = np.load(os.path.join(path_case, project, 'IX',file_version,ldb, project + '_'+task+"_ldr_info"+ldb_aux+'.npy'), allow_pickle=True).item()
#%%
#balance equality commodity level
k+=1
balance_eq_df = pd.read_csv(os.path.join(path_case,project,'IX',file_version, ldb,project+'_'+task+ldb_aux+'_bal_eq_df'+'.csv'), sep =';') 
#%%
k+=1
#LOAD LAND DATA IX INFO
if BLUES:

    with open(os.path.join(path_case, project, 'IX',file_version,ldb,project+'_' + task + '_land_ix_data' + ldb_aux), 'rb') as fp: #save list
        cat_land_tec_df, map_land_rel_df, sec_land_tec, land_tec, class_land_tec, land_generation = pickle.load(fp)
k+=1 
#%%
#LOAD EMISSION BOUND COST        
#bound_emission	node | type_emission | type_tec | type_year --> adjust emission upper
emission_bound = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_emission_bound'+ldb_aux+'_df'+'.csv'),sep=';') #exporting upper dataframe to csv
k+=1
emission_cost = pd.read_csv(os.path.join(path_case,project,'IX',file_version, ldb,project+'_' + task+'_emission_cost'+ldb_aux+'_df'+'.csv'),sep=';') #exporting upper dataframe to csv
k+=1
#LOAD RELATION UPPER, LOWER, COST 
relation_upper = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_relation_upper'+ldb_aux+'_df'+'.csv'),sep=';') #exporting upper dataframe to csv
k+=1
relation_lower = pd.read_csv(os.path.join(path_case,project,'IX',file_version,ldb, project+'_' + task+'_relation_lower'+ldb_aux+'_df'+'.csv'),sep=';') #exporting lower dataframe to csv
k+=1
relation_cost = pd.read_csv(os.path.join(path_case,project, 'IX',file_version,ldb,project+'_' + task+'_relation_cost'+ldb_aux+'_df'+'.csv'),sep=';') #exporting upper dataframe to csv

emission_bound['type_tec'] = 'all'
emission_cost['type_tec'] = 'all'
#%%
#LOADING OLD TEC DATA
k+=1 
with open(os.path.join(path_case, project,'IX', file_version,ldb,project + '_' + task+"_old_tec_list"+ldb_aux), 'rb') as fp: #save list
    old_tec_list = pickle.load(fp) #save ldr list    
    #old_tec_list = [tn_old, tn_old_out]  
#%%
#tecs f
k+=1
with open(os.path.join(path_case, project,'IX',file_version, ldb,project + '_' + task+"_tecs_f"+ldb_aux), 'rb') as fp: #save list
    tecs_f = pickle.load(fp) #save ldr list    
    #old_tec_list = [tn_old, tn_old_out]  

#%%
# =============================================================================
# 
# if new_dict_from_load:
#     if ldb_from_adb_dict:
#         BLUES_tec_df = BLUES_tec_df.loc[~(BLUES_tec_df=='nan').all(axis=1)] #removing rows in which all valuea are nan
#         tecs_f = sorted(list(BLUES_tec_df.index.get_level_values(level=0).tolist() )) #dataframe tecs
#         tecs_f = sorted(list(set(tecs_f) ) )
#         tecs_f1 = [tt for tt in tecs_f if 'Exist' not in tt] #removing existing
#         
#         BLUES_tec_df.to_csv(os.path.join(path_case,project,'IX',file_version, ldb,project+'_only_'+task+ldb_aux+'_df'+'.csv'), sep=';') #save ldb dr
#         
#         adb_tec_dict = np.load(os.path.join(path_case,project,'IX',file_version,project + "_ADB_tec_dict_ix_format.npy"), allow_pickle=True).item() #loading tec idct from adb
#         
#         new_ldb_tec = [t for t in tecs_ldb if t not in tecs_adb] #getting new ldb tecs
#     
#         new_warning_list, new_tec_dict, not_def_rel = build_adb_dfs(tpp, ldb, BLUES_tec_df, new_ldb_tec, all_time_info, old_tec_list, ldr, ldr_info, land_rel_df, rel_df_no_land,BLUES)
#         #if it work
#     
#     # =============================================================================
#         with open(os.path.join(path_case, project,'IX',file_version,ldb, project+'_'+task+"_ONLY_warning_list"+ldb_aux), 'wb') as fp: #save list
#             pickle.dump(new_warning_list, fp) #save ldr list              
#         
#         np.save(os.path.join(path_case, project,'IX', file_version,ldb, project + '_'+task+"_ONLY_tec_dict_ix_format"+ldb_aux+'.npy'),new_tec_dict)
#     # =============================================================================
#     
#         removed_adb_tecs = [t for t in tecs_adb if t not in tecs_ldb] #removing tecs from adb dict bc they are not considered in ldb tecs
#         for r in removed_adb_tecs: #LOOPING THROUGH TECS THAT MUST BE REMOVED
#             if adb_tec_dict.has_key(r): #CHECK IF DICT CONTAINS R TEC
#                 del adb_tec_dict[r] #REMOVE
#                 
#         ntec_dict = copy.deepcopy(adb_tec_dict)#.copy()     #APPENDING NEW TEC INTO ADB DICT 
#         ntec_dict.update(new_tec_dict)
#     
#     # =============================================================================
#     #     #REVIEW BUILD LDB  --> new old tech approach
#     #     warning_list, tec_dict, not_def_rel = build_ldb_dfs(tpp, ldb, BLUES_tec_df, ntec_dict, tecs_f1, all_time_info, tn_old,  ldr, land_rel_df, rel_df_no_land, ldb_ldr, BLUES=True)
#     # 
#     # =============================================================================
#       
#     else:    
#         #BLUES_tec_df.to_csv(os.path.join(path,project, ldb,project+'_'+task+ldb_aux+'_df'+'.csv'), sep =';') #exporting adb_df
#     
#         warning_list, tec_dict, not_def_rel = build_adb_dfs(tpp, ldb, BLUES_tec_df, tecs_f, all_time_info, old_tec_list, ldr, ldr_info, land_rel_df, rel_df_no_land,BLUES)
#         
#     # =============================================================================
#     #     IF EVERYTHING WORKS FINE TAKE THE FOLLOWING LINES off IF ELSE INDENT
#     # =============================================================================
#     # =============================================================================
#     #     with open(os.path.join(path, project,ldb, project+'_'+task+"_warning_list"+ldb_aux), 'wb') as fp: #save list
#     #         pickle.dump(warning_list, fp) #save ldr list              
#     #     
#     #     np.save(os.path.join(path, project, ldb, project + '_'+task+"_tec_dict_ix_format"+ldb_aux+'.npy'),tec_dict)
#     # =============================================================================
# 
#     #%%
#     with open(os.path.join(path_case, project,'IX',file_version,ldb, project+'_'+task+"_warning_list"+ldb_aux + version_), 'wb') as fp: #save list
#         pickle.dump(warning_list, fp) #save ldr list              
#     
#     np.save(os.path.join(path_case, project,'IX', file_version,ldb, project + '_'+task+"_tec_dict_ix_format"+ldb_aux + version_+'.npy'),tec_dict)
#     
#     
#     tecs_f2 = sorted( list(tec_dict.keys()) )   
# 
#     tec_inp_out = {}
#     for tt in tecs_f2:
#         tec_inp_out[tt] = {}
#         inpl = list(set([(l,c) for l,c in zip(tec_dict[tt]['input']['level'], tec_dict[tt]['input']['commodity'])]))
#         outl = list(set([(l,c) for l,c in zip(tec_dict[tt]['output']['level'], tec_dict[tt]['output']['commodity'])]))
#         tec_inp_out[tt]['input'] = inpl
#         tec_inp_out[tt]['output'] = outl
#         
#     np.save(os.path.join(path_case, project,'IX',file_version, ldb, project + '_'+task+"_tec_inp-out"+ldb_aux+version_+'.npy'),tec_inp_out)
# 
# =============================================================================
    #%%

#LOAD warning list
k+=1
with open(os.path.join(path_case, project,'IX',file_version, ldb,project+'_'+task+"_warning_list"+ldb_aux + '_'+file_version), 'rb') as fp: #loading warning list
    warning_list = pickle.load(fp) #load warning  list          
k+=1
#LOAD DICTIONARY
tec_dict = np.load(os.path.join(path_case,project,'IX',file_version,ldb,project + '_'+task+ "_tec_dict_ix_format"+ldb_aux + '_'+file_version+".npy"), allow_pickle=True).item() #loading tec idct from adb
tecs_f2 = sorted( list( tec_dict.keys() ) )
tec_inp_out = np.load(os.path.join(path_case, project,'IX',file_version, ldb, project + '_'+task+"_tec_inp-out"+ldb_aux+'_'+file_version+'.npy'), allow_pickle=True).item() #loading tec idct from adb

#%%
dt_finalg =  datetime.now() - startTimeg
print(task+ldb_aux+' Done.')
print('Runtime of load data script:    '+ str(dt_finalg))
