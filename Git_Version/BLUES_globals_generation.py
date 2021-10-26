# -*- coding: utf-8 -*-
"""
Created on Sat May  9 18:53:28 2020

@author: FD
"""
#%%

# =============================================================================
# SETTING BASICS
# =============================================================================
startTimeg = datetime.now()
#os.chdir(path) #setting the path of your working drive as the working drive directory
#sys.path.append('D:\\Fabio\\Diuana_Pdrive_IIASA\\Models\\NEST-BLUES\\Input_data_scripts') #set path into the system to avoid issues #REPLACE THIS ONE
sys.path.append(os.getcwd()) #set paths into the system to avoid issues
sys.path.append(path_script) #set paths into the system to avoid issues
sys.path.append(path_case) #set paths into the system to avoid issues
#from xfuncs import * #importing auxiliary functions
from BLUES_building_technologies_FD_r164 import * #import from BLUES conversion functios    
#from BLUES_ldb_building_technologies_FD_r13 import * #import from BLUES conversion functios    
from BLUES_conversion_r170_functions import * #import from BLUES conversion functios

try:
    # Create target Directory
    os.mkdir(os.path.join(path_case, project, 'IX'))
    print( "IX Directory " ,  " Created ") 
except FileExistsError:
    print("IX Directory " ,  " already exists")

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
if period_interp == 'recursive': #
    start_year = {year_act[y]:years[y] +1   for y in range(len(year_act))}
    #start_year[vtgs[-1]] = vtgs[-1] 
    start_year[vtgs[0]] = vtgs[0] - (first_model_year-start_year[first_model_year])
else:
    start_year = {year_act[y]:year_act[y]  for y in range(len(year_act))}
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
# =============================================================================
# GET DATA FROM .ADB or .LDB FILE
# =============================================================================
#to guarantee that ldb is empty in case task = 'ADB'
if task == 'ADB':
    ldb = ""
    ldb_aux = ldb
    ldb_from_adb_dict = False
elif task == 'LDB':
    if ldb == "":
        ldb = input("Set your ldb")
    ldb_aux = '_'+ldb
    try:
        # Create target Directory
        os.mkdir(os.path.join(path_case, project, 'IX', ldb))
        print(ldb + " Directory " ,  " Created ") 
    except FileExistsError:
        print(ldb + " Directory " ,  " already exists")
else:
    raise Exception("Task must be 'ADB' or 'LDB'")    

tpp = (task, path_case, project)    

if new_df == False:
    if file_version == 'last':
        files = os.listdir(os.path.join(path_case, project, 'IX'))
        files = [ll for ll in files if not os.path.isfile(os.path.join(path_case, project, 'IX', ll))]    
        file_version = str(max([int(ll) for ll in files]))     
    shutil.copy(os.path.join(path_case,project,'IX',file_version, ldb,project+'_'+task+"_original_df.csv"), os.path.join(path_case,project,'IX', ldb,project+'_'+task+"_original_df.csv") )
    shutil.copy(os.path.join(path_case,project,'IX',file_version,ldb,project+'_'+task+'_energy_forms_dict'+ldb_aux+'.npy'), os.path.join(path_case,project,'IX',ldb,project+'_'+task+'_energy_forms_dict'+ldb_aux+'.npy' ) ) 

    shutil.copy(os.path.join(path_case, project,'IX',file_version, project+'_'+task+'_land_relation_list'), os.path.join(path_case, project,'IX',project+'_'+task+'_land_relation_list'))
    shutil.copy(os.path.join(path_case, project,'IX',file_version, project+'_'+task+'_emission_list') , os.path.join(path_case, project,'IX',project+'_'+task+'_emission_list'))
    shutil.copy(os.path.join(path_case, project,'IX',file_version, project+'_'+task+'_relation_list'), os.path.join(path_case, project,'IX', project+'_'+task+'_relation_list'))


df = task_function(tpp, adb_files,  vtgs, year_act, BLUES, emissions_constraints, ldb, new_df)

#%%
# =============================================================================
# SETTING DATAFRAME THAT WILL BE DEFINED AS WORKING DF
# =============================================================================
if task =='LDB':
    ldb_df = df[-1].copy()     #ldb df
    adb_df = df[0].copy() #adb df   
    adb_ldb_df = df[1].copy() #adb + ldb df  
    
    tecs_adb0 = adb_df.index.get_level_values(level=0).tolist() #adb tecs
    tecs_adb = sorted(list(set(tecs_adb0)))
    tecs_adb_ldb0 = adb_ldb_df.index.get_level_values(level=0).tolist() #adb tecs
    tecs_adb_ldb = sorted(list(set(tecs_adb_ldb0)))
    #tec_dict_adb = np.load(os.path.join(path,project, project+"_tec_dict_ix_format_ADB.npy"), allow_pickle=True).item() #loading tec idct from adb
    #tecs_adb_dict = list(tec_dict_adb.keys()) #list o f tecs in the ab dictionary
    
    tecs_ldb = ldb_df.index.get_level_values(level=0).tolist() #ldb tecs
    tecs_ldb = sorted(list(set(tecs_ldb))) #removing duplicate tec values

    if ldb_from_adb_dict:
        BLUES_tec_df = ldb_df.copy() #setting ldb df as the working df
        tecs = BLUES_tec_df.index.get_level_values(level=0).tolist() #dataframe tecs      
    else:
        BLUES_tec_df = adb_ldb_df.copy() #setting ADB_ldb df as the working df
        tecs = BLUES_tec_df.index.get_level_values(level=0).tolist() #dataframe tecs    
        
elif task =='ADB':
    adb_df = df.copy()
    BLUES_tec_df = adb_df.copy()
    tecs = BLUES_tec_df.index.get_level_values(level=0).tolist() #dataframe tecs    


#%%
# =============================================================================
#         GET DEMAND RELATIONS, EMMISSIONS AND LAND RELATIONS
# =============================================================================
        
#demand
demand_df = demand_conv(tpp, adb_files, year_act, count_y, ldb)    

#get the resource info from .adb files --> OK
resource_df = resource_conv(tpp, adb_files, ldb) #--> BUILD RESOURCE TECS

#relations
rel_df_no_land = Relation_constraints(tpp,adb_files, ldb )

#emissions
emission_rel_df = Relation_constraints(tpp,adb_files,  ldb, land=False, emission=True) #get emission relations

#land relations
if BLUES:
    land_rel_df = Relation_constraints(tpp,adb_files,ldb, land=True, emission=False) #get land relations
    


#%%
BT_BUP = BLUES_tec_df.copy() #backup blues df
#BLUES_tec_df = BT_BUP.copy()
#%%

    
# =============================================================================
# ADJUSTING LIFETIME AND INVESTMENT COSTS ACCORDING TO NECESSITIES
# =============================================================================
#NO LIFETIME TECS adjustment
tec_L =  [s for s in tecs if s.startswith('Land') or s.startswith('Conv')]#getting Land and Conv tecs
tec_L =  [s for s in tec_L if 'Emissions' not in s]#getting Land and Conv tecs
for nod in adb_files:
    print(nod)
    
    aux_BLUES = BLUES_tec_df.loc[idx[:, nod, :], :].copy()
    
    NO_PLL = aux_BLUES[aux_BLUES.loc[idx[:, :, :], idx['lifetime', '1']] == 'nan'].copy() #df with no pll tecs
    NO_PLL_tecs = sorted( list(set(NO_PLL.index.tolist())) )
    
    #check if it has any parameter associated with cap
    bdc_tup = aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['bound_new_capacity_up', '1']] != 'nan'].index.get_level_values(level=0).tolist()
    bdc_t = bdc_tup + aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['bound_new_capacity_lo', '1']] != 'nan'].index.get_level_values(level=0).tolist()
    
    bdi_t = aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['bound_total_capacity_up', '1']] != 'nan'].index.get_level_values(level=0).tolist()
    bdi_t = bdi_t + aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['bound_total_capacity_lo', '1']] != 'nan'].index.get_level_values(level=0).tolist()
    
    hisc_t = aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['historical_new_capacity', '1']] != 'nan'].index.get_level_values(level=0).tolist()
    
    muf_t = aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['minimum_utilization_factor', '1']] != 'nan'].index.get_level_values(level=0).tolist()
    cf_t = aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['capacity_factor', '1']] != 'nan'].index.get_level_values(level=0).tolist()

    #adding ctime into cap tec to check missing inv info
    ctime_t = aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['construction_time', '1']] != 'nan'].index.get_level_values(level=0).tolist()
    
    #inv_t = aux_BLUES[aux_BLUES.loc[idx[:,:,:], idx['investment_cost', '1']] != 'nan'].index.get_level_values(level=0).tolist()

    cap_tec = bdc_t + bdi_t + hisc_t + muf_t + cf_t + ctime_t
    cap_tec = sorted(list(set(cap_tec)) ) # TECS WITH at least one parameter associated with capacity

    if BLUES:    
        NO_PLL_tecs = [ll[0] for ll in NO_PLL_tecs if ll[0] not in tec_L]
        
        cap_tec = [ll for ll in cap_tec if ll not in tec_L]
        
    real_no_pll_tec = [t for t in NO_PLL_tecs if t not in cap_tec] #tecs with no pll and none parameter associated to capacity
    # len(real_no_pll_tec)
    no_pll_tec = [t for t in NO_PLL_tecs if t in cap_tec] #tecs with no pll and AT LEAST ONE parameter associated to capacity
    # len(no_pll_tec)
#len(real_no_pll_tec)
#len(no_pll_tec)
#len(NO_PLL_tecs)

    real_no_pll_df = aux_BLUES.loc[idx[real_no_pll_tec, :, :] ].copy() #DF with no pll tecs with none parameter associated with capacity
    no_pll_df = aux_BLUES.loc[idx[no_pll_tec, :, :] ].copy()  #DF with no pll tecs with at least one parameter associated with capacity

    #no pll tecs with no inv and no fom
    #real_no_pll_tec_no_inv_fom = list(set(real_no_pll_df[(real_no_pll_df.loc[idx[:,:,:], idx['investment_cost', '1']] == 'nan') & (real_no_pll_df.loc[idx[:,:,:], idx['fixed_cost', '1']] == 'nan')].index.get_level_values(level=0).tolist() ))
    #no pll tecs with no inv
    real_no_pll_tec_no_inv = list(set(real_no_pll_df[(real_no_pll_df.loc[idx[:,:,:], idx['investment_cost', '1']] == 'nan')].index.get_level_values(level=0).tolist() ))
    check = BLUES_tec_df.loc[idx[real_no_pll_tec_no_inv, :,:]].copy()

    #setting lifetime and inv cost for tecs with no pll no inv, no fom and none paramter associated to capacity
    BLUES_tec_df.loc[idx[real_no_pll_tec_no_inv, nod, :], idx['lifetime', '1']] = str([main_d_period])
    BLUES_tec_df.loc[idx[real_no_pll_tec_no_inv, nod, :], idx['investment_cost', '1']] = 'nan'

    # GET NO PLL WITH AT LEAST ONE CAPACITY PARAMETER AND NO INV
    NO_PLL_NO_INV = no_pll_df[no_pll_df.loc[idx[:, :, :], idx['investment_cost', '1']] == 'nan'].copy() 
    NO_PLL_NO_INV = NO_PLL_NO_INV.index.get_level_values(level=0).tolist()
    NO_PLL_NO_INV = list(set(NO_PLL_NO_INV))

    #TECS WITH NO INV AND NO LFT BUT WITH AT LEAST ONE CAP PARAMETER
    BLUES_tec_df.loc[idx[NO_PLL_NO_INV, nod, :], idx['investment_cost', '1']] = '[0.1]' #ADDING SMALL INV COST FOR ALL TECS WITH NO INV 
    BLUES_tec_df.loc[idx[NO_PLL_NO_INV  , nod, :], idx['lifetime', '1']] = str([main_d_period])

    #CHECK TECS WITH NO LIFETIME BUT WITH INV COST
    no_pll_with_inv = list(set(aux_BLUES[(aux_BLUES.loc[idx[:, nod, :], idx['lifetime', '1']] == 'nan') & (aux_BLUES.loc[idx[:, :, :], idx['investment_cost', '1']] != 'nan')].index.get_level_values(level=0).tolist()))
    BLUES_tec_df.loc[idx[no_pll_with_inv , nod, :], idx['lifetime', '1']] = '[25]'

    #CHECK TECS WITH NO LIFETIME BUT WITH FOM
    no_pll_with_fom = list(set(aux_BLUES[(aux_BLUES.loc[idx[:, :, :], idx['lifetime', '1']] == 'nan') & (aux_BLUES.loc[idx[:, :, :], idx['fixed_cost', '1']] != 'nan')].index.get_level_values(level=0).tolist()))
    no_pll_with_fom = [ll for ll in no_pll_with_fom if ll not in no_pll_with_inv] #if it has inv it is 25 lft ow it is 2
    if len(no_pll_with_fom)>0:
        BLUES_tec_df.loc[idx[no_pll_with_fom, nod, :], idx['lifetime', '1']] = str([main_d_period])

    #list(set(BLUES_tec_df[BLUES_tec_df.loc[idx[:, :, :], idx['lifetime', '1']] == 'nan'].index.get_level_values(level=0).tolist()))
    nan_inv = list(set(aux_BLUES[aux_BLUES.loc[idx[:, :, :], idx['investment_cost', '1']] == 'nan'].index.get_level_values(level=0).tolist()))
    #BLUES_tec_df.loc[idx[nan_inv, nod, :], idx['investment_cost', '1']] = 'nan'

    #adjusting BDC
    add_bdc_tecs = [t for t in tecs if t not in nan_inv] #tecs with inv
    add_bdc_tecs = [t for t in add_bdc_tecs if t not in bdc_tup] #tecs with inv and no bdc --> bdc for tecswith bdc were adjusted in the BLUES TEC DF GENERATION

    BLUES_tec_df.loc[idx[add_bdc_tecs, nod, :], idx['add_bdc', '1']] = True # set bdc fyear =0 for tecs with inv
#%%
# =============================================================================
#             #GETTING ONLY NECESSARY DATA FOR ADB
# =============================================================================
if task == 'LDB':
    pass
else:
    if not adb_full_adjustment:
        print("\nADB data was collected")
        sys.exit()
    else:
        pass

#%%
# =============================================================================
# #ADJUSTING TO GET THE DIFFERENT ELMENTS FROM EACH LAND CONV YEAR TECHNOLOGY
# =============================================================================
        
#convert into function BLUES_tec_df, tecs 
if BLUES:
    tec_L =  [s for s in tecs if s.startswith('Land') or s.startswith('Conv')]#getting Land and Conv tecs
    tec_L =  [s for s in tec_L if 'Emissions' not in s]#getting Land and Conv tecs
    
    land_data = land_yr_tec_ajd(BLUES_tec_df.iloc[:,:-3], tecs, tec_L, year_act, count_y)
     
    # =============================================================================

    
    #UPDATE COLUMNS  OF LAND_DATA AUX INTO LAND DATA DF FOR LAND/CONV TEC
    #land_data[idx['lifetime', '1']] = '[100]'
    #land_data[idx['investment_cost', '1']] = '[0.1]'
    
    # =============================================================================
    tec_L1 = [s for s in tec_L if hasNumbers(s.split(' ')[0].split('\t')[0][-2:])] #getting technologies with number in the last two characters
    tec_L2 = [s for s in tec_L1 if int(s.split(' ')[0].split('\t')[0][-2:]) >10] #keep only the values higher than 10
    #     #tec_L10 = [s for s in tec_L1 if int(s.split(' ')[0][-2:]) == 10] #keep only the values higher than 10                
    tec_la = [t for t in tecs if t not in tec_L2] #removing values higher than 10                
    # =============================================================================
    
    BLUES_tec_df = BLUES_tec_df.loc[idx[tec_la, :,:]].copy() #subsetting based on df
    land_data = land_data.loc[idx[tec_la, :,:]].copy() #subsetting based on df
        
    BLUES_tec_df.update(land_data) #UPDATE COLUMNS  OF LAND_DATA AUX INTO BLUES
    tecs_1 = BLUES_tec_df.index.get_level_values(level=0).tolist() #dataframe tecs  
    tecs_1 = sorted(list(set(tecs_1))) #removing duplicate tec values
    
    tec_L1 =  [s for s in tecs_1 if s.startswith('Land')]#getting Land 
    tec_L1 =  [s for s in tec_L1 if 'Emissions' not in s]#getting Land and Conv tecs
    BLUES_tec_df.loc[idx[tec_L1  , :, :], idx['capacity_factor', '1']] = '[1]' #'minimum_utilization_factor'
    
    tec_conv1 =  [s for s in tecs_1 if s.startswith('Conv_')]#getting Land 
    
# =============================================================================
#    #setting conv capacity equal zero in 2010 for conv tecs
#     for ttconv in tec_conv1: 
#         for Reg in set(BLUES_tec_df.loc[idx[tec_conv1  , :, :], idx['bound_total_capacity_up', '1']].index.get_level_values(level=1).tolist()):
#             auxvv = BLUES_tec_df.loc[idx[ttconv  , Reg, 1], idx['bound_total_capacity_up', '1']]
#             auxvv2 = [float(v.strip().strip("'")) for v in auxvv.strip('[]').split(',')]  
#             auxvv2[0] = 0
#             BLUES_tec_df.loc[idx[ttconv  , Reg, 1], idx['bound_total_capacity_up', '1']] = str(auxvv2)
# =============================================================================
            
    land_data.to_csv(os.path.join(path_case,project,'IX',ldb, project+'_'+task+'_'+'land_data_df'+'.csv'), sep=';') #exporting adb_df 
    if ldb_from_adb_dict:    
        btec_L =  [s for s in tecs_adb_ldb if s.startswith('Land') or s.startswith('Conv')]#getting Land and Conv tecs
        # =============================================================================
        btec_L1 = [s for s in btec_L if hasNumbers(s.split(' ')[0].split('\t')[0][-2:])] #getting technologies with number in the last two characters
        btec_L2 = [s for s in btec_L1 if int(s.split(' ')[0].split('\t')[0][-2:]) >10] #keep only the values higher than 10
        #     #tec_L10 = [s for s in tec_L1 if int(s.split(' ')[0][-2:]) == 10] #keep only the values higher than 10                
        btec_la = [t for t in tecs_adb_ldb if t not in btec_L2] #removing values higher than 10                
        # =============================================================================
        adb_ldb_df = adb_ldb_df.loc[idx[btec_la, :,:]].copy() #subsetting based on df
        btecs_1 = adb_ldb_df.index.get_level_values(level=0).tolist() #dataframe tecs  
        btecs_1 = sorted(list(set(btecs_1))) #removing duplicate tec values
#%% 
# =============================================================================
# #removing input and output columns in which all values are nan
# =============================================================================
nan_cols = BLUES_tec_df.loc[:, (BLUES_tec_df == 'nan').all()].columns.tolist()
out_inp_cols = BLUES_tec_df[['input_values', 'input_levels', 'input_commodities', 'output_values', 'output_levels', 'output_commodities']].columns.tolist()
nan_cols_oi = [cc for cc in nan_cols if cc in out_inp_cols]   
new_cols = [c for c in BLUES_tec_df.columns.tolist() if c not in nan_cols_oi] 
BLUES_tec_df = BLUES_tec_df.loc[idx[:,:,:], new_cols]
#%%
#remove all deficit tecs
#tecs_1 = [t for t in tecs_1 if not t.startswith('Deficit')]
#BLUES_tec_df = BLUES_tec_df.loc[idx[tecs_1, :, :]]
#%%
# =============================================================================
# ADJUSTING LAND USE TECS #test it for ldb
# =============================================================================
if BLUES: #if BLUES apply land functions
    BLUES_tec_df, tec_land1_new, tecs_2, land_nodes = land_adj(BLUES_tec_df, tecs_1, land_rel_df, first_land_tec, main_d_period) #this one convert names and ajust inputs and outputs
    land_generation = [t for t in tecs_2 if t.startswith('Land_Available')]
    land_classes = [t.split('_')[-1] for t in land_generation]
    if not ldb_from_adb_dict:
        tec_land1_new2 = [t for t in tec_land1_new if 'Available' not in t.split("_")[1:]]  
        
        #adjust the input of primary Land_ technlogies to match the output of land generationt tec Land_A, Land_B....
        for tt in sorted(list(set(tec_land1_new2))): #LOOPING THROUGH LAND TECS
            for n in land_nodes: #looping thorugh nodes
                BLUES_tec_df.at[idx[tt, n, 1], idx['input_commodities','1']] = '_'.join([tt.split('_')[0], tt.split('_')[-1]]).lower()       #adjust input of primary land tec to match land generation
                
    else: # if ldb_from_adb_dict is True
        adb_ldb_df, btec_land1_new, btecs_2, bland_nodes = land_adj(adb_ldb_df, btecs_1, land_rel_df, first_land_tec)
        tec_land1_new2 = [t for t in tec_land1_new if 'Available' not in t.split("_")[1:]]     
        
        #adjust the input of primary Land_ technlogies to match the output of land generationt tec Land_A, Land_B....
        for tt in sorted(list(set(tec_land1_new2))): #LOOPING THROUGH LAND TECS
            for n in land_nodes: #looping thorugh nodes
                adb_ldb_df.at[idx[tt, n, 1], idx['input_commodities','1']] = '_'.join([tt.split('_')[0], tt.split('_')[-1]]).lower()       
#%%
# =============================================================================
#                 GETTING LINK BETWEEN TECS
# =============================================================================
if ldb_from_adb_dict:
    link_t, lvl_com_t = tec_link(adb_ldb_df, btecs_2, tpp, ldb)   #GET TECHNOLOGIES CHAIN BASED ON ALL TECHNOLOGIES INPUT AND OUTPUT
else:
    link_t, lvl_com_t = tec_link(BLUES_tec_df, tecs_2, tpp, ldb)   #GET TECHNOLOGIES CHAIN BASED ON ALL TECHNOLOGIES INPUT AND OUTPUT
#%%
# =============================================================================
#     GETTING LDR INFO
# =============================================================================
if task == "ADB":    
    if has_ldr:#if database has ldr
    # get ldr info
        ldr_info = capacity_factor_conversion( tpp, adb_files, ldr_file, time1, ldb)    #ldr info function     
        ldr = ldr_tecs(ldr_info, BLUES_tec_df, tpp, adb_files, link_t, ldb) #GETTING LDR INFO
        ldr_original_tecs = sorted ( list(set( itertools.chain(*[list(ldr_info[lt].keys()) for lt in ldr_info.keys()]) ) ))
    else: #there is no ldr
        ldr_info = {Reg+'_cf_norm': {} for Reg in adb_files}
        ldr = [[] for a in range(5)] #empty ldr list
        ldr_original_tecs = []
    ldr = list(ldr)
    ldr.append(ldr_original_tecs)    
    ldr = tuple(ldr)
    #exporting ldr list
    with open(os.path.join(path_case, project,'IX',ldb,  project+"_ADB_ldr"), 'wb') as fp: #save list
        pickle.dump([ldr], fp) #save ldr list

else:
    if ldb_ldr:      #if ldr is differnt in ldb
        ldr_info = capacity_factor_conversion(tpp, adb_files, ldr_file, time1, ldb) #GETTING LDR INFO
        ldr = ldr_tecs(ldr_info, BLUES_tec_df, tpp, adb_files, link_t, ldb) #GETTING LDR INFO
        ldr_original_tecs = sorted ( list(set( itertools.chain(*[list(ldr_info[lt].keys()) for lt in ldr_info.keys()]) ) ))
        
    else: #if ldb does not change ldr info and/or technologies --> load ldr info from adb
        ldr_info = capacity_factor_conversion( ('ADB', path_case, project) , adb_files, ldr_file, time1, ldb)    #ldr info function  
        ldr_original_tecs = sorted ( list(set( itertools.chain(*[list(ldr_info[lt].keys()) for lt in ldr_info.keys()]) ) ))
        
        if ldb_from_adb_dict:
            ldr = ldr_tecs(ldr_info, adb_ldb_df, ('ADB', path_case, project), adb_files, link_t, ldb) #GETTING LDR INFO        
        else:  
            ldr = ldr_tecs(ldr_info, BLUES_tec_df, ('ADB', path_case, project), adb_files, link_t, ldb) #GETTING LDR INFO
            
    ldr = list(ldr)
    ldr.append(ldr_original_tecs)    
    ldr = tuple(ldr)
    with open(os.path.join(path_case, project, ldb, project+"_LDB_ldr_"+ldb), 'wb') as fp: #save list
        pickle.dump([ldr], fp) #save ldr list

np.save(os.path.join(path_case, project,'IX', ldb, project + '_'+task+"_ldr_info"+ldb_aux+'.npy'),ldr_info) #saving dictionary
#%%
# =============================================================================
# #balance equality commodity level
# =============================================================================
balance_eq_df = bal_eq(tpp, adb_files, ldb)
#%%
# =============================================================================
# APLLYING LAND FUNCTIONS TO GET LAND PARAMETERS FOR MESSAGE IX
# =============================================================================

if BLUES: #if BLUES apply land functions
    land_chain = land_tec_link(link_t, land_generation) #this one get the link between land tecs
# =============================================================================
#     if task ==  'LDB':
#         adb_ldb_df, tec_land1_new, tecs_2, land_nodes = land_adj(adb_ldb_df, tecs_1, land_rel_df, first_land_tec) #this one convert names and ajust inputs and outputs        
#     else: #if task == 'ADB'
# =============================================================================


    cat_land_tec_df, map_land_rel_df, sec_land_tec, land_tec, class_land_tec = land_ix_adj(land_nodes, tec_land1_new, land_chain) #organize the data according to the ix pattern
    
    land_ix_data = cat_land_tec_df, map_land_rel_df, sec_land_tec, land_tec, class_land_tec, land_generation
    
    with open(os.path.join(path_case, project, 'IX', ldb,project+'_' + task + '_land_ix_data' + ldb_aux), 'wb') as fp: #save list
        pickle.dump(land_ix_data, fp)
#%%
# =============================================================================
# TRY TO MAKE IT BETTER ADD OLD TECS INTO BLUES TEC DF =- REMOVE THE ORGINAL TECS AHD ADJUST BUILDING TEC FUNCTION
# =============================================================================
# =============================================================================
# IT IS NOT BEING USED ANYMORE
# =============================================================================
#adding dummy old technologies into the list of tecs
if old_t:
    print("with old tecs")
    print("it needs to be improved see version _22") # bda it not being accounted for old_tec e.g.:oil_postsalt_2000 
    sys.exit()
else:  # with no old tecs
    tn_old_out = []
    tn_old_out1 = [t for t in tecs_2 if 'existing' in t.lower() or 'exist' in t.lower()] #check if it is an existing technology
    for tn in tn_old_out1 :
        nodes = list(set(BLUES_tec_df.loc[tn].index.get_level_values('node').values.tolist())) #getting nodes of the technology
        for nod in nodes:    #looping through nodes
            if BLUES_tec_df.loc[idx[tn, [nod], :], ('historical_new_capacity_years')].values.tolist()[0][0] != 'nan':
                tn_old_out.append(tn)
    
    tn_old_out = sorted(list(set( tn_old_out )))
            
    old_tec_list = [[], tn_old_out, []]      
    tecs_n = tecs_2[:]            #adding old tecs that represent technologies installed in 2005 or earlier
    
    with open(os.path.join(path_case, project,'IX',ldb, project+'_'+task+'_'+'old_tec_list'+ldb_aux), 'wb') as fp: #save list
        pickle.dump(old_tec_list, fp)    
        
#%%
# =============================================================================
# GET EMISSION AND RELATION BOUNDARIES AND COSTS
# =============================================================================
#RELATION AND EMISSION BOUND VALUES
relation_upper, relation_lower = relation_values(tpp, rel_df_no_land, year_act, count_y, ldb) #getting relations boundaries
emission_bound = relation_values(tpp, emission_rel_df, year_act, count_y, ldb, emission_values=True) #getting emissions boundaries
emission_bound['type_tec'] = 'all'

#RELATION AND EMISSION COSTS
relation_cost = relation_cost_calculaion(tpp, rel_df_no_land, year_act, count_y, ldb) #getting relations boundaries
emission_cost = relation_cost_calculaion(tpp, emission_rel_df, year_act, count_y, ldb, emission_values=True) #getting emissions boundaries
emission_cost['type_tec'] = 'all'
#%%
# =============================================================================
# CHANGE WITHDRAW COMMODITY TO HAVE A DIFFERENT NAME THAN WITHDRAW TEC
# =============================================================================
if BLUES:
    kilen = list(range(1,1+len(BLUES_tec_df['input_commodities'].columns.values.tolist()) ))
    kolen = list(range(1,1+len(BLUES_tec_df['output_commodities'].columns.values.tolist()) ))
    #adjust withdraw commodity have a different name of withdraw tec 
    for i in kilen:
        
        BLUES_tec_df.loc[idx[:, :, :], idx['input_commodities', str(i)]] = np.where(BLUES_tec_df.loc[idx[:, :, :], idx['input_commodities', str(i)]].str.startswith('withdrawn'), 
                                                                                              'water_' + BLUES_tec_df.loc[idx[:, :, :], idx['input_commodities', str(i)]],
                                                                                              BLUES_tec_df.loc[idx[:, :, :], idx['input_commodities', str(i)]])
        
        
    for i in kolen:
        BLUES_tec_df.loc[idx[:, :, :], idx['output_commodities', str(i)]] = np.where(BLUES_tec_df.loc[idx[:, :, :], idx['output_commodities', str(i)]].str.startswith('withdrawn'), 
                                                                                              'water_' +BLUES_tec_df.loc[idx[:, :, :], idx['output_commodities', str(i)]],
                                                                                              BLUES_tec_df.loc[idx[:, :, :], idx['output_commodities', str(i)]])

    balance_eq_df['commodity'] =  np.where(balance_eq_df['commodity'].str.startswith('withdrawn'), 
                                                                                              'water_' +balance_eq_df['commodity'],
                                                                                              balance_eq_df['commodity'])

    balance_eq_df.to_csv(os.path.join(path_case,project,'IX', ldb,project+'_'+task+ldb_aux+'_bal_eq_df'+'.csv'), sep =';', index=False)   
#%%
# =============================================================================
#     ADJUST LAND TEC LIFETIME
# =============================================================================
if BLUES:
    BLUES_tec_df.loc[idx[land_tec,:,:], idx['lifetime','1']] = '[100]'
    BLUES_tec_df.loc[idx[land_tec,:,:], idx['investment_cost','1']] = '[0.1]'
    
    conv_tecs = [t for t in tecs_n if t.startswith('Conv_')]
    
    BLUES_tec_df.loc[idx[conv_tecs,:,:], idx['lifetime','1']] = str([main_d_period])
    #BLUES_tec_df.loc[idx[conv_tecs,:,:], idx['investment_cost','1']] = '[1]'    
    tec_L1 =  [s for s in tecs_n if s.startswith('Land')]#getting Land 
    tec_L1 =  [s for s in tec_L1 if 'Emissions' not in s]#getting Land and Conv tecs
    tec_L1 =  [s for s in tec_L1 if 'Land_Available' not in s]#
    

#%%    
# =============================================================================
#     ADJUSTING TECS ANG GET LINKS BETWEEN TECS AGAIN
# =============================================================================
tecs_f = sorted(list(set(tecs_n))) #removing duplicate

link_t2, lvl_com_t2 = tec_link(BLUES_tec_df, tecs_f, tpp, ldb)   #GET TECHNOLOGIES CHAIN BASED ON ALL TECHNOLOGIES INPUT AND OUTPUT

    
with open(os.path.join(path_case, project,'IX',ldb, project+'_'+task+'_tecs_f'+ldb_aux), 'wb') as fp: #save list
    pickle.dump(tecs_f, fp)     
#%%
# =============================================================================
# ADJUSTING LAND USE EMISSION
# =============================================================================
#if 'Landuse_Emissions' in tecs_f:
#%%
# =============================================================================
# GENERATING DICTIONAY THAT AGGREGATE ALL INFO
# =============================================================================
    
BLUES_tec_df.to_csv(os.path.join(path_case,project,'IX', ldb,project+'_'+task+ldb_aux+'_df_'+VV+'.csv'), sep =';') #exporting adb_df

if ldb_from_adb_dict:
# =============================================================================
#     REVIEW LDB BUILDING TECS
# =============================================================================
    BLUES_tec_df = BLUES_tec_df.loc[~(BLUES_tec_df=='nan').all(axis=1)] #removing rows in which all valuea are nan
    tecs_f = sorted(list(BLUES_tec_df.index.get_level_values(level=0).tolist() )) #dataframe tecs
    tecs_f = sorted(list(set(tecs_f) ) )
    tecs_f1 = [tt for tt in tecs_f if 'Exist' not in tt] #removing existing
    
    BLUES_tec_df.to_csv(os.path.join(path_case,project,'IX', ldb,project+'_only_'+task+ldb_aux+'_df'+'.csv'), sep=';') #save ldb dr
    
    adb_tec_dict = np.load(os.path.join(path_case,project,'IX',project + "_ADB_tec_dict_ix_format.npy"), allow_pickle=True).item() #loading tec idct from adb
    
    new_ldb_tec = [t for t in tecs_ldb if t not in tecs_adb] #getting new ldb tecs

    new_warning_list, new_tec_dict, not_def_rel = build_adb_dfs(tpp, ldb, BLUES_tec_df, new_ldb_tec, all_time_info, old_tec_list, ldr, ldr_info, land_rel_df, rel_df_no_land,BLUES)
    #if it work

# =============================================================================
    with open(os.path.join(path_case, project,'IX',ldb, project+'_'+task+"_ONLY_warning_list"+ldb_aux), 'wb') as fp: #save list
        pickle.dump(new_warning_list, fp) #save ldr list              
    
    np.save(os.path.join(path_case, project,'IX', ldb, project + '_'+task+"_ONLY_tec_dict_ix_format"+ldb_aux+'.npy'),new_tec_dict)
# =============================================================================

    removed_adb_tecs = [t for t in tecs_adb if t not in tecs_ldb] #removing tecs from adb dict bc they are not considered in ldb tecs
    for r in removed_adb_tecs: #LOOPING THROUGH TECS THAT MUST BE REMOVED
        if adb_tec_dict.has_key(r): #CHECK IF DICT CONTAINS R TEC
            del adb_tec_dict[r] #REMOVE
            
    ntec_dict = copy.deepcopy(adb_tec_dict)#.copy()     #APPENDING NEW TEC INTO ADB DICT 
    ntec_dict.update(new_tec_dict)


else:    #it is building tecs based on the whole BLUES_tec_df
    
    warning_list, tec_dict, not_def_rel = build_adb_dfs(tpp, ldb, BLUES_tec_df, tecs_f, all_time_info, old_tec_list, ldr, ldr_info, land_rel_df, rel_df_no_land,BLUES)
    

#%%
with open(os.path.join(path_case, project,'IX',ldb, project+'_'+task+"_warning_list"+ldb_aux+'_'+VV), 'wb') as fp: #save list
    pickle.dump(warning_list, fp) #save ldr list              

np.save(os.path.join(path_case, project,'IX', ldb, project + '_'+task+"_tec_dict_ix_format"+ldb_aux+'_'+VV+'.npy'),tec_dict)


tecs_f2 = sorted( list(tec_dict.keys()) )
#%%
tec_inp_out = {}
for tt in tecs_f2:
    tec_inp_out[tt] = {}
    inpl = list(set([(l,c) for l,c in zip(tec_dict[tt]['input']['level'], tec_dict[tt]['input']['commodity'])]))
    outl = list(set([(l,c) for l,c in zip(tec_dict[tt]['output']['level'], tec_dict[tt]['output']['commodity'])]))
    tec_inp_out[tt]['input'] = inpl
    tec_inp_out[tt]['output'] = outl
    
np.save(os.path.join(path_case, project,'IX', ldb, project + '_'+task+"_tec_inp-out"+ldb_aux+'_'+VV+'.npy'),tec_inp_out)
#%%
dt_finalf =  datetime.now() - startTimeg
print(task+ldb_aux+' Done.')
print('\nRuntime load generation script:    '+ str(dt_finalf)) 
