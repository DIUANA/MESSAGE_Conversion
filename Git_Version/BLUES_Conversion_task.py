# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 08:55:02 2020

@author: FD
"""
"""

    If it's the first time you are working on this conversion you sould set task = 'ADB' and new_df = 'True'
     The idea is to get all adb database
        then start to work on the ldb based on the original adb like in the old Message/Blues framework 
        
        You can also set task == LDB to work on ldb database or load the data that have been generated before
        
        OBS: This conversion is defind to get only one historical year 
        message ix still presents weird behaviors if more than 1 historical year is set
        I still need to analyze it better
    
"""
#%% defining exogenous parameters.........................                 
#odir = dir()
model_name = 'BLUES' #name of your model
project = 'Brasil_1.99'# 'LB_NO_LDR'# #set the name of your project 


adb_files = ['Brasil', 'CO',  'N', 'NE', 'S','SE'] #set regions names/folders/files | it should match the name of the region folders in the original MESSAGE

task = 'ADB' #set task as 'LDB' or 'ADB' -->LDB version is not working at the moment

main_region = 'Brasil' #it should be one of the regions in adb_files # it is not working for   multi parent regions

new_df = False  #if you would like to generate a new df from adb or ldb files
adb_full_adjustment = True #Set True if you wish to do all necessary adjustments to build a dictionary from ADB ow set False (Useful to load only adb data that will be complemented by ldb data)


time1 = 'original' #it must be set as  'original' [according to the original model] or 'months' [1-12], 'month-hours' [24 hours for each month which means 24*12=288 time steps] depending on model original load regions

has_ldr = False  # Set as True if your database has ldr and False if it is not or if it must be ignored
ldr_file = False  # set True if you are using a .ldr file to get load curve info and False if you are using it direct from .adb/.ldb file



first_model_year = 2010 #first year that should be modelled (included in years_act) #values should be equal or higher than 2010 for BLUES
emissions_constraints = ['CH4L', 'GHG', 'CO2P', 'CH4', 'CO2', 'CO2L', 'N2OL','CH4', 'N2O', 'TCO2', 'cGHG'] #setting relation that are emissions

BLUES = True     # Set True IF THE CONVERSION IS FOR THE BLUES MODEL
first_land_tec = 'Land_Bal_Total' #first land tec in BLUES model || It is not necessary if it is not BLUES model (if BLUES is set as False)


Load_info = True #SET AS True IF YOU WOULD LIKE to load info that have been genrerated before




file_version = 'last' #define if it is going to be based on the 'last' file version or any other

old_t = False #it must be set as False for now --> still need to work on this to make True available

period_interp = 'progressive' #  'recursive' # how the conversion is going to consider year time-steps





#%% importing packages

#required packages
from datetime import datetime
import matplotlib.pyplot as plt
import os                                  
import pandas as pd
import numpy as np
import re
import shutil
from itertools import cycle
import itertools                                     
import string
import statistics as stat
import sys
from functools import reduce
import pickle
import collections
import copy
#os.chdir(path)

path_case = os.path.join(os.getcwd(),'Cases') #define the path in which there is all the data/functions for conversion
path_script = os.path.join(os.getcwd(),"Scripts")
path_model = os.path.join("\\".join(os.getcwd().split('\\')[:-1]),"model")
idx = pd.IndexSlice #pandas option to slice list
pd.set_option('display.max_columns', 10)

#%% check if its land
if New_Land:
    sp_opt_n = '_'+sp_opt
    if Norm:
        sp_opt_n = sp_opt_n+'_Norm'
else:
    sp_opt_n = sp_opt
    
#%% new value name
 
lf = [ll for ll in os.listdir(os.path.join(path_model, 'data')) if 'MsgData_'+model_name+'_'+project+'_IX' in ll]
lf2 = [ ll for ll in lf if ll.split('IX_')[-1].split('_')[0].split('.')[0].isnumeric() ]
lf3 = [int(ll.split('IX_')[-1].split('_')[0].split('.')[0]) for ll in lf2]
if len(lf3) >0:
    
    if New_Land:
        VV = str(max(lf3))    
    elif climate_case:
        VV = str(max(lf3)) 
    else:
        VV = str(max(lf3) +1)
else:
    VV ='01'
    
if len(VV) <=1:
    VV = '0'+VV
    


sc_name = project+'_IX_'+VV



#%%    generating or loading parameters and dict
if not Load_info:
    file = 'BLUES_globals_generation_r62.py'
    global_file = os.path.join(path_script, file) #I have previously developed a global ful that was loading all data, the short has more functions
    try:
        exec(compile(open(global_file).read(), global_file, 'exec'))  
    except SystemExit:
        print("Finishing: "+file)
    #print("Continuing")

    #CHECK DATA THAT MUST BE LOADED -->Make it better 
else:
    print('Loading info...')
    if file_version == 'last':
        files = os.listdir(os.path.join(path_case, project, 'IX'))
        files = [ll for ll in files if not os.path.isfile(os.path.join(path_case, project, 'IX', ll))]    
        file_version = str(max([int(ll) for ll in files]))
    file = 'BLUES_globals_load_data_r14.py'
    global_file = os.path.join(path_script, file ) #I have previously developed a global ful that was loading all data, the short has more functions
    try:
        exec(compile(open(global_file).read(), global_file, 'exec'))  
    except: #make the warning better
        raise Exception("It was not possible to load: \n\n\t\t"+ str(Load_variables[k]) + "\n\n\tCheck if it is the correct path")
    #print("Continuing")
    

#%%    running set par def
set_file = 'BLUES_setpardef_r125_minor_changes.py'
exec(compile(open(os.path.join(path_script,set_file )).read(), set_file, 'exec'))  
#exec(compile(open(os.path.join(os.getcwd(), 'BLUES_setpardef_r55.py')).read(), 'BLUES_setpardef_r55.py', 'exec'))  


#%% moving files
# =============================================================================
# Moving files
# =============================================================================
if not Load_info:

    try:
        # Create target Directory
        os.mkdir(os.path.join(path_case, project, 'IX', VV))
        print( "IX Directory " ,  " Created ") 
    except FileExistsError:
        print("IX Directory "+VV ,  " already exists")
        fvv = os.listdir(os.path.join(path_case, project, 'IX', VV))
        for fll in fvv:
            os.remove(os.path.join(path_case, project, 'IX', VV, fll))
        
    files = os.listdir(os.path.join(path_case, project, 'IX'))
    files = [ll for ll in files if os.path.isfile(os.path.join(path_case, project, 'IX', ll))]
    
    for ff in files: #move files
        shutil.move(os.path.join(path_case, project, 'IX', ff), os.path.join(path_case, project, 'IX',VV, ff))
        
else: #if Laof_info copy and rename

    try:
        # Create target Directory
        os.mkdir(os.path.join(path_case, project, 'IX', VV))
        print( "IX Directory " ,  " Created ") 
    except FileExistsError:
        print("IX Directory "+VV ,  " already exists")
    
    files = os.listdir(os.path.join(path_case, project, 'IX', file_version))
    files = [ll for ll in files if os.path.isfile(os.path.join(path_case, project, 'IX',file_version, ll))]

    for ff in files:
        shutil.copy(os.path.join(path_case, project, 'IX', file_version, ff), os.path.join(path_case, project, 'IX',VV, ff))    
        if ff.split(project)[-1].split('.')[0][-2:] == file_version:
            nff = ff.replace(file_version, VV)
            os.rename(os.path.join(path_case, project, 'IX',VV, ff), os.path.join(path_case, project, 'IX',VV, nff))

#%% running ixmp file
#mp.close_db()    
ixmp_file = 'BLUES_ixmp_r24.py'

exec(compile(open(os.path.join( path_script, ixmp_file)).read(), ixmp_file, 'exec'))  

file_version = VV
#print(nn + ' is done!')
# =============================================================================
#     import gc
#     gc.collect()
#     for name in dir():
#         if name not in odir:
#             del globals()[name]
# 
#     for name in dir():
#         if name not in odir:
#             del locals()[name]
# =============================================================================

#%%

#check if building functions are wroking
# =============================================================================
# #try to make them better. They are demanding many parameters
# ============================================================================= 

