# AMSTERDAM FOOT MODEL KINEMATICS
# Wouter Van Caekenberghe, march 2023 

import numpy as np
import pandas as pd
import os
from Static import Static_AFM
from Dynamic_ext import Dynamic_AFM
from plot_ext import plot_AFM_R, plot_AFM_L
from plot_group_ext import plot_group_AFM_R, plot_group_AFM_L
from plot_conditions_ext import plot_conditions_AFM_R, plot_conditions_AFM_L
from Trial_info import trial_paths, trial_gaitcycle, excel_gaitinfo
from combine import avg_groupdata_right, avg_groupdata_left
from df_to_excel import *

############### INPUT VARIABLES ################################################################################################
# Define directory
directory = (r'C:\Users\WVanCaeckenberge\Desktop\AFM_final\Wedge') # c3d files
directory_gs = (r'C:\Users\WVanCaeckenberge\.spyder-py3\Gait_Events\Wedge_AFM.xlsx') # Excel file containing gait cycles 
# What side do you want to process 
footside = 2 #(1=right, 2=left, 3=both)
forefootanalysis = 2 #(1=full forefoot, 2=full forefoot + lat/med seperately)
# What data do you want
indivplots = 1 #(1=individual plots for every trial)
indivexcel = 0 #(1=Excel file with data for individual trials)
groupedplots = 0 #(1=grouped plots per condition)
conditionplots = 0 #(1=plot conditions together)
conditionexcel = 0 #(1=Excel file with data per condition)
normalize = 0 #(1=normalize)
# Please adjust marker names Markernames.py
#################################################################################################################################

if indivexcel==1 or conditionexcel==1:
    wd = os.getcwd()
    path = os.path.join(wd, 'AFM_Kinematics','Excel_data')
    if not os.path.isdir(path):
        os.makedirs(path)

Static_paths, Dynamic_paths, conditions, filename = trial_paths(directory)

gaitinfo_right, gaitinfo_left = excel_gaitinfo(directory_gs)

# Process and plot data
groupdata = dict()
count = -1
for idx, c in enumerate(conditions):
    tempr = list()
    templ = list()
    for idxx, trial in enumerate(Dynamic_paths[c]):
        count += 1 
        data_path_static = Static_paths[c]
        data_path_dynamic = trial
        trialname = filename[idx][idxx]
        print('Processing condition: ' + c + ', Trial: ' + trialname)  
        # Right foot
        if footside==1 or footside==3:
            gc = trial_gaitcycle(count, gaitinfo_right)
            staticangles_right, staticinfo = Static_AFM(data_path_static,1,forefootanalysis)
            dynamicangles_right = Dynamic_AFM(data_path_dynamic,staticinfo,gc,1,forefootanalysis,trialname,c)
            
            # Save data for later
            tempr.append(dynamicangles_right)
            anglenames = list(dynamicangles_right['gaitcycle'].keys())

            if indivexcel == 1:
                indiv_data_to_excel_R(wd, dynamicangles_right, staticangles_right, trialname, c)

        # Left foot
        if footside==2 or footside==3:
            gc = trial_gaitcycle(count, gaitinfo_left)
            staticangles_left, staticinfo = Static_AFM(data_path_static,2,forefootanalysis)
            dynamicangles_left = Dynamic_AFM(data_path_dynamic,staticinfo,gc,2,forefootanalysis,trialname,c)
            
            # Save data for later
            templ.append(dynamicangles_left)
            anglenames = list(dynamicangles_left['gaitcycle'].keys())

            if indivexcel==1:
                indiv_data_to_excel_L(wd, dynamicangles_left, staticangles_left, trialname, c)

        # Plot data
        if indivplots==1:
            if footside==1 or footside==3:
                plot_AFM_R(dynamicangles_right, forefootanalysis, normalize)
                

            if footside==2 or footside==3:
                plot_AFM_L(dynamicangles_left, forefootanalysis, normalize)
                
    if footside==1 or footside==3:   
        groupdata[c + '_right'] = tempr
    if footside==2 or footside==3:
        groupdata[c + '_left'] = templ

# Average data per condition 
if footside==1 or footside==3:
    avgdata_right = avg_groupdata_right(groupdata, conditions, anglenames)
    if conditionexcel==1:
        avg_data_to_excel_R(avgdata_right, conditions, wd)

if footside==2 or footside==3:
    avgdata_left = avg_groupdata_left(groupdata, conditions, anglenames)
    if conditionexcel==2:
        avg_data_to_excel_L(avgdata_left, conditions, wd)

# Plot grouped data
if groupedplots==1:
    for idx, c in enumerate(conditions):
        if footside==1 or footside==3:
            xx = groupdata[c + '_right']
            plot_group_AFM_R(xx, forefootanalysis, normalize, c, avgdata_right, filename[idx])

        if footside==2 or footside==3:
            xx = groupdata[c + '_left']
            plot_group_AFM_L(xx, forefootanalysis, normalize, c, avgdata_left, filename[idx])

# Plot conditions 
if conditionplots==1:
    if footside==1 or footside==3:
        plot_conditions_AFM_R(avgdata_right, forefootanalysis, conditions)
    
    if footside==2 or footside==3:
        plot_conditions_AFM_L(avgdata_left, forefootanalysis, conditions)

print('Finished')