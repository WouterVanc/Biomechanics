# Detect Gait Events: BigFootLab, Paal
# Wouter Van Caekenberghe, 4th of may, 2023

# The script currently only takes the first gaitcycle for each foot, if you wish otherwise, adjust line 111-114 (note: not every trial has 2 gaitcycles per foot)

import pandas as pd
from DGE_BFL import detect_gait_events
from Trial_info import trial_paths
import os

# Please make sure the directory(data) listed below, is the same directory(data) used as input for the 'AFM_script_extended'
directory = (r'C:\Users\WVanCaeckenberge\Desktop\AFM_final\Wedge')
# How do you want your output excel to be named:
excelfile_name = 'Wedge_AFM'
# Enter name of the heel markers that were used
lhm = 'LCALD'
rhm = 'RCALD'

# Create folders
wd = os.getcwd()
path = os.path.join(wd, 'Gait_Events')
path2 = os.path.join(path, f'Visualizations_{excelfile_name}')
if not os.path.isdir(path):
    os.makedirs(path)
if not os.path.isdir(path2):
    os.makedirs(path2)

# Find trial paths
Static_paths, Dynamic_paths, conditions, filename = trial_paths(directory)

# Detect gait events and write to excel
detect_gait_events(Dynamic_paths, conditions, wd, excelfile_name, filename, rhm, lhm)

print('Finished')
