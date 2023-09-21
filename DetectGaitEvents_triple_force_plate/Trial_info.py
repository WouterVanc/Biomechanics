# Find paths to all trials
# directory should look like this:
# Main folder = input path
# subfolders are named to each conditions (i.e. barefoot, shoe, shoe+insoles, ...)
# Each of those folders should contain a 'Static' and a 'Dynamic' folder
# Static folder contains one static file for specific condition
# Dynamic folder contains all trials for specific condition

import os 
import pandas as pd
import numpy as np
                                                                                                                                
def trial_paths(main_folder_path):
    conditions = os.listdir(main_folder_path)
    
    # Find paths to all static and dynamic files in specififed folder
    Static_paths = dict()
    for v in conditions:  
        temp = os.listdir(os.path.join(main_folder_path,v,'Static'))
        Static_paths[v] = os.path.join(main_folder_path,v,'Static',temp[0])
        
    Dynamic_paths = dict()
    filename = list()
    for v in conditions:
        temp = os.listdir(os.path.join(main_folder_path,v,'Dynamic'))
        filename.append(os.listdir(os.path.join(main_folder_path,v,'Dynamic')))
        temp2 = list()
        for n in temp:
            temp2.append(os.path.join(main_folder_path,v,'Dynamic',n))
        Dynamic_paths[v] = temp2

    return Static_paths, Dynamic_paths, conditions, filename

def excel_gaitinfo(excel_path):
    xls = pd.ExcelFile(excel_path)
    gaitinfo_right = pd.read_excel(xls, sheet_name='Right')
    gaitinfo_left = pd.read_excel(xls, sheet_name='Left')
    
    return gaitinfo_right, gaitinfo_left

def trial_gaitcycle(count, gaitinfo):
    hs = gaitinfo['HS'].iloc[count]
    to = gaitinfo['TO'].iloc[count]

    gaitcycle = np.array([hs, to])

    return gaitcycle 