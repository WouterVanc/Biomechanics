# Final AFM script
# Wouter Van Caekenberghe, march 2023

# forefootanalysis niet gelijk aan 2 eens nachecken of de plotting nog werkt
import numpy as np
from Static import Static_AFM
from Dynamic import Dynamic_AFM
from plot import plot_AFM_R

############### INPUT VARIABLES ################################################################################################
# Define path
data_path_static = (r'C:\Users\WVanCaeckenberge\Desktop\MOCAP_data\Reference_barefoot\Static.c3d') # C3D-file static trial
data_path_dynamic = (r'C:\Users\WVanCaeckenberge\Desktop\MOCAP_data\Reference_barefoot\GaitTrial.c3d') # C3D-file dynamic trial
# What side do you want to process 
footside = 1 #(1=right, 2=left, 3=both)
forefootanalysis = 2 #(1=full forefoot, 2=full forefoot + lat/med seperately)
# Gaitcycles
startframe = 168 # First frame of C3D file (if trimmed)
HS_to_HS = np.array([411,524]) # Heelstrike to heelstrike frames i.e. one full gaitcycle to be analyzed
# Do you want to normalize gaitcycle to 101 data points
normalize = 1 #(1=normalize)
#################################################################################################################################

# Right foot
if footside==1 or footside==3: 
    staticangles_right, staticinfo = Static_AFM(data_path_static,1,forefootanalysis)
    dynamicangles_right = Dynamic_AFM(data_path_dynamic,staticinfo,HS_to_HS,startframe,footside,forefootanalysis)

# Left foot
if footside==2:
    staticangles_left, staticinfo = Static_AFM(data_path_static,2,forefootanalysis)
    dynamicangles_left = Dynamic_AFM(data_path_dynamic,staticinfo,HS_to_HS,startframe,footside,forefootanalysis)

# Plot data
if footside==1 or footside==3:
    plot_AFM_R(dynamicangles_right, forefootanalysis, normalize)

if footside==2 or footside==3:
    plot_AFM_R(dynamicangles_right, forefootanalysis, normalize)
