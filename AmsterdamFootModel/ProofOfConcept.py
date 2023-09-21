# Amsterdal Foot Model kinematics
# Wouter Van Caekenberghe, 14/03/2023
# Import packages
from ezc3d import c3d
import statistics as st
import sys
import numpy as np
import math 
from datafilter import filter_withnan
from normalize import norm_to_101
import matplotlib.pyplot as plt
# Coordinates system functions
import Shank_cs 
import Hindfoot_cs
import Midfoot_cs
import Forefoot_cs
import Hallux_cs
import Foot_cs
from RDecompose import RDecompose
from TTA import TTA
from Projection import projectonplane

############### INPUT VARIABLES ################################################################################################
# Define path
data_path_static = (r'C:\Users\WVanCaeckenberge\Desktop\MOCAP_data\Reference_barefoot\Static.c3d') # C3D-file static trial
data_path_dynamic = (r'C:\Users\WVanCaeckenberge\Desktop\MOCAP_data\Reference_barefoot\GaitTrial.c3d') # C3D-file dynamic trial 
# What side do you want to process 
footside = 2 #(1=right, 2=left, 3=both)
forefootanalysis = 1 #(1=full forefoot, 2=full forefoot + lat/med seperately)
# Gaitcycles
startframe = 168 # First frame of C3D file (if trimmed)
HS_to_HS = np.array([411,524]) # Heelstrike to heelstrike frames i.e. one full gaitcycle to be analyzed
# Do you want to normalize gaitcycle to 101 data points
normalize = 1 #(1=normalize)
#################################################################################################################################

#### Load in data
c = c3d(data_path_static)
c_dyn = c3d(data_path_dynamic)

# Extract variables static trial
point_data = c['data']['points']
labels = list(c['parameters']['POINT']['LABELS']['value'])
points_residuals = c['data']['meta_points']['residuals']
analog_data = c['data']['analogs']
framerate = int(c['header']['points']['frame_rate'])
lastframe = int(c['header']['points']['last_frame'])

# Extract variables dynamic trials
point_data_dyn= c_dyn['data']['points']
labels_dyn = list(c_dyn['parameters']['POINT']['LABELS']['value'])
points_residuals_dyn = c_dyn['data']['meta_points']['residuals']
analog_data_dyn = c_dyn['data']['analogs']
framerate_dyn = int(c_dyn['header']['points']['frame_rate'])
lastframe_dyn = int(c_dyn['header']['points']['last_frame'])

# Create time list
t = list(range(0,lastframe_dyn+1))
t = [x/100 for x in t]

# Low pass filter data (4th order, 6hz cut off, butterworth) leaving NaNs in place
point_data = filter_withnan(point_data, framerate)
point_data_dyn = filter_withnan(point_data_dyn, framerate_dyn)

# Allocate data to marker names for static trial
m = dict()
for ind, name in enumerate(labels):
    m[name + 'x'] = point_data[0,ind,:]
    m[name + 'y'] = point_data[1,ind,:]
    m[name + 'z'] = point_data[2,ind,:]

# Allocate data to marker names for dynamic trial
m_dyn = dict()
for ind, name in enumerate(labels_dyn):
    m_dyn[name + 'x'] = point_data_dyn[0,ind,:]
    m_dyn[name + 'y'] = point_data_dyn[1,ind,:]
    m_dyn[name + 'z'] = point_data_dyn[2,ind,:]

#### AFM STATIC FUNCTION ####
#### Static Coordinate System
nanmatrix = np.array([float('nan')]*9).reshape(3,3)
# SHANK
if 'ASHN' in labels and 'LSHN' in labels and 'TT' in labels and 'FH' in labels and 'LM' in labels and 'MM' in labels:
    # Mean x,y,z coordinate per marker
    ASHN = np.array([st.mean(m['ASHNx']), st.mean(m['ASHNy']), st.mean(m['ASHNz'])])
    LSHN = np.array([st.mean(m['LSHNx']), st.mean(m['LSHNy']), st.mean(m['LSHNz'])])
    FH = np.array([st.mean(m['FHx']), st.mean(m['FHy']), st.mean(m['FHz'])])
    TT = np.array([st.mean(m['TTx']), st.mean(m['TTy']), st.mean(m['TTz'])])
    LM = np.array([st.mean(m['LMx']), st.mean(m['LMy']), st.mean(m['LMz'])])
    MM = np.array([st.mean(m['MMx']), st.mean(m['MMy']), st.mean(m['MMz'])])
    
    Shank_anat, Oshanka = Shank_cs.anat(ASHN,LSHN,FH,TT,LM,MM,footside)
    Shank_tech, Oshankt = Shank_cs.tech(ASHN, LSHN, TT)
                              
else:
    sys.exit('ERROR: marker name not found in static file [shank]')
    
# Hindfoot
if 'CALD' in labels and 'CALP' in labels and 'ST' in labels and 'PT' in labels:
    CALD = np.array([st.mean(m['CALDx']), st.mean(m['CALDy']), st.mean(m['CALDz'])])
    CALP = np.array([st.mean(m['CALPx']), st.mean(m['CALPy']), st.mean(m['CALPz'])])
    ST = np.array([st.mean(m['STx']), st.mean(m['STy']), st.mean(m['STz'])])
    PT = np.array([st.mean(m['PTx']), st.mean(m['PTy']), st.mean(m['PTz'])])
    
    Hindfoot_anat, Ohindfoota = Hindfoot_cs.anat(CALD, CALP, ST, PT)
    Hindfoot_tech, Ohindfoott = Hindfoot_cs.tech(CALD, ST, PT)

else:
    sys.exit('ERROR: marker name not found in static file [hindfoot]')

# Midfoot
if 'NAV' in labels and 'BM2' in labels and 'BM5' in labels:
    NAV = np.array([st.mean(m['NAVx']), st.mean(m['NAVy']), st.mean(m['NAVz'])])
    BM2 = np.array([st.mean(m['BM2x']), st.mean(m['BM2y']), st.mean(m['BM2z'])])
    BM5 = np.array([st.mean(m['BM5x']), st.mean(m['BM5y']), st.mean(m['BM5z'])])
    
    Midfoot_anat, Omidfoota = Midfoot_cs.anat(NAV, BM2, BM5, footside)
    Midfoot_tech, Omidfoott = Midfoot_cs.tech(NAV,BM2,BM5)
    
else:
    sys.exit('ERROR: marker name not found in static file [Midfoot]')    
    
# Forefoot
if 'BM1' in labels and 'BM2' in labels and 'BM5' in labels and 'HM1' in labels and 'HM2' in labels and 'HM5' in labels:
    BM1 = np.array([st.mean(m['BM1x']), st.mean(m['BM1y']), st.mean(m['BM1z'])])
    BM2 = np.array([st.mean(m['BM2x']), st.mean(m['BM2y']), st.mean(m['BM2z'])])
    BM5 = np.array([st.mean(m['BM5x']), st.mean(m['BM5y']), st.mean(m['BM5z'])])
    HM1 = np.array([st.mean(m['HM1x']), st.mean(m['HM1y']), st.mean(m['HM1z'])])
    HM2 = np.array([st.mean(m['HM2x']), st.mean(m['HM2y']), st.mean(m['HM2z'])])
    HM5 = np.array([st.mean(m['HM5x']), st.mean(m['HM5y']), st.mean(m['HM5z'])])
    
    Forefoot_anat, Oforefoota = Forefoot_cs.anat(BM1, BM2, BM5, HM1, HM2, HM5, footside)
    Forefoot_tech, Oforefoott = Forefoot_cs.tech(BM1, BM5, HM2)
    
    if forefootanalysis==2:
        Forefootmed_anat, Oforefootmed = Forefoot_cs.med(BM1, BM2, HM1, HM2, footside)
        Forefootlat_anat, Oforefootlat = Forefoot_cs.lat(BM2, BM5, HM2, HM5, footside)
    elif forefootanalysis==1:
        Forefootlat_anat = nanmatrix
        Forefootmed_anat = nanmatrix
        Oforefootlat = nanmatrix
        Oforefootmed = nanmatrix
else:
    sys.exit('ERROR: marker name not found in static file [Forefoot]')
    
# Hallux
if 'HLX' in labels and 'HM1' in labels and 'HM2' in labels and 'HM5' in labels:
    HLX = np.array([st.mean(m['HLXx']), st.mean(m['HLXy']), st.mean(m['HLXz'])])
    HM1 = np.array([st.mean(m['HM1x']), st.mean(m['HM1y']), st.mean(m['HM1z'])])
    HM2 = np.array([st.mean(m['HM2x']), st.mean(m['HM2y']), st.mean(m['HM2z'])])
    HM5 = np.array([st.mean(m['HM5x']), st.mean(m['HM5y']), st.mean(m['HM5z'])])
    
    if forefootanalysis==2:
        HalluxFFm_anat, OhalluxFFm = Hallux_cs.anat(HLX, HM1, HM2, HM5, Forefootmed_anat)
        Hallux_anat, Ohalluxa = Hallux_cs.anat(HLX, HM1, HM2, HM5, Forefoot_anat) 
    elif forefootanalysis==1:
        HalluxFFm_anat = nanmatrix
        OhalluxFFm = nanmatrix
        Hallux_anat, Ohalluxa = Hallux_cs.anat(HLX, HM1, HM2, HM5, Forefoot_anat)
else:
    sys.exit('ERROR: marker name not found in static file [Hallux/forefoot]')
    
# Foot
if 'CALD' in labels and 'HM1' in labels and 'HM2' in labels and 'HM5' in labels:
    CALD = np.array([st.mean(m['CALDx']), st.mean(m['CALDy']), st.mean(m['CALDz'])])
    HM1 = np.array([st.mean(m['HM1x']), st.mean(m['HM1y']), st.mean(m['HM1z'])])
    HM2 = np.array([st.mean(m['HM2x']), st.mean(m['HM2y']), st.mean(m['HM2z'])])
    HM5 = np.array([st.mean(m['HM5x']), st.mean(m['HM5y']), st.mean(m['HM5z'])])
    
    Foot_anat, Ofoot = Foot_cs.anat(CALD, HM1, HM2, HM5, footside)

else:
    sys.exit('ERROR: marker name not found in static file [Foot segment]')

#### STATIC JOINT ANGLES #######
static_angles = dict()   
# Hindfoot - Shank angle
hfsk = np.dot(np.linalg.inv(Shank_anat),Hindfoot_anat)
[a,b,c] = RDecompose(hfsk)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b
    
static_angles['hfsk'] = np.array([[a,b,c]])

# Forefoot - Hindfoot angle
ffhf = np.dot(np.linalg.inv(Hindfoot_anat),Forefoot_anat)
[a,b,c] = RDecompose(ffhf)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b
    
static_angles['ffhf'] = np.array([[a,b,c]])

# Forefoot - Shank angle
ffsk = np.dot(np.linalg.inv(Shank_anat),Forefoot_anat)
[a,b,c] = RDecompose(ffsk)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b

static_angles['ffsk'] = np.array([[a,b,c]])

# Midfoot - Hindfoot angle
mfhf = np.dot(np.linalg.inv(Hindfoot_anat),Midfoot_anat)
[a,b,c] = RDecompose(mfhf)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b

static_angles['mfhf'] = np.array([[a,b,c]])

# Forefoot - Midfoot angle
ffmf = np.dot(np.linalg.inv(Midfoot_anat),Forefoot_anat)
[a,b,c] = RDecompose(ffmf)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b

static_angles['ffmf'] = np.array([[a,b,c]])

# Forefootmed - Midfoot angle
ffMmf = np.dot(np.linalg.inv(Midfoot_anat),Forefootmed_anat)
[a,b,c] = RDecompose(ffMmf)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b

static_angles['ffMmf'] = np.array([[a,b,c]])

# Forefootlat - Midfoot angle
ffLmf = np.dot(np.linalg.inv(Midfoot_anat),Forefootlat_anat)
[a,b,c] = RDecompose(ffLmf)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b

static_angles['ffLmf'] = np.array([[a,b,c]])

# Hallux - Forefoot angle
hxff = np.dot(np.linalg.inv(Forefoot_anat),Hallux_anat)
[a,b,c] = RDecompose(hxff)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b

static_angles['hxff'] = np.array([[a,b,c]])

# Hallux - Forefootmed angle
hxffM = np.dot(np.linalg.inv(Forefootmed_anat),HalluxFFm_anat)
[a,b,c] = RDecompose(hxffM)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b

static_angles['hxffM'] = np.array([[a,b,c]])

# Foot - Shank angle
fosk = np.dot(np.linalg.inv(Shank_anat),Foot_anat)
[a,b,c] = RDecompose(fosk)
angle_temp = np.array([[a,b,c]])
if footside==2: # left foot
    a = -a
    b = -b

static_angles['fosk'] = np.array([[a,b,c]])

# Rotation matrices for dynamic trials
RT_to_RA = dict()

RT_to_RA['Shank'] = np.linalg.solve(Shank_tech,Shank_anat)
RT_to_RA['Hindfoot'] = np.linalg.solve(Hindfoot_tech,Hindfoot_anat)
RT_to_RA['Midfoot'] = np.linalg.solve(Midfoot_tech,Midfoot_anat)
RT_to_RA['Forefoot'] = np.linalg.solve(Forefoot_tech,Forefoot_anat)
RT_to_RA['Forefootmed'] = np.linalg.solve(Forefootmed_anat,Forefootmed_anat)
RT_to_RA['Forefootlat'] = np.linalg.solve(Forefootlat_anat,Forefootlat_anat)
RT_to_RA['Hallux'] = np.linalg.solve(Hallux_anat,Hallux_anat)
RT_to_RA['Hallux'] = np.linalg.solve(HalluxFFm_anat,HalluxFFm_anat)
RT_to_RA['Foot'] = np.linalg.solve(Foot_anat,Foot_anat)

#### PA_MLA_static functie 
midSTPT = np.divide(PT+ST,2)
CAm = np.divide(CALP+midSTPT,2)
CAp = projectonplane(CAm, np.array([0,0,1]), np.array([0,0,0])) # projection on ground, lab coord syst = x=ML, y=Anterior, z=vertical
    
V1 = np.subtract(NAV,CAp)

# Second vector is between TN to the projection on the ground of MT1

HM1p = projectonplane(HM1, np.array([0,0,1]), np.array([0,0,0]))
V2 = np.subtract(NAV, HM1p)

# 3D angle between two vectors
temp = np.divide(np.dot(V1,V2), np.multiply(np.linalg.norm(V1), np.linalg.norm(V2)))
MLA = math.degrees(math.acos(temp))

# Express virtual marker position in relevent segment local reference frame 

#### LC_xform functie zit in de PA_MLA_static functie ####
# CAplocal
TT = np.subtract(CAp,Ohindfoott)

RR = Hindfoot_tech
temp = np.linalg.solve(RR, TT)
XD = np.array(temp)
CAplocal = XD

# HM1plocal
TT = np.subtract(HM1p,Oforefoott)

RR = Forefoot_tech
temp = np.linalg.solve(RR, TT)
XD = np.array(temp)
HM1plocal = XD

# Calculate thr transverse tarsal arch (TTA) for the Amsterdam Foot Model
TTA_var = TTA(BM1, BM2, BM5)

# Combining all static info in dicts for dynamic trial
localmarkers = dict()
staticinfo = dict()

localmarkers['CAplocal'] = CAplocal
localmarkers['HM1plocal'] = HM1plocal

static_angles['MLA'] = MLA
static_angles['TTA'] = TTA_var

staticinfo['RT_to_RA'] = RT_to_RA
staticinfo['localmarkers'] = localmarkers

#### AFM DYNAMIC FUNCTION ####
# Define gaitcycle
gc = np.subtract(HS_to_HS, startframe)

#### Initialize 3D matrices
if 'framerate_dyn' in globals():
    shank_tech = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    shank_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    hindfoot_tech = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    hindfoot_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    ohindfoott = np.zeros((len(m_dyn['CALDx']),1,3), dtype=object)
    midfoot_tech = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    midfoot_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    omidfoott = np.zeros((len(m_dyn['CALDx']),1,3), dtype=object)
    forefoot_tech = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    forefoot_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    oforefoott = np.zeros((len(m_dyn['CALDx']),1,3), dtype=object)
    forefootmed_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    forefootlat_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    hallux_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    halluxffm_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    foot_anat = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    hfsk = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    a_hfsk = list()
    b_hfsk = list()
    c_hfsk = list()
    a_ffhf = list()
    b_ffhf = list()
    c_ffhf = list()
    a_ffsk = list()
    b_ffsk = list()
    c_ffsk = list()
    a_mfhf = list()
    b_mfhf = list()
    c_mfhf = list()
    a_ffmf = list()
    b_ffmf = list()
    c_ffmf = list()
    a_ffMmf = list()
    b_ffMmf = list()
    c_ffMmf = list()
    a_ffLmf = list()
    b_ffLmf = list()
    c_ffLmf = list()
    a_hxff = list()
    b_hxff = list()
    c_hxff = list()
    a_hxffm = list()
    b_hxffm = list()
    c_hxffm = list()
    a_fosk = list()
    b_fosk = list()
    c_fosk = list()
    mla = list()
    tta = list()
    hfsk = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    ffhf = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    ffsk = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    mfhf = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    ffmf = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    ffMmf = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    ffLmf = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    hxff = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    hxffm = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
    fosk = np.zeros((len(m_dyn['CALDx']),3,3), dtype=object)
####

#### Dynamic coordinate system #### 
for i in range(0,len(m_dyn['CALDx'])):
    # Shank
    if 'TT' in labels_dyn and 'ASHN' in labels_dyn and 'LSHN' in labels_dyn:
        ASHN = np.array([m_dyn['ASHNx'][i], m_dyn['ASHNy'][i], m_dyn['ASHNz'][i]])
        TT = np.array([m_dyn['TTx'][i], m_dyn['TTy'][i], m_dyn['TTz'][i]])
        LSHN = np.array([m_dyn['LSHNx'][i], m_dyn['LSHNy'][i], m_dyn['LSHNz'][i]])
        
        shank_tech[i] =  Shank_cs.tech(ASHN, LSHN, TT)[0]
        shank_anat[i] = np.dot(shank_tech[i], RT_to_RA['Shank'])

    else:
        sys.exit('ERROR: marker name not found in dynamic file [Shank]')
        
    # Hindfoot
    if 'CALD' in labels_dyn and 'PT' in labels_dyn and 'ST' in labels_dyn:
        CALD = np.array([m_dyn['CALDx'][i], m_dyn['CALDy'][i], m_dyn['CALDz'][i]])
        PT = np.array([m_dyn['PTx'][i], m_dyn['PTy'][i], m_dyn['PTz'][i]])
        ST = np.array([m_dyn['STx'][i], m_dyn['STy'][i], m_dyn['STz'][i]])
        
        hindfoot_tech[i], ohindfoott[i] = Hindfoot_cs.tech(CALD, ST, PT)
        hindfoot_anat[i] = np.dot(hindfoot_tech[i], RT_to_RA['Hindfoot'])

    else:
        sys.exit('ERROR: marker name not found in dynamic file [Hindfoot]')
    
    # Midfoot
    if 'NAV' in labels_dyn and 'BM2' in labels_dyn and 'BM5' in labels_dyn:
       NAV = np.array([m_dyn['NAVx'][i], m_dyn['NAVy'][i], m_dyn['NAVz'][i]])
       BM2 = np.array([m_dyn['BM2x'][i], m_dyn['BM2y'][i], m_dyn['BM2z'][i]])
       BM5 = np.array([m_dyn['BM5x'][i], m_dyn['BM5y'][i], m_dyn['BM5z'][i]])
       
       midfoot_tech[i], omidfoott[i] = Midfoot_cs.tech(NAV, BM5, BM2)
       midfoot_anat[i] = np.dot(midfoot_tech[i], RT_to_RA['Midfoot'])

    else:
       sys.exit('ERROR: marker name not found in dynamic file [Midfoot]')
    
    # Forefoot
    if 'BM1' in labels and 'BM2' in labels and 'BM5' in labels and 'HM1' in labels and 'HM2' in labels and 'HM5' in labels:
        BM1 = np.array([m_dyn['BM1x'][i], m_dyn['BM1y'][i], m_dyn['BM1z'][i]])
        BM2 = np.array([m_dyn['BM2x'][i], m_dyn['BM2y'][i], m_dyn['BM2z'][i]])
        BM5 = np.array([m_dyn['BM5x'][i], m_dyn['BM5y'][i], m_dyn['BM5z'][i]])
        HM1 = np.array([m_dyn['HM1x'][i], m_dyn['HM1y'][i], m_dyn['HM1z'][i]])
        HM2 = np.array([m_dyn['HM2x'][i], m_dyn['HM2y'][i], m_dyn['HM2z'][i]])
        HM5 = np.array([m_dyn['HM5x'][i], m_dyn['HM5y'][i], m_dyn['HM5z'][i]])
        
        forefoot_tech[i], oforefoott[i] = Forefoot_cs.tech(BM1,BM5,HM2)
        forefoot_anat[i] = np.dot(forefoot_tech[i], RT_to_RA['Forefoot'])
        
        # Forefoot lat/med
        if forefootanalysis==2:
            forefootmed_anat[i] = Forefoot_cs.med(BM1, BM2, HM1, HM2, footside)[0]
            forefootlat_anat[i] = Forefoot_cs.lat(BM2, BM5, HM2, HM5, footside)[0]
        else:
            forefootmed_anat[i] = nanmatrix
            forefootlat_anat[i] = nanmatrix
    else:
       sys.exit('ERROR: marker name not found in dynamic file [Forefoot]')
            
    # Hallux
    if 'HLX' in labels and 'HM1' in labels and 'HM2' in labels and 'HM5' in labels:
        HLX = np.array([m_dyn['HLXx'][i], m_dyn['HLXy'][i], m_dyn['HLXz'][i]])
        HM1 = np.array([m_dyn['HM1x'][i], m_dyn['HM1y'][i], m_dyn['HM1z'][i]])
        HM2 = np.array([m_dyn['HM2x'][i], m_dyn['HM2y'][i], m_dyn['HM2z'][i]])
        HM5 = np.array([m_dyn['HM5x'][i], m_dyn['HM5y'][i], m_dyn['HM5z'][i]])
        
        if forefootanalysis==2:
            hallux_anat[i] = Hallux_cs.anat(HLX, HM1, HM2, HM5, forefoot_anat[i])[0]
            halluxffm_anat[i] = Hallux_cs.anat(HLX, HM1, HM2, HM5, forefootmed_anat[i])[0]
        elif forefootanalysis==1:
            hallux_anat[i] = Hallux_cs.anat(HLX, HM1, HM2, HM5, forefoot_anat[i])[0]
            halluxffm_anat[i] = nanmatrix
    else:
       sys.exit('ERROR: marker name not found in dynamic file [Hallux/Forefoot]')
        
    # Foot
    if 'CALD' in labels and 'HM1' in labels and 'HM2' in labels and 'HM5' in labels:
        CALD = np.array([m_dyn['CALDx'][i], m_dyn['CALDy'][i], m_dyn['CALDz'][i]])
        HM1 = np.array([m_dyn['HM1x'][i], m_dyn['HM1y'][i], m_dyn['HM1z'][i]])
        HM2 = np.array([m_dyn['HM2x'][i], m_dyn['HM2y'][i], m_dyn['HM2z'][i]])
        HM5 = np.array([m_dyn['HM5x'][i], m_dyn['HM5y'][i], m_dyn['HM5z'][i]])
        
        foot_anat[i] = Foot_cs.anat(CALD, HM1, HM2, HM5, footside)[0]
        
    else:
       sys.exit('ERROR: marker name not found in dynamic file [Foot segment]')
       
    #### DYNAMIC JOINT ANGLES #### 
    # Hindfoot - Shank angle
    hfsk[i] = np.dot(np.linalg.inv(shank_anat.astype(float)[i]),hindfoot_anat.astype(float)[i])
    [a,b,c] = RDecompose(hfsk[i])
    
    a_hfsk.append(a)
    b_hfsk.append(b)
    c_hfsk.append(c)
    
    # Forefoot - Hindfoot angle
    ffhf[i] = np.dot(np.linalg.inv(hindfoot_anat.astype(float)[i]),forefoot_anat.astype(float)[i])
    [a,b,c] = RDecompose(ffhf[i])
    
    a_ffhf.append(a)
    b_ffhf.append(b)
    c_ffhf.append(c)
    
    # Forefoot - Shank angle
    ffsk[i] = np.dot(np.linalg.inv(shank_anat.astype(float)[i]),forefoot_anat.astype(float)[i])
    [a,b,c] = RDecompose(ffsk[i])
    
    a_ffsk.append(a)
    b_ffsk.append(b)
    c_ffsk.append(c)
    
    # Midfoot - Hindfoot angle 
    mfhf[i] = np.dot(np.linalg.inv(hindfoot_anat.astype(float)[i]),midfoot_anat.astype(float)[i])
    [a,b,c] = RDecompose(mfhf[i])
    
    a_mfhf.append(a)
    b_mfhf.append(b)
    c_mfhf.append(c)

    # Forefoot - midfoot angle
    ffmf[i] = np.dot(np.linalg.inv(midfoot_anat.astype(float)[i]),forefoot_anat.astype(float)[i])
    [a,b,c] = RDecompose(ffmf[i])
    
    a_ffmf.append(a)
    b_ffmf.append(b)
    c_ffmf.append(c)

    # Forefootmed - midfoot angle
    ffMmf[i] = np.dot(np.linalg.inv(midfoot_anat.astype(float)[i]),forefootmed_anat.astype(float)[i])
    [a,b,c] = RDecompose(ffMmf[i])
    
    a_ffMmf.append(a)
    b_ffMmf.append(b)
    c_ffMmf.append(c)
    
    # Forefootlat - midfoot angle
    ffLmf[i] = np.dot(np.linalg.inv(midfoot_anat.astype(float)[i]),forefootlat_anat.astype(float)[i])
    [a,b,c] = RDecompose(ffLmf[i])
    
    a_ffLmf.append(a)
    b_ffLmf.append(b)
    c_ffLmf.append(c)

    # Hallux forefoot angle
    hxff[i] = np.dot(np.linalg.inv(forefoot_anat.astype(float)[i]),hallux_anat.astype(float)[i])
    [a,b,c] = RDecompose(hxff[i])
    
    a_hxff.append(a)
    b_hxff.append(b)
    c_hxff.append(c)

    # Hallux - forefootmed angle
    hxffm[i] = np.dot(np.linalg.inv(forefootmed_anat.astype(float)[i]),halluxffm_anat.astype(float)[i])
    [a,b,c] = RDecompose(hxffm[i])
    
    a_hxffm.append(a)
    b_hxffm.append(b)
    c_hxffm.append(c)    
    
    # Foot - shank angle
    fosk[i] = np.dot(np.linalg.inv(shank_anat.astype(float)[i]),foot_anat.astype(float)[i])
    [a,b,c] = RDecompose(fosk[i])
    
    a_fosk.append(a)
    b_fosk.append(b)
    c_fosk.append(c)

    # Planar angles
    #### PA_MLA_dynamic functie
    caplocal = staticinfo['localmarkers']['CAplocal']
    hm1plocal = staticinfo['localmarkers']['HM1plocal']

    cap = np.dot(hindfoot_tech[i] , caplocal) + ohindfoott[i]
    hm1p = np.dot(forefoot_tech[i] , hm1plocal) + oforefoott[i]
    
    v1 = NAV - cap
    v2 = NAV - hm1p
    
    u = np.dot(v1.flatten(),v2.flatten())
    l = np.dot(np.linalg.norm(v1.flatten()), np.linalg.norm(v2.flatten()))
    ul_temp = np.divide(u,l)
    ul = math.degrees(math.acos(ul_temp))
    
    mla.append(ul)
    
    tta.append(TTA(BM1, BM2, BM5))

gaitcycle = dict()
# Slice angle data for gaitcycle only 
tgc = t[gc[0]:gc[1]]
if footside==1: # Right foot
    gaitcycle['HFSKx'] = a_hfsk[gc[0]:gc[1]]
    gaitcycle['HFSKy'] = b_hfsk[gc[0]:gc[1]]
    gaitcycle['HFSKz'] = c_hfsk[gc[0]:gc[1]]
    gaitcycle['FFHFx'] = a_ffhf[gc[0]:gc[1]]
    gaitcycle['FFHFy'] = b_ffhf[gc[0]:gc[1]]
    gaitcycle['FFHFz'] = c_ffhf[gc[0]:gc[1]]
    gaitcycle['FFSKx'] = a_ffsk[gc[0]:gc[1]]
    gaitcycle['FFSKy'] = b_ffsk[gc[0]:gc[1]]
    gaitcycle['FFSKz'] = c_ffsk[gc[0]:gc[1]]
    gaitcycle['MFHFx'] = a_mfhf[gc[0]:gc[1]]
    gaitcycle['MFHFy'] = b_mfhf[gc[0]:gc[1]]
    gaitcycle['MFHFz'] = c_mfhf[gc[0]:gc[1]]
    gaitcycle['FFMFx'] = a_ffmf[gc[0]:gc[1]]
    gaitcycle['FFMFy'] = b_ffmf[gc[0]:gc[1]]
    gaitcycle['FFMFz'] = c_ffmf[gc[0]:gc[1]]
    gaitcycle['FFmMFx'] = a_ffMmf[gc[0]:gc[1]]
    gaitcycle['FFmMFy'] = b_ffMmf[gc[0]:gc[1]]
    gaitcycle['FFmMFz'] = c_ffMmf[gc[0]:gc[1]]
    gaitcycle['FFlMFx'] = a_ffLmf[gc[0]:gc[1]]
    gaitcycle['FFlMFy'] = b_ffLmf[gc[0]:gc[1]]
    gaitcycle['FFlMFz'] = c_ffLmf[gc[0]:gc[1]]
    gaitcycle['HXFFx'] = a_hxff[gc[0]:gc[1]]
    gaitcycle['HXFFy'] = b_hxff[gc[0]:gc[1]]
    gaitcycle['HXFFz'] = c_hxff[gc[0]:gc[1]]
    gaitcycle['HXFFMx'] = a_hxffm[gc[0]:gc[1]]
    gaitcycle['HXFFMy'] = b_hxffm[gc[0]:gc[1]]
    gaitcycle['HXFFMz'] = c_hxffm[gc[0]:gc[1]]
    gaitcycle['FOSKx'] = a_fosk[gc[0]:gc[1]]
    gaitcycle['FOSKy'] = b_fosk[gc[0]:gc[1]]
    gaitcycle['FOSKz'] = c_fosk[gc[0]:gc[1]]
    gaitcycle['MLA'] = mla[gc[0]:gc[1]]
    gaitcycle['TTA'] = tta[gc[0]:gc[1]]
    
if footside==2:
    gaitcycle['HFSKx'] = -(np.asarray(a_hfsk[gc[0]:gc[1]]))
    gaitcycle['HFSKy'] = -(np.asarray(b_hfsk[gc[0]:gc[1]]))
    gaitcycle['HFSKz'] = (np.asarray(c_hfsk[gc[0]:gc[1]]))
    gaitcycle['FFHFx'] = -(np.asarray(a_ffhf[gc[0]:gc[1]]))
    gaitcycle['FFHFy'] = -(np.asarray(b_ffhf[gc[0]:gc[1]]))
    gaitcycle['FFHFz'] = (np.asarray(c_ffhf[gc[0]:gc[1]]))
    gaitcycle['FFSKx'] = -(np.asarray(a_ffsk[gc[0]:gc[1]]))
    gaitcycle['FFSKy'] = -(np.asarray(b_ffsk[gc[0]:gc[1]]))
    gaitcycle['FFSKz'] = (np.asarray(c_ffsk[gc[0]:gc[1]]))
    gaitcycle['MFHFx'] = -(np.asarray(a_mfhf[gc[0]:gc[1]]))
    gaitcycle['MFHFy'] = -(np.asarray(b_mfhf[gc[0]:gc[1]]))
    gaitcycle['MFHFz'] = (np.asarray(c_mfhf[gc[0]:gc[1]]))
    gaitcycle['FFMFx'] = -(np.asarray(a_ffmf[gc[0]:gc[1]]))
    gaitcycle['FFMFy'] = -(np.asarray(b_ffmf[gc[0]:gc[1]]))
    gaitcycle['FFMFz'] = (np.asarray(c_ffmf[gc[0]:gc[1]]))
    gaitcycle['FFmMFx'] = -(np.asarray(a_ffMmf[gc[0]:gc[1]]))
    gaitcycle['FFmMFy'] = -(np.asarray(b_ffMmf[gc[0]:gc[1]]))
    gaitcycle['FFmMFz'] = (np.asarray(c_ffMmf[gc[0]:gc[1]]))
    gaitcycle['FFlMFx'] = -(np.asarray(a_ffLmf[gc[0]:gc[1]]))
    gaitcycle['FFlMFy'] = -(np.asarray(b_ffLmf[gc[0]:gc[1]]))
    gaitcycle['FFlMFz'] = (np.asarray(c_ffLmf[gc[0]:gc[1]]))
    gaitcycle['HXFFx'] = -(np.asarray(a_hxff[gc[0]:gc[1]]))
    gaitcycle['HXFFy'] = -(np.asarray(b_hxff[gc[0]:gc[1]]))
    gaitcycle['HXFFz'] = (np.asarray(c_hxff[gc[0]:gc[1]]))
    gaitcycle['HXFFMx'] = -(np.asarray(a_hxffm[gc[0]:gc[1]]))
    gaitcycle['HXFFMy'] = -(np.asarray(b_hxffm[gc[0]:gc[1]]))
    gaitcycle['HXFFMz'] = (np.asarray(c_hxffm[gc[0]:gc[1]]))
    gaitcycle['FOSKx'] = -(np.asarray(a_fosk[gc[0]:gc[1]]))
    gaitcycle['FOSKy'] = -(np.asarray(b_fosk[gc[0]:gc[1]]))
    gaitcycle['FOSKz'] = (np.asarray(c_fosk[gc[0]:gc[1]]))
    gaitcycle['MLA'] = mla[gc[0]:gc[1]]
    gaitcycle['TTA'] = tta[gc[0]:gc[1]]

# Normalize gaitcycle to 101 datapoints
gaitcycle_normalized = dict()
if normalize==1:
    # markerlist = list(gaitcycle.keys())
    for key in gaitcycle:
        gaitcycle_normalized[key] = norm_to_101(gaitcycle[key])

# Calculate ROM 
gaitcycle_ROM = dict()
for key in gaitcycle:
    gaitcycle_ROM[key] = [min(gaitcycle[key]), max(gaitcycle[key])]
    
# Combine data 
Dynamic_data = dict()
Dynamic_data['gaitcycle'] = gaitcycle
Dynamic_data['gaitcycle_normalized'] = gaitcycle_normalized
Dynamic_data['gaitcycle_ROM'] = gaitcycle_ROM

#### PLOTTING JONT ANGLES ####
if normalize==1: 
    yy = Dynamic_data['gaitcycle_normalized']
else:
    yy = Dynamic_data['gaitcycle']

footlabels = ['FOSK','HFSK','FFHF','HXFFM']
midfootlabels = ['MFHF','FFMF','FFmMF','FFlMF']
archlabels = ['MLA', 'TTA']
planes = ['z','x','y']

if forefootanalysis==1:
    footlabels[3]='HXFF'
    midfootlabels = ['MFHF','FFMF']

# Subplots: forefoot
f1, ax = plt.subplots(4,3, figsize=(15,10))
f1.subplots_adjust(hspace=0.7, wspace=0.3)
f1.delaxes(ax[3,1])
f1.suptitle('Amsterdam Foot Model - Right Foot', fontsize = 20)

for row, seg in enumerate(footlabels):
    for col, plane in enumerate(planes):
        numofpoints = list(range(1,len(yy[seg + plane])+1))
        temp_mean = np.mean(yy[seg + plane])
        temp_std = np.std(yy[seg + plane])
        temp_ub = temp_mean + temp_std
        temp_lb = temp_mean - temp_std
        if row==0 and col==0:
            ax[row,col].fill_between(numofpoints,temp_ub,temp_lb, facecolor = 'lightblue', label = '± 1 std range')
            ax[row,col].plot(yy[seg + plane], label = 'Dynamic trial')
            ax[row,col].set_xlabel('Gait cycle [%]')
            ax[row,col].set_xlim(numofpoints[0], numofpoints[-1])
        else:
            ax[row,col].fill_between(numofpoints,temp_ub,temp_lb, facecolor = 'lightblue')
            ax[row,col].plot(yy[seg + plane])
            ax[row,col].set_xlabel('Gait cycle [%]')
            ax[row,col].set_xlim(numofpoints[0], numofpoints[-1])

# Y labels
for row in range(0,4):
    ax[row,0].set_ylabel('Pfl(-) / Dfl(+)', fontsize=9)
for row in [0,2]:
    ax[row,1].set_ylabel('Ev(-) / Inv(+)', fontsize=9)
for row in [1,2]:
    if row==1:
        col = 1
        ax[row,col].set_ylabel('Valg(-) / Var(+)', fontsize=9)
    elif row==2:
        col = 2
        ax[row,col].set_ylabel('Valg(-) / Var(+)', fontsize=9)
for row in [0,1]:
        ax[row,2].set_ylabel('Exo(-) / Endo(+)', fontsize=9)
ax[2,2].set_ylabel('Abd(-) / Add(+)', fontsize=9)

# Row labels
for ax, row in zip(ax[:,0], footlabels):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - 10, 0), xycoords=ax.yaxis.label, textcoords='offset points', size='large', ha='right', va='center')

f1.tight_layout           
f1.legend(loc = (0.45,0.18))
f1.align_ylabels()

# Subplots: midfoot
if forefootanalysis==1:
    numrows = 2
elif forefootanalysis==2:
    numrows = 4
f2, ax2 = plt.subplots(numrows,3, figsize=(15,10))
f2.subplots_adjust(hspace=0.7, wspace=0.4)
f2.suptitle('Amsterdam Foot Model - Right Midfoot', fontsize = 20)

for row, seg in enumerate(midfootlabels):
    for col, plane in enumerate(planes):
        numofpoints = list(range(1,len(yy[seg + plane])+1))
        temp_mean = np.mean(yy[seg + plane])
        temp_std = np.std(yy[seg + plane])
        temp_ub = temp_mean + temp_std
        temp_lb = temp_mean - temp_std  
        if row==0 and col==0:
            ax2[row,col].fill_between(numofpoints,temp_ub,temp_lb, facecolor = 'lightblue', label = '± 1 std range')
            ax2[row,col].plot(yy[seg + plane], label = 'Dynamic trial')
            ax2[row,col].set_xlabel('Gait cycle [%]')
            ax2[row,col].set_xlim(numofpoints[0], numofpoints[-1])
        else:
            ax2[row,col].fill_between(numofpoints,temp_ub,temp_lb, facecolor = 'lightblue')
            ax2[row,col].plot(yy[seg + plane])
            ax2[row,col].set_xlabel('Gait cycle [%]')
            ax2[row,col].set_xlim(numofpoints[0], numofpoints[-1])

# Y labels
for row in range(0,numrows):
    ax2[row,0].set_ylabel('Pfl(-) / Dfl(+)', fontsize=9)
    ax2[row,1].set_ylabel('Ev(-) / Inv(+)', fontsize=9)
    ax2[row,2].set_ylabel('Abd(-) / Add(+)', fontsize=9)

# Row labels
for ax2, row in zip(ax2[:,0], midfootlabels):
    ax2.annotate(row, xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - 10, 0), xycoords=ax2.yaxis.label, textcoords='offset points', size='large', ha='right', va='center')

f2.tight_layout           
f2.align_ylabels()

# Subplots: arches
f3, (ax31, ax32) = plt.subplots(2, sharex=True, sharey=True, figsize=(10,5))

MLA_mean = np.mean(yy['MLA'])
TTA_mean = np.mean(yy['TTA'])
MLA_std = np.std(yy['MLA'])
TTA_std = np.std(yy['TTA'])

UB_mla = MLA_mean + MLA_std
LB_mla = MLA_mean - MLA_std
UB_TTA = TTA_mean + TTA_std
LB_TTA = TTA_mean - TTA_std

# ax31
ax31.plot(yy['MLA'])
ax31.set_ylabel('MLA', fontsize=12)
ax31.set_ylim(80,140)
ax31.set_xlim(0,len(yy['MLA'])-1)
ax31.fill_between(range(0,len(yy['MLA'])), UB_mla, LB_mla, facecolor='lightblue')

# ax32
ax32.plot(yy['TTA'])
ax32.set_ylabel('TTA', fontsize=12)
ax32.fill_between(range(0,len(yy['TTA'])), UB_TTA, LB_TTA, facecolor='lightblue')

f3.suptitle('Amsterdam Foot Model - Right Arches', fontsize = 20)
f3.supxlabel('Gait cycle [%]')
plt.show()

import savefig as save

print(type(f1))

save.jpeg(f1,'AFM_Kinematics','foot')
save.jpeg(f2,'AFM_Kinematics','midfoot')
save.jpeg(f3,'AFM_Kinematics','arches')

    