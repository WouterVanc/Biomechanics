# Static AFM function

# Import packages
from ezc3d import c3d
from datafilter import nanfilter
import statistics as st
import numpy as np
import math 
import sys

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
from Markernames import adjust_marker_names_L, adjust_marker_names_R

def Static_AFM(path,footside,forefootanalysis):
    # load in data
    c = c3d(path)   

    # Extract variables 
    point_data_orig = c['data']['points']
    labels = list(c['parameters']['POINT']['LABELS']['value'])
    framerate = int(c['header']['points']['frame_rate'])

    # Remove NaN values from data
    logical = ~np.isnan(point_data_orig[0]).any(axis=0)
    point_data = point_data_orig[:,:,logical]

    if footside==1:
        labels = adjust_marker_names_R(labels)
    if footside==2:
        labels = adjust_marker_names_L(labels)

    # Low pass filter data
    point_data = nanfilter(point_data, framerate)

    # Allocate data to marker names for static trial
    m = dict()
    for ind, name in enumerate(labels):
        m[name + 'x'] = point_data[0,ind,:]
        m[name + 'y'] = point_data[1,ind,:]
        m[name + 'z'] = point_data[2,ind,:]

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
    if footside==2: # left foot
        a = -a
        b = -b
        
    static_angles['hfsk'] = [a,b,c]

    # Forefoot - Hindfoot angle
    ffhf = np.dot(np.linalg.inv(Hindfoot_anat),Forefoot_anat)
    [a,b,c] = RDecompose(ffhf)
    if footside==2: # left foot
        a = -a
        b = -b
        
    static_angles['ffhf'] = [a,b,c]

    # Forefoot - Shank angle
    ffsk = np.dot(np.linalg.inv(Shank_anat),Forefoot_anat)
    [a,b,c] = RDecompose(ffsk)
    if footside==2: # left foot
        a = -a
        b = -b

    static_angles['ffsk'] = [a,b,c]

    # Midfoot - Hindfoot angle
    mfhf = np.dot(np.linalg.inv(Hindfoot_anat),Midfoot_anat)
    [a,b,c] = RDecompose(mfhf)
    if footside==2: # left foot
        a = -a
        b = -b

    static_angles['mfhf'] = [a,b,c]

    # Forefoot - Midfoot angle
    ffmf = np.dot(np.linalg.inv(Midfoot_anat),Forefoot_anat)
    [a,b,c] = RDecompose(ffmf)
    if footside==2: # left foot
        a = -a
        b = -b

    static_angles['ffmf'] = [a,b,c]

    # Forefootmed - Midfoot angle
    ffMmf = np.dot(np.linalg.inv(Midfoot_anat),Forefootmed_anat)
    [a,b,c] = RDecompose(ffMmf)
    if footside==2: # left foot
        a = -a
        b = -b

    static_angles['ffMmf'] = [a,b,c]

    # Forefootlat - Midfoot angle
    ffLmf = np.dot(np.linalg.inv(Midfoot_anat),Forefootlat_anat)
    [a,b,c] = RDecompose(ffLmf)
    if footside==2: # left foot
        a = -a
        b = -b

    static_angles['ffLmf'] = [a,b,c]

    # Hallux - Forefoot angle
    hxff = np.dot(np.linalg.inv(Forefoot_anat),Hallux_anat)
    [a,b,c] = RDecompose(hxff)
    if footside==2: # left foot
        a = -a
        b = -b

    static_angles['hxff'] = [a,b,c]

    # Hallux - Forefootmed angle
    hxffM = np.dot(np.linalg.inv(Forefootmed_anat),HalluxFFm_anat)
    [a,b,c] = RDecompose(hxffM)
    if footside==2: # left foot
        a = -a
        b = -b

    static_angles['hxffM'] = [a,b,c]

    # Foot - Shank angle
    fosk = np.dot(np.linalg.inv(Shank_anat),Foot_anat)
    [a,b,c] = RDecompose(fosk)
    if footside==2: # left foot
        a = -a
        b = -b

    static_angles['fosk'] = [a,b,c]

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

    #### PA_MLA_static function 
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

    return static_angles, staticinfo