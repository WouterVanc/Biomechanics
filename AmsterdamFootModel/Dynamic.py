# Dynamic AFM function

# Import packages
from ezc3d import c3d
from datafilter import filter_withnan
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
from normalize import norm_to_101

def Dynamic_AFM(path,staticinfo,HS_to_HS,startframe,footside,forefootanalysis):
    # load in data
    c_dyn = c3d(path)

    # Extract variables 
    point_data_dyn= c_dyn['data']['points']
    labels_dyn = list(c_dyn['parameters']['POINT']['LABELS']['value'])
    points_residuals_dyn = c_dyn['data']['meta_points']['residuals']
    analog_data_dyn = c_dyn['data']['analogs']
    framerate_dyn = int(c_dyn['header']['points']['frame_rate'])
    lastframe_dyn = int(c_dyn['header']['points']['last_frame'])

    # Create time list
    t = list(range(0,lastframe_dyn+1))
    t = [x/100 for x in t]

    # Low pass filter data
    point_data_dyn = filter_withnan(point_data_dyn, framerate_dyn)

    # Allocate data to marker names for dynamic trial
    m_dyn = dict()
    for ind, name in enumerate(labels_dyn):
        m_dyn[name + 'x'] = point_data_dyn[0,ind,:]
        m_dyn[name + 'y'] = point_data_dyn[1,ind,:]
        m_dyn[name + 'z'] = point_data_dyn[2,ind,:]

    #### AFM DYNAMIC FUNCTION ####
    # Define gaitcycle
    gc = np.subtract(HS_to_HS, startframe)
    nanmatrix = np.array([float('nan')]*9).reshape(3,3)

    #### Dynamic coordinate system #### 
    RT_to_RA = staticinfo['RT_to_RA']

    for i in range(0,len(m_dyn['CALDx'])):
        if i == 0:
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
        if 'BM1' in labels_dyn and 'BM2' in labels_dyn and 'BM5' in labels_dyn and 'HM1' in labels_dyn and 'HM2' in labels_dyn and 'HM5' in labels_dyn:
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
            elif forefootanalysis==1:
                forefootmed_anat[i] = nanmatrix
                forefootlat_anat[i] = nanmatrix

        else:
            sys.exit('ERROR: marker name not found in dynamic file [Forefoot]')
                
        # Hallux
        if 'HLX' in labels_dyn and 'HM1' in labels_dyn and 'HM2' in labels_dyn and 'HM5' in labels_dyn:
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
        if 'CALD' in labels_dyn and 'HM1' in labels_dyn and 'HM2' in labels_dyn and 'HM5' in labels_dyn:
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
    
    return Dynamic_data 