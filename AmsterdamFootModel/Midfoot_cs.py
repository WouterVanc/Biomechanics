# Anatomical/technical Midfoot Coordinate System

import numpy as np

def anat(NAV,BM2,BM5,footside):
    if footside==1: #right foot
        # Origin
        Omidfoot = np.divide(NAV+BM5,2)
        # Anterior-posterior 
        x = BM2-Omidfoot
        # Verical
        y = np.cross(BM5-NAV,BM2-NAV)
        # Medio-latero
        z = np.cross(x,y)
        
    elif footside==2: #left foot
        # Origin
        Omidfoot = np.divide(NAV+BM5,2)
        # Anterior-posterior 
        x = BM2-Omidfoot
        # Verical
        y = np.cross(NAV-BM5,BM2-NAV)
        # Medio-latero
        z = np.cross(x,y)
        
    # Convert to unit vectors
    X = x / np.linalg.norm(x)
    Y = y / np.linalg.norm(y)
    Z = z / np.linalg.norm(z)
    
    # Final coordinate system
    Hindfoot = np.transpose(np.array([X,Y,Z]))
    
    return Hindfoot, Omidfoot

def tech(NAV,BM2,BM5):
    
    # Origin
    Omidfoot = np.divide(NAV+BM5,2)
    # Axis
    x = BM2-Omidfoot
    y = np.cross(BM5-NAV,BM2-NAV)
    z = np.cross(x,y)
    
    # Convert to unit vectors
    X = x / np.linalg.norm(x)
    Y = y / np.linalg.norm(y)
    Z = z / np.linalg.norm(z)
    
    # Final coordinate system
    Midfoot = np.transpose(np.array([X,Y,Z]))
    
    return Midfoot, Omidfoot
    
    
    
    

    