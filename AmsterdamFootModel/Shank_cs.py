# Anatomical/technical Shank Coordinate System

import numpy as np
from Projection import projectonplane 

def anat(ASHN,LSHN,FH,TT,LM,MM,footside):
    # Define axes
    if footside==1: # right foot
        # Origin
        Oshank = np.divide(LM + MM,2)
        # Anterior posterior axis (x)
        x = np.cross(MM - LM, FH - LM)
        # Vertical axis (y)
        TT_proj = projectonplane(TT, np.divide(x, np.linalg.norm(x)), LM)
        y = TT_proj-Oshank
        # Medial lateral axis (z)
        z = np.cross(x,y)
        
    elif footside==2: # left foot
        # Origin
        Oshank = np.divide(LM + MM,2)
        # Anterior posterior axis (x)
        x = np.cross(MM - LM, FH - LM)
        # Vertical axis (y)
        TT_proj = projectonplane(TT, np.divide(x, np.linalg.norm(x)), LM)
        y = TT_proj-Oshank
        # Medial lateral axis (z)
        z = np.cross(x,y)
        
    # Convert to unit vectors
    X = x / np.linalg.norm(x)
    Y = y / np.linalg.norm(y)
    Z = z / np.linalg.norm(z)
    
    # Final coordinate system
    Shank = np.transpose(np.array([X,Y,Z]))
    
    return Shank, Oshank

def tech(ASHN,LSHN,TT):
    
    # Origin
    Oshank = TT
    # Axis
    first_axis = LSHN-TT
    third_axis_temp = ASHN-TT
    
    second_axis = np.cross(third_axis_temp, first_axis)
    third_axis = np.cross(first_axis, second_axis)
    
    # Convert to unit vectors
    axis1 = first_axis / np.linalg.norm(first_axis)
    axis2 = second_axis / np.linalg.norm(second_axis)
    axis3 = third_axis / np.linalg.norm(third_axis)
    
    # Final coordinate system
    Shank = np.transpose(np.array([axis1, axis2, axis3]))
    
    return Shank, Oshank


