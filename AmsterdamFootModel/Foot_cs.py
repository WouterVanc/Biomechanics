# Calculate local coordinate system for the foot according to Cappazo et al (1995)
# R = local axis
# O = origin

import numpy as np
from Projection import projectonplane

def anat(CALD, HM1, HM2, HM5, footside):
    # Origin
    O = np.transpose(np.array(CALD))
    
    a1 = HM1-O
    a2 = HM5-O
    if footside==1: # right foot
        ytemp = np.cross(a2,a1)
    elif footside==2: # left foot
        ytemp = np.cross(a1,a2)
    
    ynorm = ytemp / np.linalg.norm(ytemp)
    
    projHM2 = projectonplane(HM2, ynorm, HM5)
    
    # axis
    x_axis = projHM2-O
    z_axis = np.cross(x_axis,ytemp)
    y_axis = np.cross(z_axis,x_axis)
    
    # Convert to unit vectors
    X = x_axis / np.linalg.norm(x_axis)
    Y = y_axis / np.linalg.norm(y_axis)
    Z = z_axis / np.linalg.norm(z_axis)
    
    R = np.transpose(np.array([X,Y,Z]))
    
    return R,O
    