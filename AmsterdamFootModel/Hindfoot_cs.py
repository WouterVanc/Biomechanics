# Anatomical/technical Hindfoot Coordinate System

import numpy as np

def anat(CALD,CALP,ST,PT):
    
    # Origin
    Ohindfoot = np.divide(CALD+CALP, 2)
    # Anterior Posterior axis
    x = np.divide(ST+PT, 2) - Ohindfoot
    # Medial Lateral axis
    z = np.cross(x,CALP-CALD)
    # Vertical axis
    y = np.cross(z,x)

    # Convert to unit vectors
    X = x / np.linalg.norm(x)
    Y = y / np.linalg.norm(y)
    Z = z / np.linalg.norm(z)
    
    # Final coordinate system
    Hindfoot = np.transpose(np.array([X,Y,Z]))
    
    return Hindfoot, Ohindfoot

def tech(CALD,ST,PT):
    
    # Origin
    Ohindfoot = CALD
    # Axis
    first_axis = ST-CALD
    third_axis_temp = PT-CALD
    
    second_axis = np.cross(third_axis_temp, first_axis)
    third_axis = np.cross(first_axis, second_axis)
    
    # Convert to unit vectors
    axis1 = first_axis / np.linalg.norm(first_axis)
    axis2 = second_axis / np.linalg.norm(second_axis)
    axis3 = third_axis / np.linalg.norm(third_axis)
    
    # Final coordinate system
    Hindfoot = np.transpose(np.array([axis1, axis2, axis3]))
    
    return Hindfoot, Ohindfoot
        

