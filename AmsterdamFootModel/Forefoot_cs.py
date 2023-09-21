# Anatomical/technical Forefoot Coordinate System

import numpy as np

def anat(BM1,BM2,BM5,HM1,HM2,HM5,footside):
    
    # Create virtual markers
    MH12=np.divide(HM1+HM2,2)
    MH25=np.divide(HM2+HM5,2)
    MB12=np.divide(BM1+BM2,2)
    MB25=np.divide(BM2+BM5,2)
    
    if footside==1: # right foot
        # Origin
        Oforefoot = np.divide(MB12+MB25, 2)
        # Anterior Posterior
        x = np.divide(MH25+MH12,2)-Oforefoot
        # Vertical
        y = np.cross(MH25-Oforefoot,MH12-Oforefoot)
        # Medio Latero
        z = np.cross(x,y)
        
    elif footside==2:
        # Origin
        Oforefoot = np.divide(MB12+MB25, 2)
        # Anterior Posterior
        x = np.divide(MH25+MH12,2)-Oforefoot
        # Vertical
        y = np.cross(Oforefoot-MH25,MH12-Oforefoot)
        # Medio Latero
        z = np.cross(x,y)
        
    # Convert to unit vectors
    X = x / np.linalg.norm(x)
    Y = y / np.linalg.norm(y)
    Z = z / np.linalg.norm(z)
    
    # Final coordinate system
    Forefoot = np.transpose(np.array([X,Y,Z]))
    
    return Forefoot, Oforefoot

def tech(BM1,BM5,HM2):
    
    # Origin
    Oforefoot = BM1
    # Axis
    first_axis = BM5-Oforefoot
    third_axis_temp = HM2-Oforefoot
    
    second_axis = np.cross(third_axis_temp, first_axis)
    third_axis = np.cross(first_axis, second_axis)
    
    # Convert to unit vectors
    axis1 = first_axis / np.linalg.norm(first_axis)
    axis2 = second_axis / np.linalg.norm(second_axis)
    axis3 = third_axis / np.linalg.norm(third_axis)
    
    # Final coordinate system
    Forefoot = np.transpose(np.array([axis1, axis2, axis3]))
    
    return Forefoot, Oforefoot
    
def med(BM1,BM2,HM1,HM2,footside):
    
    MH12=np.divide(HM1+HM2,2)
    MB12=np.divide(BM1+BM2,2)
    
    if footside==1: # right foot
        # Origin
        Oforefootmed = MB12
        # Anterior Posterior
        x = MH12-Oforefootmed
        # Verical
        y = np.cross(BM2-BM1,x)
        # Medio latero
        z = np.cross(x,y)
        
    elif footside==2: # left foot
        # Origin
        Oforefootmed = MB12
        # Anterior Posterior
        x = MH12-Oforefootmed
        # Verical
        y = np.cross(BM1-BM2,x)
        # Medio latero
        z = np.cross(x,y)
    
    # Convert to unit vectors
    X = x / np.linalg.norm(x)
    Y = y / np.linalg.norm(y)
    Z = z / np.linalg.norm(z)
    
    # Final coordinate system
    forefootmed = np.transpose(np.array([X,Y,Z]))
    
    return forefootmed, Oforefootmed

def lat(BM2,BM5,HM2,HM5,footside):
    
    MH25=np.divide(HM2+HM5,2)
    MB25=np.divide(BM2+BM5,2)
    
    if footside==1: # right foot
        # Origin
        Oforefootmed = MB25
        # Anterior Posterior
        x = MH25-MB25
        # Verical
        y = np.cross(BM5-BM2,x)
        # Medio latero
        z = np.cross(x,y)
        
    elif footside==2: # left foot
        # Origin
        Oforefootmed = MB25
        # Anterior Posterior
        x = MH25-MB25 
        # Verical
        y = np.cross(BM2-BM5,x)
        # Medio latero
        z = np.cross(x,y)
    
    # Convert to unit vectors
    X = x / np.linalg.norm(x)
    Y = y / np.linalg.norm(y)
    Z = z / np.linalg.norm(z)
    
    # Final coordinate system
    forefootmed = np.transpose(np.array([X,Y,Z]))
    
    return forefootmed, Oforefootmed
    
   
    
    
    
