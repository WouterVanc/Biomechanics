# Anatomical/technical Hallux Coordinate System

import numpy as np

def anat(HLX,HM1,HM2,HM5,forefoot):
    
    # MH12=np.divide(HM1+HM2,2);
    # MH25=np.divide(HM2+HM5,2);
    forefootfloat = forefoot.astype(float)
    
    # Origin
    Ohallux = HM1
    # Anterior posterior
    x = HLX-Ohallux
    # Vertical
    # y = np.cross(MH25-MH12,x) 
    y = np.cross(forefootfloat[:,2],x) # When directly using the ML axis of the forefoot
    # Medio Latero
    z = np.cross(x,y)
    
    # Convert to unit vectors
    X = x / np.linalg.norm(x)
    Y = y / np.linalg.norm(y)
    Z = z / np.linalg.norm(z)
    
    # Final coordinate system
    Hallux = np.transpose(np.array([X,Y,Z]))
    
    return Hallux, Ohallux
    
    
