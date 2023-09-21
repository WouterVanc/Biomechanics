# Calculates the transverse tarsal arch (TTA) for the Amsterdam Foot Model
import numpy as np
import math as m

def TTA(BM1, BM2, BM5):
    V1 = np.subtract(BM2, BM1)
    V2 = np.subtract(BM2, BM5)
    
    #3d angle between two vectors
    temp = np.divide(np.dot(V1,V2),np.multiply(np.linalg.norm(V1), np.linalg.norm(V2)))
    
    TTA = m.degrees(m.acos(temp))
    
    return TTA

