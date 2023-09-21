# Rdecompose
# Decomposition of rotation matrix

import math as m

def RDecompose(R):
    R = R.astype(float)
    alpha = m.asin(R[2,1])
    gamma = m.atan2(-R[0,1], R[1,1])
    beta = m.atan2(-R[2,0], R[2,2])
    
    alpha = m.degrees(alpha)
    gamma = m.degrees(gamma)
    beta = m.degrees(beta)
    
    return alpha,beta,gamma