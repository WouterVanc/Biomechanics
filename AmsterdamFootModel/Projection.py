# Function: Projection of marker on plane 
# Normvector should be perpendicular on plane
# Point should be on plane

import numpy as np

def projectonplane(marker,normvector,point):
    
        normvectortrans = np.transpose(normvector)
        
        vector = np.subtract(marker , point)
        
        dotproduct = np.dot(vector, normvectortrans)
        
        temp = dotproduct * normvectortrans
        
        projected_marker = np.subtract(marker , temp)
        
        return projected_marker