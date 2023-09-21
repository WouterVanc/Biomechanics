# Function to normalize to 101 datapoints

from scipy.interpolate import InterpolatedUnivariateSpline as ius
import numpy as np
import sys
def norm_to_101(y):
    if len(y)==0:
        sys.exit('ERROR: input data to normalize function is empty')
    datapoints = list(range(1,len(y)+1))
    
    c = ius(datapoints,y,k=3)
    xq = np.linspace(1, len(y)+1, 101)
    
    c_interp = c(xq)
    
    return c_interp