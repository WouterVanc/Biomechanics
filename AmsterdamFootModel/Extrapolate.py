# Extrapolate to same frequency


from scipy.interpolate import InterpolatedUnivariateSpline as ius
import numpy as np

def norm_to_fp(y):
    datapoints = list(range(1,len(y)+1))
    
    
    c = ius(datapoints,y,k=3)
    xq = np.linspace(1, len(y)+1, int(len(y)/10))
    
    c_interp = c(xq)
    
    return c_interp