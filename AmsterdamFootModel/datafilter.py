# Filter data
# Low pass filter data while ignoring NaN values

# Input: 3D array of time series data marker (dimension one, two and three represent x, y and z coord), framerate
# Ouput: Identical 3D array that has been low pass filtered with nans in the same place
import numpy as np
from ButterFilterLow import LowpassFilter
import sys

def nanfilter(data, framerate):
    size = np.shape(data)
    for dim in range(0,3):
        # Check for full unlabed frames
        for i in range(0,size[2]):
            col = data[dim,:,i]
            if np.isnan(col).any():
                togical = np.isnan(col)
                if all(togical):
                    sys.exit('C3D file contains frame(s) with no labeled markers')
        # Apply filter 
        for idx, row in enumerate(data[dim]):
            if np.isnan(row).any():
                temp = row
                logical = np.isnan(temp)
                if all(logical):
                    continue
                filter_no_nans = LowpassFilter(temp[logical==False], framerate)
                temp[logical==False] = filter_no_nans
                
                data[dim,idx,:] = temp
            else:
                data[dim,idx,:] = LowpassFilter(data[dim,idx,:], framerate)
                
    return data


    
