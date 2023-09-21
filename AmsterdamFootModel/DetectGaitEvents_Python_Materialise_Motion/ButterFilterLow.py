# Wouter Van Caekenberghe, 16/03/2023

from scipy.signal import butter
import scipy

# Function to low pass filter data
# INPUT: array of values
# OUTPUT: filtered array
def LowpassFilter(data,SampleFreq):
    # Initialize filter variables
    order = 4
    Fs = SampleFreq # sample frequency, Hz = fps
    Fc = 6 # critical frequency, Hz = cut off
    Wn = (2*Fc)/Fs
    
    # Apply butterworth filter
    b,a = butter(order,Wn,'low', analog=False)
    y = scipy.signal.filtfilt(b,a,data, padlen=len(data)-1)
    
    return y 

