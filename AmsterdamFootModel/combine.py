# Find average and standev across trials per condition
import numpy as np
import pandas as pd

def avg_groupdata_right(groupdata, conditions, anglenames):
    avgdata = dict()
    for c in conditions:
        numoftrials = len(groupdata[c + '_right'])
        averages = dict()
        standev = dict()
        for name in anglenames:    
            tempi = np.zeros([1,len(groupdata[c + '_right'][0]['gaitcycle_normalized']['FOSKx'])])
            tempy = np.zeros([numoftrials,len(groupdata[c + '_right'][0]['gaitcycle_normalized']['FOSKx'])])
            for i in range(0,numoftrials):
                tempi = np.add(tempi, groupdata[c + '_right'][i]['gaitcycle_normalized'][name])
                tempy[i,:] = groupdata[c + '_right'][i]['gaitcycle_normalized'][name]
                tempy_std = list()
            for ii in range(0,tempy.shape[1]):
                tempy_std.append(np.std(tempy[:,ii]))
            tempi_avg = np.divide(tempi,numoftrials)
            tempi_avg = tempi_avg.flatten()
                
            averages[name] = tempi_avg
            standev[name] = tempy_std
        avgdata[c + '_avg'] = averages
        avgdata[c + '_std'] = standev
    
    for c in conditions:
        zero_data = np.zeros(shape=(1,len(groupdata[c + '_right'][0]['gaitcycle_ROM'])))
        ROM = pd.DataFrame(zero_data, columns=list(groupdata[c + '_right'][0]['gaitcycle_ROM'].keys()))
        ROMS = pd.DataFrame(zero_data, columns=list(groupdata[c + '_right'][0]['gaitcycle_ROM'].keys()))
        ROM_std = pd.DataFrame(zero_data, columns=list(groupdata[c + '_right'][0]['gaitcycle_ROM'].keys()))
        numoftrials = len(groupdata[c + '_right'])
        for i in range(0,numoftrials):
            df = pd.DataFrame([groupdata[c + '_right'][i]['gaitcycle_ROM']])
            ROM = ROM.add(df)
            ROMS = pd.concat([ROMS, df], ignore_index=True, sort=False)
        
        ROM = ROM.divide(numoftrials)
        ROMS = ROMS.iloc[1:]
        for (colnames, colvalues) in ROMS.items():
            ROM_std.iloc[0].at[colnames] = np.std(colvalues)
                                                                                                    
        avgdata[c + '_avg_ROM'] = ROM
        avgdata[c + '_std_ROM'] = ROM_std

    return avgdata 

def avg_groupdata_left(groupdata, conditions, anglenames):
    avgdata = dict()
    for c in conditions:
        numoftrials = len(groupdata[c + '_left'])
        averages = dict()
        standev = dict()
        for name in anglenames:
            tempi = np.zeros([1,len(groupdata[c + '_left'][0]['gaitcycle_normalized']['FOSKx'])])
            tempy = np.zeros([numoftrials,len(groupdata[c + '_left'][0]['gaitcycle_normalized']['FOSKx'])])
            for i in range(0,numoftrials):
                tempi = np.add(tempi, groupdata[c + '_left'][i]['gaitcycle_normalized'][name])
                tempy[i,:] = groupdata[c + '_left'][i]['gaitcycle_normalized'][name]
                tempy_std = list()
            for ii in range(0,tempy.shape[1]):
                tempy_std.append(np.std(tempy[:,ii]))
            tempi_avg = np.divide(tempi,numoftrials)
            tempi_avg = tempi_avg.flatten()
                
            averages[name] = tempi_avg
            standev[name] = tempy_std
        avgdata[c + '_avg'] = averages
        avgdata[c + '_std'] = standev

        for c in conditions:
            zero_data = np.zeros(shape=(1,len(groupdata[c + '_left'][0]['gaitcycle_ROM'])))
            ROM = pd.DataFrame(zero_data, columns=list(groupdata[c + '_left'][0]['gaitcycle_ROM'].keys()))
            ROMS = pd.DataFrame(zero_data, columns=list(groupdata[c + '_right'][0]['gaitcycle_ROM'].keys()))
            ROM_std = pd.DataFrame(zero_data, columns=list(groupdata[c + '_right'][0]['gaitcycle_ROM'].keys()))
            numoftrials = len(groupdata[c + '_left'])
            for i in range(0,numoftrials):
                df = pd.DataFrame([groupdata[c + '_left'][i]['gaitcycle_ROM']])
                ROM = ROM.add(df)
                ROMS = pd.concat([ROMS, df], ignore_index=True, sort=False)

            ROM = ROM.divide(numoftrials)
            ROMS = ROMS.iloc[1:]
            for (colnames, colvalues) in ROMS.items():
                ROM_std.iloc[0].at[colnames] = np.std(colvalues)

        avgdata[c + '_avg_ROM'] = ROM
        avgdata[c + '_std_ROM'] = ROM_std
            
    return avgdata