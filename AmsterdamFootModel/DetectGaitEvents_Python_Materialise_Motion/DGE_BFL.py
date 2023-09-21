# Function to detect gaitevents (footcontact) in BigFootLab, Paal and write to excel file as input for AFM script.
# Wouter Van Caekenberghe, may, 2023

# Forceplate 1 = first forceplate, closest to the door
# Forceplate 2 = middle force plate
# Forceplate 3 = large force plate

from ezc3d import c3d
import numpy as np
from ButterFilterLow import LowpassFilter
import sys
import pandas as pd
import os   

def detect_gait_events(Dynamic_paths, conditions, wd, exname,filename, rhm, lhm):
    RHC = list()
    RTO = list()
    LHC = list()
    LTO = list()
    count = 0   
    for idx, cond in enumerate(conditions):
        for idxx, trial in enumerate(Dynamic_paths[cond]):
            trialname = filename[idx][idxx]
            print(trial)
            leftside_hc = list()
            leftside_to = list()
            rightside_hc = list()
            rightside_to = list()
            count += 1
            # Retrieve data from c3d file
            c = c3d(trial)
            analog_data = c['data']['analogs']
            point_data_dyn= c['data']['points']
            labels_dyn = list(c['parameters']['POINT']['LABELS']['value'])

            # Adjust name of heel marker if necessary
            if rhm in labels_dyn and lhm in labels_dyn:
                Leftheel = point_data_dyn[2,labels_dyn.index(lhm),:]
                Rightheel = point_data_dyn[2,labels_dyn.index(rhm),:]
            else:
                sys.exit('ERROR: heelmarkers not found in detect gait events script')
                                    
            Data = dict()                                                           
            Data['fp1'] = abs(analog_data[0,2,:])
            Data['fp2'] = abs(analog_data[0,8,:])
            Data['fp3'] = abs(analog_data[0,14,:])

            # Find heelstrikes and toe offs 
            Gaitevents = dict()
            for key in Data:
                logical = (Data[key] > 10) * 1
                difflist = np.diff(logical)
                hc = np.where(difflist==1)[0]
                to = np.where(difflist==-1)[0]

                if not to.size > 0:
                    to = np.array(len(logical) - 1)
                if not hc.size > 0:
                    hc = np.array(0)
                if hc.size > 1 and hc.size==to.size:
                    gaits = [[z,z1] for z, z1 in zip(hc,to)]
                    difflist = list()
                    for xx in gaits:
                        temp = xx[1]-xx[0]
                        difflist.append(temp)
                    m = max(difflist)
                    hc, to = gaits[difflist.index(m)]     
                if hc.size > 1:
                    hc = hc[0]
                if to.size > 1:
                    to = to[0]                

                Gaitevents[key] = [int(hc), int(to)]

            # Check left or right side for first forceplate contact
            for value in Gaitevents.values():
                if Rightheel[int(value[0]/10)] > Leftheel[int(value[0]/10)]:
                    leftside_hc.append(int(value[0]/10))
                    leftside_to.append(int(value[1]/10))
                else:
                    rightside_hc.append(int(value[0]/10))
                    rightside_to.append(int(value[1]/10))

            rightside_hc.sort()
            rightside_to.sort()
            leftside_hc.sort()
            leftside_to.sort()

            # Visualize
            from Extrapolate import norm_to_fp
            import matplotlib.pyplot as plt 

            plt.figure(figsize=(15,9))
            plt.plot(norm_to_fp(Data['fp1']), label = 'Forceplate 1', color='black', linestyle = 'dotted', alpha=0.2, zorder=1)
            plt.plot(norm_to_fp(Data['fp2']), label = 'Forceplate 2', color='black', alpha=0.2, zorder=1)
            plt.plot(norm_to_fp(Data['fp3']), label = 'Forceplate 3', color='black', linestyle = 'dashed', alpha=0.2, zorder=1)
            plt.plot(Rightheel, color = 'r', zorder=3, label = 'Right Heel Marker')
            plt.plot(Leftheel, color = 'g', zorder=3, label = 'Left Heel Marker')
            plt.title('DetecGaitEvents', pad=30, fontsize = 15)
            plt.xlabel('Time [0.01s]', labelpad=20)
            plt.ylabel('Newton [N]', labelpad=20)
            plt.ylim([0, max(Data['fp1'])+100])
            plt.plot(rightside_hc[0], Rightheel[rightside_hc[0]], marker='X', markerfacecolor='red', markeredgecolor='red', label = 'Right Heelstrike', markersize = 12, zorder=4, linestyle='None')
            plt.plot(rightside_to[0], Rightheel[rightside_to[0]], marker='D', markerfacecolor='red', markeredgecolor='red', label = 'Right Toe-off', markersize = 9, zorder=4, linestyle='None')
            plt.plot(leftside_hc[0], Leftheel[leftside_hc[0]], marker='X', markerfacecolor='green', markeredgecolor='green', label = 'Left Heelstrike', markersize = 12, zorder=4, linestyle='None')
            plt.plot(leftside_to[0], Leftheel[leftside_to[0]], marker='D', markerfacecolor='green', markeredgecolor='green', label = 'Left Toe-off', markersize = 9, zorder=4, linestyle='None')
            plt.legend(loc='upper right', ncols=2)
            plt.savefig(os.path.join(wd, 'Gait_Events', f'Visualizations_{exname}', trialname[:-4]))
            plt.close()
        
            # Take only first footstrike per side
            RHC.append(rightside_hc[0])
            RTO.append(rightside_to[0])
            LHC.append(leftside_hc[0])
            LTO.append(leftside_to[0])

            # Prepare data for excel
            dataR = dict()
            dataL = dict()
            
            dataR['HS'] = RHC
            dataR['TO'] = RTO
            dataL['HS'] = LHC
            dataL['TO'] = LTO

            Right = pd.DataFrame(dataR)
            Left = pd.DataFrame(dataL)

            # Write to excel
            with pd.ExcelWriter(os.path.join(wd, 'Gait_Events', f'{exname}.xlsx')) as writer:
                Right.to_excel(writer, sheet_name='Right')
                Left.to_excel(writer, sheet_name='Left')

