# Detect gait events, BigFootLab Materialise Motion - Paal
# Wouter Van Caekenberghe, 27/04/2023

from ezc3d import c3d
import numpy as np
from Extrapolate import norm_to_fp
import matplotlib.pyplot as plt

path = (r'C:\Users\WVanCaeckenberge\Desktop\Pilot_1_AFM_orig\HeelRaise\Dynamic\HeelRaise_Dynamic_01.c3d')
c = c3d(path)

# Forceplate 1 = first forceplate, closest to the door
# Forceplate 2 = middle force plate
# Forceplate 3 = large force plate

analog_data = c['data']['analogs']
point_data_dyn= c['data']['points']
labels_dyn = list(c['parameters']['POINT']['LABELS']['value'])
first_frame = c['header']['points']['first_frame']

# Adjust name of heel marker if necessary 
Leftheel = point_data_dyn[2,labels_dyn.index('LCALD'),:]
Rightheel = point_data_dyn[2,labels_dyn.index('RCALD'),:]

Data = dict()
Data['fp1'] = abs(analog_data[0,2,:])
Data['fp2'] = abs(analog_data[0,8,:])
Data['fp3'] = abs(analog_data[0,14,:])

Gaitevents = dict()
for key in Data:
    logical = (Data[key] > 10) * 1
    difflist = np.diff(logical)
    hc = np.where(difflist==1)[0]
    to = np.where(difflist==-1)[0]
    if not to.size > 0:
        to = len(logical) - 1

    Gaitevents[key] = [int(hc), int(to)]

# Check left or right side for first forceplate contact
leftside_hc = list()
leftside_to = list()
rightside_hc = list()
rightside_to = list()
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
plt.figure(figsize=(15,15))
plt.plot(norm_to_fp(Data['fp1']), label = 'Forceplate 1', color='black', linestyle = 'dotted', alpha=0.2, zorder=1)
plt.plot(norm_to_fp(Data['fp2']), label = 'Forceplate 2', color='black', alpha=0.2, zorder=1)
plt.plot(norm_to_fp(Data['fp3']), label = 'Forceplate 3', color='black', linestyle = 'dashed', alpha=0.2, zorder=1)
plt.plot(Rightheel, color = 'r', zorder=3, label = 'Right Heel Marker')
plt.plot(Leftheel, color = 'g', zorder=3, label = 'Left Heel Marker')
plt.title('DetecGaitEvents', pad=30, fontsize = 15)
plt.xlabel('Time [s]', labelpad=20)
plt.ylabel('Newton [N]', labelpad=20)
plt.ylim([0, max(Data['fp1'])+100])
plt.vlines((Gaitevents['fp1'][0]/10,Gaitevents['fp2'][0]/10,Gaitevents['fp3'][0]/10), ymin = 0, ymax = max(Data['fp1'])+100, colors='k', label='Heelstrikes', zorder=4)
plt.legend(loc='upper right')

plt.show()