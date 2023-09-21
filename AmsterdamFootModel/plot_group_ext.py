# Plot AFM function for groups

import matplotlib.pyplot as plt
import numpy as np
import savefig as save
import os 

def plot_group_AFM_R(xx, forefootanalysis, normalize, condition, avgdata, filename):
    
    footlabels = ['FOSK','HFSK','FFHF','HXFFM']
    midfootlabels = ['MFHF','FFMF','FFmMF','FFlMF']
    archlabels = ['MLA', 'TTA']
    planes = ['z','x','y']  

    if forefootanalysis==1:
        footlabels[3]='HXFF'
        midfootlabels = ['MFHF','FFMF']

    # Subplots: forefoot
    f1, ax = plt.subplots(4,3, figsize=(15,10))
    f1.subplots_adjust(hspace=0.7, wspace=0.3)
    f1.delaxes(ax[3,1])
    f1.suptitle('Amsterdam Foot Model - Right Foot \n' + condition, fontsize = 20)

    for i in range(0,len(xx)):
        if normalize==1:
            yy = xx[i]['gaitcycle_normalized']
            ixlabel = 'Gait cycle [%]'
        else:
            yy = xx[i]['gaitcycle']
            ixlabel = 'Gait cycle'
            # Find longest trial to set xlim values
            if i == 0:
                trial_length = list()
                for ii in range(0,len(xx)):
                    trial_length.append(len(xx[ii]['gaitcycle']['HFSKx']))
                maxlen = max(trial_length)
                numofpoints = list(range(1, maxlen))
        for row, seg in enumerate(footlabels):
            for col, plane in enumerate(planes):
                if normalize==1:
                    numofpoints = list(range(1,len(yy[seg + plane])+1))
                upperbound = np.add(avgdata[condition + '_avg'][seg + plane], avgdata[condition + '_std'][seg + plane])
                lowerbound = np.subtract(avgdata[condition + '_avg'][seg + plane], avgdata[condition + '_std'][seg + plane])
                if row==0 and col==0:
                    if normalize==1:
                        ax[row,col].fill_between(numofpoints, upperbound, lowerbound, facecolor = 'lightblue')
                    ax[row,col].plot(yy[seg + plane], label = filename[i])
                    ax[row,col].set_xlabel(ixlabel)
                    ax[row,col].set_xlim(numofpoints[0], numofpoints[-1])
                else:
                    if normalize==1:
                        ax[row,col].fill_between(numofpoints, upperbound, lowerbound, facecolor = 'lightblue')
                    ax[row,col].plot(yy[seg + plane])
                    ax[row,col].set_xlabel(ixlabel)
                    ax[row,col].set_xlim(numofpoints[0], numofpoints[-1])

    # Y labels
    for row in range(0,4):
        ax[row,0].set_ylabel('Pfl(-) / Dfl(+)', fontsize=9)
    for row in [0,2]:
        ax[row,1].set_ylabel('Ev(-) / Inv(+)', fontsize=9)
    for row in [1,2]:
        if row==1:
            col = 1
            ax[row,col].set_ylabel('Valg(-) / Var(+)', fontsize=9)
        elif row==2:
            col = 2
            ax[row,col].set_ylabel('Valg(-) / Var(+)', fontsize=9)
    for row in [0,1]:
            ax[row,2].set_ylabel('Exo(-) / Endo(+)', fontsize=9)
    ax[2,2].set_ylabel('Abd(-) / Add(+)', fontsize=9)

    # Row labels
    for ax, row in zip(ax[:,0], footlabels):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - 10, 0), xycoords=ax.yaxis.label, textcoords='offset points', size='large', ha='right', va='center')

    f1.tight_layout           
    f1.legend(loc = (0.45,0.12))
    f1.align_ylabels()
    plt.close()

    # Subplots: midfoot
    if forefootanalysis==1:
        numrows = 2
    elif forefootanalysis==2:
        numrows = 4
    f2, ax2 = plt.subplots(numrows,3, figsize=(15,10))
    f2.subplots_adjust(hspace=0.7, wspace=0.4)
    f2.suptitle('Amsterdam Foot Model - Right Midfoot \n' + condition, fontsize = 20)

    for i in range(0,len(xx)):
        if normalize==1:
            yy = xx[i]['gaitcycle_normalized']
        else:
            yy = xx[i]['gaitcycle']
            # Find longest trial to set xlim values
            if i == 0:
                trial_length = list()
                for i in range(0,len(xx)):
                    trial_length.append(len(xx[i]['gaitcycle']['HFSKx']))
                maxlen = max(trial_length)
                numofpoints = list(range(1, maxlen))            
        for row, seg in enumerate(midfootlabels):
            for col, plane in enumerate(planes):
                if normalize==1:
                    numofpoints = list(range(1,len(yy[seg + plane])+1))
                upperbound = np.add(avgdata[condition + '_avg'][seg + plane], avgdata[condition + '_std'][seg + plane])
                lowerbound = np.subtract(avgdata[condition + '_avg'][seg + plane], avgdata[condition + '_std'][seg + plane]) 
                if row==0 and col==0:
                    if normalize==1:
                        ax2[row,col].fill_between(numofpoints, upperbound, lowerbound, facecolor = 'lightblue')
                    ax2[row,col].plot(yy[seg + plane], label = filename[i])
                    ax2[row,col].set_xlabel(ixlabel)
                    ax2[row,col].set_xlim(numofpoints[0], numofpoints[-1])
                else:
                    if normalize==1:
                        ax2[row,col].fill_between(numofpoints, upperbound, lowerbound, facecolor = 'lightblue')
                    ax2[row,col].plot(yy[seg + plane])
                    ax2[row,col].set_xlabel(ixlabel)
                    ax2[row,col].set_xlim(numofpoints[0], numofpoints[-1])
        

    # Y labels
    for row in range(0,numrows):
        ax2[row,0].set_ylabel('Pfl(-) / Dfl(+)', fontsize=9)
        ax2[row,1].set_ylabel('Ev(-) / Inv(+)', fontsize=9)
        ax2[row,2].set_ylabel('Abd(-) / Add(+)', fontsize=9)

    # Row labels
    for ax2, row in zip(ax2[:,0], midfootlabels):
        ax2.annotate(row, xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - 10, 0), xycoords=ax2.yaxis.label, textcoords='offset points', size='large', ha='right', va='center')

    f2.tight_layout           
    f2.align_ylabels()
    plt.close()

    # Subplots: arches
    f3, (ax31, ax32) = plt.subplots(2, sharex=True, sharey=True, figsize=(10,5))

    upperboundmla = np.add(avgdata[condition + '_avg']['MLA'], avgdata[condition + '_std']['MLA'])
    lowerboundmla = np.subtract(avgdata[condition + '_avg']['MLA'], avgdata[condition + '_std']['MLA'])
    
    upperboundtta = np.add(avgdata[condition + '_avg']['TTA'], avgdata[condition + '_std']['TTA'])
    lowerboundtta = np.subtract(avgdata[condition + '_avg']['TTA'], avgdata[condition + '_std']['TTA'])

    for i in range(0,len(xx)):
        if normalize==1:
            yy = xx[i]['gaitcycle_normalized']
        else:
            yy = xx[i]['gaitcycle']

        # ax31
        if normalize==1:
            ax31.fill_between(numofpoints, upperboundmla, lowerboundmla, facecolor='lightblue')
        ax31.plot(yy['MLA'])
        ax31.set_ylabel('MLA', fontsize=12)
        ax31.set_ylim(80,140)
        ax31.set_xlim(0,len(yy['MLA'])-1)

        # ax32
        if normalize==1:
            ax32.fill_between(numofpoints, upperboundtta, lowerboundtta, facecolor='lightblue')
        ax32.plot(yy['TTA'])
        ax32.set_ylabel('TTA', fontsize=12)
    
    f3.suptitle('Amsterdam Foot Model - Right Arches \n' + condition, fontsize = 15)
    f3.supxlabel(ixlabel)
    plt.close()
    
    save.jpeg(f1,os.path.join('AFM_Kinematics',(condition + '_grouped')),'foot_right')
    save.jpeg(f2,os.path.join('AFM_Kinematics',(condition + '_grouped')),'midfoot_right')
    save.jpeg(f3,os.path.join('AFM_Kinematics',(condition + '_grouped')),'arches_right')



def plot_group_AFM_L(xx, forefootanalysis, normalize, condition, avgdata, filename):


    footlabels = ['FOSK','HFSK','FFHF','HXFFM']
    midfootlabels = ['MFHF','FFMF','FFmMF','FFlMF']
    archlabels = ['MLA', 'TTA']
    planes = ['z','x','y']

    if forefootanalysis==1:
        footlabels[3]='HXFF'
        midfootlabels = ['MFHF','FFMF']

    # Subplots: forefoot
    f1, ax = plt.subplots(4,3, figsize=(15,10))
    f1.subplots_adjust(hspace=0.7, wspace=0.3)
    f1.delaxes(ax[3,1])
    f1.suptitle('Amsterdam Foot Model - Left Foot \n' + condition, fontsize = 20)

    for i in range(0,len(xx)):
        if normalize==1:
            yy = xx[i]['gaitcycle_normalized']
            ixlabel = 'Gait cycle [%]'
        else:
            yy = xx[i]['gaitcycle']
            ixlabel = 'Gait cycle'
            # Find longest trial to set xlim values
            if i == 0:
                trial_length = list()
                for i in range(0,len(xx)):
                    trial_length.append(len(xx[i]['gaitcycle']['HFSKx']))
                maxlen = max(trial_length)
                numofpoints = list(range(1, maxlen))            
        for row, seg in enumerate(footlabels):
            for col, plane in enumerate(planes):
                if normalize==1:
                    numofpoints = list(range(1,len(yy[seg + plane])+1))
                upperbound = np.add(avgdata[condition + '_avg'][seg + plane], avgdata[condition + '_std'][seg + plane])
                lowerbound = np.subtract(avgdata[condition + '_avg'][seg + plane], avgdata[condition + '_std'][seg + plane])                
                if row==0 and col==0:
                    if normalize==1:
                        ax[row,col].fill_between(numofpoints, upperbound, lowerbound, facecolor = 'lightblue')
                    ax[row,col].plot(yy[seg + plane], label = filename[i])
                    ax[row,col].set_xlabel(ixlabel)
                    ax[row,col].set_xlim(numofpoints[0], numofpoints[-1])
                else:
                    if normalize==1:
                        ax[row,col].fill_between(numofpoints, upperbound, lowerbound, facecolor = 'lightblue')
                    ax[row,col].plot(yy[seg + plane])
                    ax[row,col].set_xlabel(ixlabel)
                    ax[row,col].set_xlim(numofpoints[0], numofpoints[-1])

    # Y labels
    for row in range(0,4):
        ax[row,0].set_ylabel('Pfl(-) / Dfl(+)', fontsize=9)
    for row in [0,2]:
        ax[row,1].set_ylabel('Ev(-) / Inv(+)', fontsize=9)
    for row in [1,2]:
        if row==1:
            col = 1
            ax[row,col].set_ylabel('Valg(-) / Var(+)', fontsize=9)
        elif row==2:
            col = 2
            ax[row,col].set_ylabel('Valg(-) / Var(+)', fontsize=9)
    for row in [0,1]:
            ax[row,2].set_ylabel('Exo(-) / Endo(+)', fontsize=9)
    ax[2,2].set_ylabel('Abd(-) / Add(+)', fontsize=9)

    # Row labels
    for ax, row in zip(ax[:,0], footlabels):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - 10, 0), xycoords=ax.yaxis.label, textcoords='offset points', size='large', ha='right', va='center')

    f1.tight_layout           
    f1.legend(loc = (0.45,0.12))
    f1.align_ylabels()
    plt.close()

    # Subplots: midfoot
    if forefootanalysis==1:
        numrows = 2
    elif forefootanalysis==2:
        numrows = 4
    f2, ax2 = plt.subplots(numrows,3, figsize=(15,10))
    f2.subplots_adjust(hspace=0.7, wspace=0.4)
    f2.suptitle('Amsterdam Foot Model - Left Midfoot \n' + condition, fontsize = 20)

    for i in range(0,len(xx)):
        if normalize==1:
            yy = xx[i]['gaitcycle_normalized']
        else:
            yy = xx[i]['gaitcycle']
            # Find longest trial to set xlim values
            if i == 0:
                trial_length = list()
                for i in range(0,len(xx)):
                    trial_length.append(len(xx[i]['gaitcycle']['HFSKx']))
                maxlen = max(trial_length)
                numofpoints = list(range(1, maxlen))            
        for row, seg in enumerate(midfootlabels):
            for col, plane in enumerate(planes):
                if normalize==1:
                    numofpoints = list(range(1,len(yy[seg + plane])+1))
                upperbound = np.add(avgdata[condition + '_avg'][seg + plane], avgdata[condition + '_std'][seg + plane])
                lowerbound = np.subtract(avgdata[condition + '_avg'][seg + plane], avgdata[condition + '_std'][seg + plane])                
                if row==0 and col==0:
                    if normalize==1:
                        ax2[row,col].fill_between(numofpoints, upperbound, lowerbound, facecolor = 'lightblue')
                    ax2[row,col].plot(yy[seg + plane], label = filename[i])
                    ax2[row,col].set_xlabel(ixlabel)
                    ax2[row,col].set_xlim(numofpoints[0], numofpoints[-1])
                else:
                    if normalize==1:
                        ax2[row,col].fill_between(numofpoints, upperbound, lowerbound, facecolor = 'lightblue')
                    ax2[row,col].plot(yy[seg + plane])
                    ax2[row,col].set_xlabel(ixlabel)
                    ax2[row,col].set_xlim(numofpoints[0], numofpoints[-1])

    # Y labels
    for row in range(0,numrows):
        ax2[row,0].set_ylabel('Pfl(-) / Dfl(+)', fontsize=9)
        ax2[row,1].set_ylabel('Ev(-) / Inv(+)', fontsize=9)
        ax2[row,2].set_ylabel('Abd(-) / Add(+)', fontsize=9)

    # Row labels
    for ax2, row in zip(ax2[:,0], midfootlabels):
        ax2.annotate(row, xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - 10, 0), xycoords=ax2.yaxis.label, textcoords='offset points', size='large', ha='right', va='center')

    f2.tight_layout           
    f2.align_ylabels()
    plt.close()

    # Subplots: arches
    f3, (ax31, ax32) = plt.subplots(2, sharex=True, sharey=True, figsize=(10,5))

    upperboundmla = np.add(avgdata[condition + '_avg']['MLA'], avgdata[condition + '_std']['MLA'])
    lowerboundmla = np.subtract(avgdata[condition + '_avg']['MLA'], avgdata[condition + '_std']['MLA'])
    
    upperboundtta = np.add(avgdata[condition + '_avg']['TTA'], avgdata[condition + '_std']['TTA'])
    lowerboundtta = np.subtract(avgdata[condition + '_avg']['TTA'], avgdata[condition + '_std']['TTA'])

    for i in range(0,len(xx)):
        if normalize==1:
            yy = xx[i]['gaitcycle_normalized']
        else:
            yy = xx[i]['gaitcycle']
    
        # ax31
        if normalize==1:
            ax31.fill_between(numofpoints, upperboundmla, lowerboundmla, facecolor='lightblue')
        ax31.plot(yy['MLA'])
        ax31.set_ylabel('MLA', fontsize=12)
        ax31.set_ylim(80,140)
        ax31.set_xlim(0,len(yy['MLA'])-1)

        # ax32
        if normalize==1:
            ax32.fill_between(numofpoints, upperboundtta, lowerboundtta, facecolor='lightblue')
        ax32.plot(yy['TTA'])
        ax32.set_ylabel('TTA', fontsize=12)

    f3.suptitle('Amsterdam Foot Model - Left Arches \n' + condition, fontsize = 15)
    f3.supxlabel(ixlabel)
    plt.close()

    save.jpeg(f1,os.path.join('AFM_Kinematics',(condition + '_grouped')),'foot_left')
    save.jpeg(f2,os.path.join('AFM_Kinematics',(condition + '_grouped')),'midfoot_left')
    save.jpeg(f3,os.path.join('AFM_Kinematics',(condition + '_grouped')),'arches_left')
