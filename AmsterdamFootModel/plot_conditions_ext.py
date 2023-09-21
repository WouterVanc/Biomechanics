# Plot AFM function for conditions

import matplotlib.pyplot as plt
import numpy as np
import savefig as save
import os 

def plot_conditions_AFM_R(avgdata, forefootanalysis, conditions):

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
    f1.suptitle('Amsterdam Foot Model - Right Foot', fontsize = 20)

    for condition in conditions:
        for row, seg in enumerate(footlabels):
            for col, plane in enumerate(planes):
                numofpoints = list(range(1,len(avgdata[condition + '_avg']['FOSKx']+1)))
                if row==0 and col==0:
                    ax[row,col].plot(avgdata[condition + '_avg'][seg + plane], label = condition)
                    ax[row,col].set_xlabel('Gait cycle [%]')
                    ax[row,col].set_xlim(numofpoints[0], numofpoints[-1])
                else:
                    ax[row,col].plot(avgdata[condition + '_avg'][seg + plane])
                    ax[row,col].set_xlabel('Gait cycle [%]')
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
    f1.legend(loc = (0.45,0.12), ncol=2)
    f1.align_ylabels()
    plt.close()

    # Subplots: midfoot
    if forefootanalysis==1:
        numrows = 2
    elif forefootanalysis==2:
        numrows = 4
    f2, ax2 = plt.subplots(numrows,3, figsize=(15,10))
    f2.subplots_adjust(hspace=0.7, wspace=0.4)
    f2.suptitle('Amsterdam Foot Model - Right Midfoot' , fontsize = 20)

    for condition in conditions:           
        for row, seg in enumerate(midfootlabels):
            for col, plane in enumerate(planes):
                numofpoints = list(range(1,len(avgdata[condition + '_avg']['MFHFx']+1))) 
                if row==0 and col==0:
                    ax2[row,col].plot(avgdata[condition + '_avg'][seg + plane], label = condition)
                    ax2[row,col].set_xlabel('Gait cycle [%]')
                    ax2[row,col].set_xlim(numofpoints[0], numofpoints[-1])
                else:
                    ax2[row,col].plot(avgdata[condition + '_avg'][seg + plane], label = condition)
                    ax2[row,col].set_xlabel('Gait cycle [%]')
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

    for condition in conditions:
        # ax31
        ax31.plot(avgdata[condition + '_avg']['MLA'])
        ax31.set_ylabel('MLA', fontsize=12)
        ax31.set_ylim(80,140)
        ax31.set_xlim(0,len(avgdata[condition + '_avg']['MLA'])-1)

        # ax32
        ax32.plot(avgdata[condition + '_avg']['TTA'])
        ax32.set_ylabel('TTA', fontsize=12)
    
    f3.suptitle('Amsterdam Foot Model - Right Arches', fontsize = 15)
    f3.supxlabel('Gait cycle [%]')
    plt.close()
    
    save.jpeg(f1,os.path.join('AFM_Kinematics', 'Compare_Conditions'),'foot_right')
    save.jpeg(f2,os.path.join('AFM_Kinematics','Compare_Conditions'),'midfoot_right')
    save.jpeg(f3,os.path.join('AFM_Kinematics','Compare_Conditions'),'arches_right')



def plot_conditions_AFM_L(avgdata, forefootanalysis, conditions):


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
    f1.suptitle('Amsterdam Foot Model - Left Foot', fontsize = 20)

    for condition in conditions:
        for row, seg in enumerate(footlabels):
            for col, plane in enumerate(planes):
                numofpoints = list(range(1,len(avgdata[condition + '_avg']['FOSKx']+1)))
                if row==0 and col==0:
                    ax[row,col].plot(avgdata[condition + '_avg'][seg + plane], label = condition)
                    ax[row,col].set_xlabel('Gait cycle [%]')
                    ax[row,col].set_xlim(numofpoints[0], numofpoints[-1])
                else:
                    ax[row,col].plot(avgdata[condition + '_avg'][seg + plane])
                    ax[row,col].set_xlabel('Gait cycle [%]')
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
    f1.legend(loc = (0.45,0.12), ncol=2)
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

    for condition in conditions:           
        for row, seg in enumerate(midfootlabels):
            for col, plane in enumerate(planes):
                numofpoints = list(range(1,len(avgdata[condition + '_avg']['MFHFx']+1))) 
                if row==0 and col==0:
                    ax2[row,col].plot(avgdata[condition + '_avg'][seg + plane], label = condition)
                    ax2[row,col].set_xlabel('Gait cycle [%]')
                    ax2[row,col].set_xlim(numofpoints[0], numofpoints[-1])
                else:
                    ax2[row,col].plot(avgdata[condition + '_avg'][seg + plane], label = condition)
                    ax2[row,col].set_xlabel('Gait cycle [%]')
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

    for condition in conditions:
        # ax31
        ax31.plot(avgdata[condition + '_avg']['MLA'])
        ax31.set_ylabel('MLA', fontsize=12)
        ax31.set_ylim(80,140)
        ax31.set_xlim(0,len(avgdata[condition + '_avg']['MLA'])-1)

        # ax32
        ax32.plot(avgdata[condition + '_avg']['TTA'])
        ax32.set_ylabel('TTA', fontsize=12)

    f3.suptitle('Amsterdam Foot Model - Left Arches \n' + condition, fontsize = 15)
    f3.supxlabel('Gait cycle [%]')
    plt.close()

    save.jpeg(f1,os.path.join('AFM_Kinematics','Compare_Conditions'),'foot_left')
    save.jpeg(f2,os.path.join('AFM_Kinematics','Compare_Conditions'),'midfoot_left')
    save.jpeg(f3,os.path.join('AFM_Kinematics','Compare_Conditions'),'arches_left')
