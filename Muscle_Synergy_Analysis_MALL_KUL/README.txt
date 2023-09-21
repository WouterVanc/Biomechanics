
  __  __                        _             _____                                                                              _                 _       
 |  \/  |                      | |           / ____|                                                     /\                     | |               (_)      
 | \  / |  _   _   ___    ___  | |   ___    | (___    _   _   _ __     ___   _ __    __ _   _   _       /  \     _ __     __ _  | |  _   _   ___   _   ___ 
 | |\/| | | | | | / __|  / __| | |  / _ \    \___ \  | | | | | '_ \   / _ \ | '__|  / _` | | | | |     / /\ \   | '_ \   / _` | | | | | | | / __| | | / __|
 | |  | | | |_| | \__ \ | (__  | | |  __/    ____) | | |_| | | | | | |  __/ | |    | (_| | | |_| |    / ____ \  | | | | | (_| | | | | |_| | \__ \ | | \__ \
 |_|  |_|  \__,_| |___/  \___| |_|  \___|   |_____/   \__, | |_| |_|  \___| |_|     \__, |  \__, |   /_/    \_\ |_| |_|  \__,_| |_|  \__, | |___/ |_| |___/
                                                       __/ |                         __/ |   __/ |                                    __/ |                
                                                      |___/                         |___/   |___/                                    |___/                 
Wouter Van Caekenberghe, November 2022
KU Leuven - FaBer - Human Movement Biomechanics Research Group
contact: wouter.vancaekenberghe@student.kuleuven.be

This folder contains a workflow to perform a muscle synergy analysis starting from EMG and MOCAP data.
There are two main scripts('EMG_to_MSA' and 'KmeansClustering') and a number of functions in the 'functions' folder.

Short tutorial

The workflow is designed for minimal input/work for the user. Run the EMG_to_MSA script for every subject, filling in the right input variables each time beforehand. 
Once a group is finished, run the KmeansClustering script, filling in the right input variables beforehand. Leave the structure of the directories/names in this folder the same.  

Full tutorial

*** EMG_to_MSA *** 

	1) Modify the input variables to your needs ('INPUT VARIABLES' line 13-40). Part 1 concerns the gait cycle detection, Part 2 the EMG filtering and Part 3 combining the trials and normalizing
 
	2) Run the program

		INPUT --> Load in *.trc *.mot *.c3d files of all the trials you want to process for ONE subject in the pop-up windows.
 
		(These should only be valid overground walking trials = proper contact with forceplate of atleast 1 foot) 
		(.mot files MUST contain BOTH ground_force_vy and 1_ground_force_vy, even if one is empty)

		OUTPUT -->  'Results_EMG_Processing...' folders in current directory, this is an intermediate ouput required for the muscle synergy analysis GUI

		GUI --> For the Muscle Synergy Analysis a GUI will pop-up automatically. Follow these steps
				
				1) Check the 'walking data' box
				2) Press load data (Select ONE side of the synergy data in the 'Results_EMG_Processing...' folder of your current subject)
				3) Press plot data
				4) Press Extract NMF muscle synergies
				5) Wait untill the 'How many synergies'-box shows a number (optional: plot the NMF)
				6) Close the GUI by pressing the x in the top right corner
				7) The GUI will pop-up again automatically. Please repeat steps 1-7 but this time for the synergy data of the OTHER side from the SAME subject

		OUTPUT -->  'Output_Matrix_NMF...' folders in current directory, containing synergy data in *.mat files for Left and Right side
				

		(All output directories will be made according to your preferences from the input variables and are used in the following script. Don't change these!) 

*** KmeansClustering ***

		INPUT --> Use the KmeansClustering script once all subjects of a group are processed with the previous script = output folders EMG_to_MSA

	1) Modify the input variables to your needs ('INPUT VARIABLES' line 12-16). 

	2) Run the script (For each group, for Left and Right)

		OUTPUT --> 'KmeansClustering_Results...' folder in current directory, containing two figures and the kmeans result data in *.mat file





The muscle synergy analysis GUI used in this workflow belongs to Ting and Chvatal, 2010. The README of their GUI can be found in the 'functions' folder.


%%% References %%%

> Ting, L.H.*, Chvatal, S.A. (2010) “Decomposing muscle activity in motor tasks: methods and interpretation”. In Motor Control: Theories, Experiments, and Applications, Danion, F., Latash, M.L. (eds), Oxford University Press, New York.



  