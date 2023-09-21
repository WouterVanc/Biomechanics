
  _____       _            _      _____       _ _     ______               _        
 |  __ \     | |          | |    / ____|     (_) |   |  ____|             | |       
 | |  | | ___| |_ ___  ___| |_  | |  __  __ _ _| |_  | |____   _____ _ __ | |_ ___  
 | |  | |/ _ \ __/ _ \/ __| __| | | |_ |/ _` | | __| |  __\ \ / / _ \ '_ \| __/ __| 
 | |__| |  __/ ||  __/ (__| |_  | |__| | (_| | | |_  | |___\ V /  __/ | | | |_\__ \ 
 |_____/ \___|\__\___|\___|\__|  \_____|\__,_|_|\__| |______\_/ \___|_| |_|\__|___/ 
                                                                                    
 
Wouter Van Caekenberghe, November 2022
KU Leuven - FaBer - Human Movement Biomechanics Research Group
contact: wouter.vancaekenberghe@student.kuleuven.be

This folder contains a script that finds the gait cycles (indices) that touched the forceplate, utilizing *.trc en *.mot files, and writes them to an excel file.
There's one main script('DetectGaitEvents') and a number of functions in the 'functions' folder.


*** DetectGaitEvents ***

1) Modify the input variables to your needs ('INPUT VARIABLES' line 15-21).
2) Run the script 
3) Load *.trc and *.mot file of multiple trials of one subject
4) Output folder in current directory

INPUT --> *.trc and *.mot files of multiple trials of one subject

	(These should only be valid overground walking trials = proper contact with forceplate of atleast 1 foot) 
	(.mot files MUST contain BOTH ground_force_vy and 1_ground_force_vy, even if one is empty)

OUTPUT --> results folder ('Results_Detect_Gait_Events') containing excel file and *.mat file with start and stop indeces of the gait cycle for each trial.
	     Note that the results are indices of the *.trc files. If you want to know the exact time, look up the index in the *.trc file or more roughly divide index by 100. 	 
 
                                                                                  
