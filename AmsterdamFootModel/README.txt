The Amsterdam Foot Model is a multi-segment foot model to track foot kinematics. 
The author of the paper included a MatLab script to analyze the kinematics of 1 trial where the AFM was used. 
During an internship at a company, my supervisor and I decided that we wanted to use the AFM for an upcoming study. However, the company did not hold any MatLab subscriptions. Therefore, I took it upon myself to translate the entire kinematics script into python. During this process I also added a number of additional features. 1. The ability to process all trials for multiple conditions at once and return visualizations and excel data for the individual trials, inside conditions, and between conditions. 2. The ability to add an excel file containing the heel strike and toe-off timings. (The detect gait events script in this repository outputs the necessary excel input, but can also be done manually)


.py files containing 'ext' refer to the main script that processes all conditions together
(e.g. Dynamic_ext)

Main files are: AFM_script_extended (multiple conditions), AFM_Script (1 trial), ProofOfConcept



> Schallig, W., van den Noort, J.C., Piening, M. et al. The Amsterdam Foot Model: a clinically informed multi-segment foot model developed to minimize measurement errors in foot kinematics. J Foot Ankle Res 15, 46 (2022). https://doi.org/10.1186/s13047-022-00543-6
