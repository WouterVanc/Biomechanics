%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% EMG data to Muscle Synergy Analysis input %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wouter Van Caekenberghe, November 2022
% Check README text file for more information

clearvars
close all
clc
tic

addpath(fullfile(cd,'functions'));
%% INPUT VARIABLES

subjectnr = '210'; 

% Part 1
crit_speed = 0.0001; % critical speed of heel marker
ge = 1;
Upperboundry = 150; % max possible duration of gaitcycle for filtering
Lowerboundry = 60; % min possible duration of gaitcycle for filtering
threshold = 10; % min Newton on forceplate to count as heel-strike 
RheelMarker = 'RHEE'; % name of right heel marker in .trc file
LheelMarker = 'LHEE'; % name of left heel marker in .trc file 

% Part 2 
Outside_lab = 0;
% Define the corresponding muscles  IN RESPECTIVE ORDER
channels_used = [1 2 3 4 5 6 7 8 17 18 19 20 21 22 24 25 26 27 28 29]; %1-16 = White system 17-32 = Blue system
muscles = {'GlutMedL', 'TibAntL' ,'GastMedL' ,'GastLatL', 'SoleusL', 'GlutMedR', 'TibAntR', 'GastMedR', 'RecFemL', 'VastLatL', 'VastMedL', 'HamLatL', 'HamMedL', 'RecFemR', 'VastLatR', 'HamLatR', 'HamMedR', 'GastLatR', 'SoleusR', 'VastMedR'};% Define the corresponding muscles  IN RESPECTIVE ORDER 

%Define filter properties
order = 4; 
cutoff_band = [20 400];
cutoff_low = 6;

% Part 3
NormalizeToMaxDuringTrials = 0; % Put 1 if you want to Normalize, 0 if you don't
RemoveChannel = 0; % Put 1 if you want to remove, 0 if you don't = optional
UnwantedMuscle = 'TeFaLaR '; % Channel/Muscle you want to remove = optional
UnwantedMuscle2 = 'TeFaLaL '; % Channel/Muscle you want to remove = optional

%% INITIALIZED VARIABLES

HeelStrike = struct;
Gaitcycle = struct;
h = char('TRC_HS_FP1', 'TRC_HS_FP3');
sides = ["Left" "Right"];
Time = zeros(1,101);
Time(1,:) = 1:101;
DATA_LST = char(muscles);
nmus = length(muscles);
EMG_results_R = zeros(101,nmus);
EMG_results_L = zeros(101,nmus);
str = {'Normalized_' 'NOT_Normalized_'};

%% DIRECTORIES

output_path = fullfile(cd,'Temporary_Storage'); % Temporary folder 
if ~exist(output_path, 'dir')
    mkdir(output_path)
end

if NormalizeToMaxDuringTrials==1
    path_out = fullfile(cd,['Results_EMG_Processing_' char(str(1)) num2str(subjectnr)]);
    if~exist(path_out, 'dir')
        mkdir(path_out);
    end
elseif NormalizeToMaxDuringTrials==0
    path_out = fullfile(cd,['Results_EMG_Processing_' char(str(2)) num2str(subjectnr)]);
    if~exist(path_out, 'dir')
        mkdir(path_out);
    end
end

%% IMPORT DATA 

% trials need to be within one participant not over multiple participants
[TRCFilenames,path_TRC] = uigetfile('*.trc',  'Select trc files to process . . .',  'MultiSelect', 'on'); 
[MOTFilenames,path_MOT] = uigetfile('*.mot',  'Select mot files to process . . .',  'MultiSelect', 'on'); 
[c3dFilenames,path_C3D] = uigetfile('*.c3d',  'Select c3d files to process . . .',  'MultiSelect', 'on'); 

%% Part 1: Detect Gait Events

for j = 1:length(TRCFilenames)   
    
    % Load in the .trc data
    filename= fullfile(path_TRC , TRCFilenames(j));
    filename = char(filename);
    v = read_trcFile(filename);
    name = num2str(cell2mat(TRCFilenames(j)));

    % Right side
    Rheel_index = find(contains(v.MarkerList,RheelMarker));
    Rheel_index_TRC = ((Rheel_index)*3); % initial two columns are time and then every marker has 3 columns for x,y,z
    Rheel = v.Data(:,Rheel_index_TRC+1);
    
    % HS_identifier_right
    ende_right=find(((Rheel(ge+1:end-3)-Rheel(ge:end-4))<-crit_speed) .*...
        ((Rheel(ge+2:end-2)-Rheel(ge+1:end-3))<-crit_speed) .*...
        ((Rheel(ge+3:end-1)-Rheel(ge+2:end-2))<-crit_speed) .*...
        ((Rheel(ge+4:end)-Rheel(ge+3:end-1)>-crit_speed)).*...
        ((Rheel(ge:end-4)<abs(47))))+(ge+2);
    
    % Remove faulty asterisks 
    for q = 1:length(ende_right)-1
        if (min(abs(ende_right(q+1,1)-ende_right(q,1))))<=25
            ende_right(q)= 0;
        else
        end
    end

    ende_right = nonzeros(ende_right);
    
    % Find actual gait cycles
    for r = 1:length(ende_right)-1
        if (ende_right(r+1,1)-ende_right(r,1))<Upperboundry && (ende_right(r+1,1)-ende_right(r,1))>Lowerboundry
            ende_right_EMG(r,:)=[ende_right(r,1), ende_right(r+1,1)];
        else %do nothing
        end
    end  
    
    ind = ende_right_EMG(:,1)==0;
    ende_right_EMG(ind,:)=[];
    
    % Plot the gaitcycles and asterisks
    figure(1)
    hold on
    plot(Rheel)
    plot(ende_right_EMG,Rheel(ende_right_EMG),'*')
    title('Right gaitcycles')
    grid on
    hold on
    
    % Left side
    Lheel_index = find(contains(v.MarkerList,LheelMarker));
    Lheel_index_TRC = ((Lheel_index)*3); % initial two columns are time and then every marker has 3 columns for x,y,z
    Lheel = v.Data(:,Lheel_index_TRC+1);

    % HS_identifier_left
    ende_left=find(((Lheel(ge+1:end-3)-Lheel(ge:end-4))<-crit_speed) .*...
        ((Lheel(ge+2:end-2)-Lheel(ge+1:end-3))<-crit_speed) .*...
        ((Lheel(ge+3:end-1)-Lheel(ge+2:end-2))<-crit_speed) .*...
        ((Lheel(ge+4:end)-Lheel(ge+3:end-1)>-crit_speed)).*...
        ((Lheel(ge:end-4)<abs(58))))+(ge+2);

    % Remove faulty asterisks
    for q = 1:length(ende_left)-1
        if (min(abs(ende_left(q+1,1)-ende_left(q,1))))<=25
            ende_left(q)= 0;
        else
        end
    end

    ende_left = nonzeros(ende_left);

    % Find actual gait cycles
    for w = 1:length(ende_left)-1
        if (ende_left(w+1,1)-ende_left(w,1))<Upperboundry && (ende_left(w+1,1)-ende_left(w,1))>Lowerboundry
            ende_left_EMG(w,:)=[ende_left(w,1), ende_left(w+1,1)];
        else 
        end
    end
    
    ind = find(ende_left_EMG(:,1)==0);
    ende_left_EMG(ind,:)=[];
    
    % Plot the gaitcycles and asterisks
    figure(2)
    hold on
    plot(Lheel)
    plot(ende_left,Lheel(ende_left),'*')
    title('Left gaitcycles')
    grid on
    hold on

    % Save data in struct
    HeelStrike.left.(name(1:end-4)) = ende_left_EMG;
    HeelStrike.right.(name(1:end-4)) = ende_right_EMG;

    clear ende_right ende_left ende_right_EMG ende_left_EMG

% Match FP contact in .mot file to gaitcycle from .trc
    
    % Load in the .mot data
    filenameMOT = append(path_MOT,MOTFilenames(j));
    filenameMOT = char(filenameMOT);
    MOT_data = ReadMotFile(filenameMOT);
    
    % Define ForcePlates
    ForcePlate1 = MOT_data.data(:,contains(MOT_data.names, 'ground_force_vy') & ~contains(MOT_data.names, '1'));
    ForcePlate3 = MOT_data.data(:,contains(MOT_data.names, '1_ground_force_vy'));
    LogicFP1 = ForcePlate1>threshold;
    LogicFP3 = ForcePlate3>threshold;
    
    % Check if there is contact with each forceplate
    if any(LogicFP1)
        HeelstrikeFP1 = find(diff(LogicFP1));
        HeelstrikeFP1 = HeelstrikeFP1(1);
        TimeHS_FP1 = round(MOT_data.data(HeelstrikeFP1, contains(MOT_data.names, 'time')),2);
        TRC_HS_FP1 = find(v.Data(:,2)==TimeHS_FP1);
    else
    end
    
    if any(LogicFP3)
        HeelstrikeFP3 = find(diff(LogicFP3));
        HeelstrikeFP3 = HeelstrikeFP3(1);
        TimeHS_FP3 = round(MOT_data.data(HeelstrikeFP3, contains(MOT_data.names, 'time')),2);
        TRC_HS_FP3 = find(v.Data(:,2)==TimeHS_FP3);
    else
    end
    
    % Match heel-strike .mot data with .trc data
    Rightside = HeelStrike.right.(name(1:end-4));
    Leftside = HeelStrike.left.(name(1:end-4));
    if exist (h(1,:), 'var') && exist (h(2,:), 'var')
        if (min(abs(Leftside- TRC_HS_FP1)))<=15
            idx_Left = find(TRC_HS_FP1+5 > Leftside(:,1) & TRC_HS_FP1+5 < Leftside(:,2));
            idx_Right = find(TRC_HS_FP3+5 > Rightside(:,1) & TRC_HS_FP3+5 < Rightside(:,2));
        else 
            idx_Left = find(TRC_HS_FP3+5 > Leftside(:,1) & TRC_HS_FP3+5 < Leftside(:,2));
            idx_Right = find(TRC_HS_FP1+5 > Rightside(:,1) & TRC_HS_FP1+5 < Rightside(:,2));
        end
    
    Gaitcycle.(name(1:end-4)).Right = Rightside(idx_Right,:);
    Gaitcycle.(name(1:end-4)).Left = Leftside(idx_Left,:);
    
    elseif ~exist (h(1,:), 'var') && exist (h(2,:), 'var')
        if (min(abs(Leftside- TRC_HS_FP3)))<=15
            idx_Left = find(TRC_HS_FP3+5 > Leftside(:,1) & TRC_HS_FP3+5 < Leftside(:,2));
            Gaitcycle.(name(1:end-4)).Left = Leftside(idx_Left,:);
            Gaitcycle.(name(1:end-4)).Right = [];
        else
            idx_Right = find(TRC_HS_FP3+5 > Rightside(:,1) & TRC_HS_FP3+5 < Rightside(:,2));
            Gaitcycle.(name(1:end-4)).Right = Rightside(idx_Right,:);
            Gaitcycle.(name(1:end-4)).Left = [];
        end
    
    elseif exist (h(1,:), 'var') && ~exist (h(2,:), 'var')
        if (min(abs(Leftside- TRC_HS_FP1)))<=15
            idx_Left = find(TRC_HS_FP1+5 > Leftside(:,1) & TRC_HS_FP1+5 < Leftside(:,2));
            Gaitcycle.(name(1:end-4)).Left = Leftside(idx_Left,:);
            Gaitcycle.(name(1:end-4)).Right = [];
        else
            idx_Right = find(TRC_HS_FP1+5 > Rightside(:,1) & TRC_HS_FP1+5 < Rightside(:,2));
            Gaitcycle.(name(1:end-4)).Right = Rightside(idx_Right,:);
            Gaitcycle.(name(1:end-4)).Left = [];
        end
    end

%% Part 2 Process EMG data
    
    % Load .c3d data
    if iscell(c3dFilenames)
    filename = char(c3dFilenames(j));
    elseif ischar(c3dFilenames)
    filename = char(c3dFilenames);  
    end   
        if Outside_lab == 0
        [Markers,MLabels,VideoFrameRate,AnalogSignals,ALabels,AUnits,AnalogFrameRate,Event,ParameterGroup,CameraInfo]=...
            readC3D(char(strcat(path_C3D,filename))); %import c3d file       
        else 
        [AnalogSignals,ALabels,AnalogFrameRate]=...
            readc3dEMGRemote(char(strcat(path_C3D,filename)));
        end
    
    % Safety measure for missing data 
    Rightside2 = Gaitcycle.(name(1:end-4)).Right;
    Leftside2 = Gaitcycle.(name(1:end-4)).Left;
    if any(Leftside2) & any(Rightside2)
       q1=0;
       q2=2;
       w = j+length(c3dFilenames);
    elseif ~any(Leftside2) & any(Rightside2) 
       q1=1;
       q2=2;
       w = j;
    elseif any(Leftside2) & ~any(Rightside2)
       q1=0;
       q2=1;
       w = j;
    end
    
    for k = 1+q1:q2
        
       %get EMG-data out of the Analog signals
       [EMG_channel_white, EMG_channel_blue] = determine_EMG(ALabels, AnalogSignals);
    
       EMG_raw = [EMG_channel_white, EMG_channel_blue];

       % Filtering
        
       % band pass filter;
       [a,b] = butter(order/2,cutoff_band./(0.5*AnalogFrameRate),'bandpass');
       EMG_band = filtfilt(a,b,EMG_raw);
        
       clear a b

       %Rectify the band pass filtered signal
       EMG_rect = abs(EMG_band);
        
       %low pass filter the rectified signal
       [a,b] = butter(order,cutoff_low./(0.5*AnalogFrameRate),'low');
       EMG_low = filtfilt(a,b,EMG_rect);
        
       %rearange EMG - channels
       EMG_out = zeros(length(EMG_raw),size(channels_used,1));        
       for channel = 1:size(channels_used,2)
                EMG_out(:,channel) = EMG_low(:,channels_used(channel));          
       end
        
       % Get GaitCycle and normalize to 101 datapoints
       ST = Gaitcycle.(name(1:end-4)).(sides(k));
       start = ST(1)/100;
       stop = ST(2)/100;
       timelist = (1:length(EMG_raw))';
       timelist = timelist./1000;
       start_ind = find(timelist(:,1)==start);
       stop_ind = find(timelist(:,1)==stop);
       EMG_out = EMG_out(start_ind:stop_ind,:);
       EMG_norm = TimeNormalize(EMG_out);
       if k==1
           EMG_out_L = EMG_norm;
           save(fullfile(output_path, (['Results_' num2str(j)])), 'EMG_out_L', "DATA_LST", 'Time') 
       elseif k==2
           EMG_out_R = EMG_norm;
           save(fullfile(output_path, (['Results_' num2str(w)])), 'EMG_out_R', 'DATA_LST', 'Time') 
       end
       clear EMG_out EMG_channel_white EMG_channel_blue
    end
end
%% Part 3 Combine Trial Results  
c = 0;
lcount = 0;
rcount = 0;
input = dir(output_path);
input = {input(3:end).name};
for x = 1:length(input)
    c = c+1;
    data = load (['./Temporary_Storage/' char(input(x))]);  
    LOR = fieldnames(data);
    if contains(LOR(2,1), 'R')
        rcount = rcount+1;
        EMG_results_R = EMG_results_R + data.EMG_out_R;
        EMG_MAX_R(c,:) = max(data.EMG_out_R);
    else
        lcount = lcount+1;
        EMG_results_L = EMG_results_L + data.EMG_out_L;
        EMG_MAX_L(c,:) = max(data.EMG_out_L);
    end
end

rmdir Temporary_Storage\ s

EMG_results_R = EMG_results_R./rcount;
EMG_results_L = EMG_results_L./lcount;

% Delete unwanted channels
MuscleNames = string(data.DATA_LST)';
if RemoveChannel==1
    
    
    EMG_results_R(:,MuscleNames(1,:)==UnwantedMuscle | MuscleNames(1,:)==UnwantedMuscle2) = [];
    EMG_results_L(:,MuscleNames(1,:)==UnwantedMuscle | MuscleNames(1,:)==UnwantedMuscle2) = [];
    
    EMG_MAX_R(:,MuscleNames(1,:)==UnwantedMuscle | MuscleNames(1,:)==UnwantedMuscle2) = [];
    EMG_MAX_L(:,MuscleNames(1,:)==UnwantedMuscle | MuscleNames(1,:)==UnwantedMuscle2) = [];
    
    DATA_LST = data.DATA_LST;
    DATA_LST(MuscleNames(1,:)==UnwantedMuscle | MuscleNames(1,:)==UnwantedMuscle2,:) = [];
    MuscleNames(MuscleNames(1,:)==UnwantedMuscle | MuscleNames(1,:)==UnwantedMuscle2) = [];
end

% Normalize to max of trials
if NormalizeToMaxDuringTrials == 1
    Norm_MAX_R = max(EMG_MAX_R);
    Norm_MAX_L = max(EMG_MAX_L);
    [x,y] = size(EMG_results_R);

    for i = 1:y
        EMG_results_R(:,i) = EMG_results_R(:,i)./Norm_MAX_R(1,i);
        EMG_results_L(:,i) = EMG_results_L(:,i)./Norm_MAX_L(1,i);
    end
end

% saving right side 
DATAsrc = EMG_results_R';
time = data.Time';

save(fullfile(path_out, 'Synergy_Data_R'), 'DATAsrc', 'DATA_LST', 'time', 'MuscleNames', 'subjectnr')

% saving left side 
DATAsrc = EMG_results_L';
time = data.Time';

save(fullfile(path_out, 'Synergy_Data_L'),'DATAsrc', 'DATA_LST', 'time', 'MuscleNames', 'subjectnr')

%% Part 4 Run Muscle Synergy Analysis for Both Sides

app = PosturalData_NMFvsPCA_GUI_July2013;
run PosturalData_NMFvsPCA_GUI_July2013.m
while isvalid(app), pause(0.1); end
options.WindowStyle = 'modal';
options.Interpreter = 'tex';
f = msgbox('Please load in data of the OPPOSITE side now',options);
while isvalid(f), pause(0.1); end
app2 = PosturalData_NMFvsPCA_GUI_July2013;
run PosturalData_NMFvsPCA_GUI_July2013.m
while isvalid(app2), pause(0.1); end
f = msgbox('When all subjects are processed, run the KmeansClustering script');


t = toc;
disp(['%% Finished processing succesfully in ' num2str(floor(t/60)) ' minutes and ' num2str(rem(t,60)) ' seconds'])
