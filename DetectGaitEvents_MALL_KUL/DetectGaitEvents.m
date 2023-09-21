%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Gait Events: find forceplate gait cycle %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wouter Van Caekenberghe, November 2022
% Check README text file for more information

clearvars
close all
clc
tic

addpath(fullfile(cd,'functions'));
%% INPUT VARIABLES

crit_speed = 0.0001; % critical speed of heel marker
ge = 1;
Upperboundry = 150; % max possible duration of gaitcycle for filtering
Lowerboundry = 60; % min possible duration of gaitcycle for filtering
threshold = 10; % min Newton on forceplate to count as heel-strike 
RheelMarker = 'RHEE'; % name of right heel marker in .trc file
LheelMarker = 'LHEE'; % name of left heel marker in .trc file 

%% INITIALIZED VARIABLES

Gaitcycle = struct;
HeelStrike = struct;
h = char('TRC_HS_FP1', 'TRC_HS_FP3');
sides = ["Left" "Right"];

%% IMPORT DATA 

% trials need to be within one participant not over multiple participants
[TRCFilenames,path_TRC] = uigetfile('*.trc',  'Select trc files to process . . .',  'MultiSelect', 'on'); 
[MOTFilenames,path_MOT] = uigetfile('*.mot',  'Select mot files to process . . .',  'MultiSelect', 'on'); 

%% DIRECTORIES

path_out = fullfile(cd,'Results_Detect_Gait_Events');

if ~exist(path_out, 'dir')
    mkdir(path_out)
end

%% Detect Gait Events

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
        if (min(abs(Leftside- TRC_HS_FP1)))<=15 % Deze conditie is niet altijd juist 
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
end

save(fullfile(path_out, 'Result_struct_Gaitcycles'),'Gaitcycle')

trial = TRCFilenames';
Start_Index_Left = zeros(length(fieldnames(Gaitcycle)),1);
Stop_Index_Left = zeros(length(fieldnames(Gaitcycle)),1);
Start_Index_Right = zeros(length(fieldnames(Gaitcycle)),1);
Stop_Index_Right = zeros(length(fieldnames(Gaitcycle)),1);

for i = 1:length(TRCFilenames)
    name = num2str(cell2mat(TRCFilenames(i)));
    Y1 = Gaitcycle.(name(1:end-4)).Left;
    Q1 = Gaitcycle.(name(1:end-4)).Right;
    if any(Y1)
        Start_Index_Left(i,1) = Gaitcycle.(name(1:end-4)).Left(1);
        Stop_Index_Left(i,1) = Gaitcycle.(name(1:end-4)).Left(2);
    else
    end
    if any(Q1)
        Start_Index_Right(i,1) = Gaitcycle.(name(1:end-4)).Right(1);
        Stop_Index_Right(i,1) = Gaitcycle.(name(1:end-4)).Right(2);
    else
    end
end

T = table(trial,Start_Index_Right, Stop_Index_Right, Start_Index_Left, Stop_Index_Left);
writetable(T, (fullfile(path_out,'Gaitcycles_results.xlsx')))

t= toc;
disp(['%% Found all gaitcycles succesfully in ' num2str(floor(t/60)) ' minutes and ' num2str(rem(t,60)) ' seconds'])