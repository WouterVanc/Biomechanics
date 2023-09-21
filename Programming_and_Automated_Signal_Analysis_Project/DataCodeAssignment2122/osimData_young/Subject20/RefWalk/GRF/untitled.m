clear all; close all;

% Import data

Gait_data = importdata('Refwalk_GRF.mot');

% % Extract time + vertical force of left and right leg
% 
% time = Gait_data.data(:,1);
% GRF_left_vy = Gait_data.data(:,3);
% GRF_right_vy = Gait_data.data(:,9);
% 
% % Detect the Initial Contact (IC) of the gait cycles
% % create logical array: contact is true (1), no contact is false (0)
% 
% threshold = 10;                                 
% logicGRF_left_vy = GRF_left_vy > threshold;
% logicGRF_right_vy = GRF_right_vy > threshold;
% 
% % Determine the indices of initial contact (IC). 
% % Indices for which the difference with the next logical number is +1. 
% 
% indIC_left = find(diff(logicGRF_left_vy)==1);
% indIC_right = find(diff(logicGRF_right_vy)==1);
% 
% TIC_Left = time(indIC_left);
% TIC_Right = time(indIC_right);
% 
% % duration stride left
% 
% StrideDuration_left = diff(TIC_Left);
% StrideDuration_left_mean = mean(StrideDuration_left);
% StrideDuration_left_std = std(StrideDuration_left);
% 
% % Number of standard deviations
% 
% nstd = 2;
% 
% % Define upper and lower limit for bad cycle detection
% 
% Upperlimit_left = StrideDuration_left_mean + nstd * StrideDuration_left_std;
% Lowerlimit_left = StrideDuration_left_mean - nstd * StrideDuration_left_std;
% 
% % Find and remove the indeces of IC that are bad cycles
% % Creating a new list with indeces of good cycles 
% 
% IndGoodGaits_Left = {}; 
% for i = 1:length(StrideDuration_left)
%     if not(StrideDuration_left(i) > Upperlimit_left) || (StrideDuration_left(i) < Lowerlimit_left)
%         if not(i == length(indIC_left))  
%             IndGoodGaits_Left{end+1} = [indIC_left(i), indIC_left(i+1)-1];
%         end
%     end
% end
% 
% % Creating the same list but with time instead of index 
% 
% TimeGoodGaits_Left = cell([1 length(IndGoodGaits_Left)]);
% for i = 1:length(IndGoodGaits_Left)
%     StartIndex_Left = IndGoodGaits_Left{i}(1);
%     EndIndex_Left = IndGoodGaits_Left{i}(2); 
%     TimeGoodGaits_Left{i} = [time(StartIndex_Left), time(EndIndex_Left)];
% end
% 
% % Repeated for right side 
%     
% StrideDuration_right = diff(TIC_Right);
% StrideDuration_right_mean = mean(StrideDuration_right);
% StrideDuration_right_std = std(StrideDuration_right);
% 
% Upperlimit_right = StrideDuration_right_mean + nstd * StrideDuration_right_std;
% Lowerlimit_right = StrideDuration_right_mean - nstd * StrideDuration_right_std;
% 
% IndGoodGaits_Right = {}; 
% for i = 1:length(StrideDuration_right)
%     if not(StrideDuration_right(i) > Upperlimit_right) || (StrideDuration_right(i) < Lowerlimit_right)
%         if not(i == length(indIC_right))
%             IndGoodGaits_Right{end+1} = [indIC_right(i), indIC_right(i+1)-1];
%         end
%     end
% end
% 
% TimeGoodGaits_Right= cell([1 length(IndGoodGaits_Right)]);
% for i = 1:length(IndGoodGaits_Right)
%     StartIndex_Right = IndGoodGaits_Right{i}(1);
%     EndIndex_Right = IndGoodGaits_Right{i}(2); 
%     TimeGoodGaits_Right{i} = [time(StartIndex_Right), time(EndIndex_Right)];
% end
% 


% function [TimeGoodGaits_Left, TimeGoodGaits_Right] = DetectGaitEvents_v2(path)
% % Import data
% 
% Gait_data = importdata(path);

% Extract time + vertical force of left and right leg

time = Gait_data.data(:,1);
GRF_left_vy = Gait_data.data(:,3);
GRF_right_vy = Gait_data.data(:,9);

% Detect the Initial Contact (IC) of the gait cycles
% create logical array: contact is true (1), no contact is false (0)

threshold = 10;                                 
logicGRF_left_vy = GRF_left_vy > threshold;
logicGRF_right_vy = GRF_right_vy > threshold;

% Determine the indices of initial contact (IC). 
% Indices for which the difference with the next logical number is +1. 

indIC_left = find(diff(logicGRF_left_vy)==1);
indIC_right = find(diff(logicGRF_right_vy)==1);

TIC_Left = time(indIC_left);
TIC_Right = time(indIC_right);

% duration stride left

StrideDuration_left = diff(TIC_Left);
StrideDuration_left_mean = mean(StrideDuration_left);
StrideDuration_left_std = std(StrideDuration_left);

% Find and remove the indeces of IC that are bad cycles
% Creating a new list with indeces of good cycles 

IndGoodGaits_Left = {}; 
[~,BadGaits] = rmoutliers(StrideDuration_left, 'quartiles');
for i = 1:length(BadGaits)
    if not(BadGaits(i)) == 1 && not(i == length(BadGaits))
        IndGoodGaits_Left{end+1} = [indIC_left(i), indIC_left(i+1)];
    end
end

% Creating the same list but with time instead of index 

TimeGoodGaits_Left = cell([1 length(IndGoodGaits_Left)]);
for i = 1:length(IndGoodGaits_Left)
    StartIndex_Left = IndGoodGaits_Left{i}(1);
    EndIndex_Left = IndGoodGaits_Left{i}(2); 
    TimeGoodGaits_Left{i} = [time(StartIndex_Left), time(EndIndex_Left)];
end

% Repeated for right side 
    
StrideDuration_right = diff(TIC_Right);
StrideDuration_right_mean = mean(StrideDuration_right);
StrideDuration_right_std = std(StrideDuration_right);

IndGoodGaits_Right = {}; 
[~,BadGaits] = rmoutliers(StrideDuration_right, 'quartiles');
for i = 1:length(BadGaits)
    if not(BadGaits(i)) == 1 && not(i == length(BadGaits))
        IndGoodGaits_Right{end+1} = [indIC_right(i), indIC_right(i+1)];
    end
end

TimeGoodGaits_Right= cell([1 length(IndGoodGaits_Right)]);
for i = 1:length(IndGoodGaits_Right)
    StartIndex_Right = IndGoodGaits_Right{i}(1);
    EndIndex_Right = IndGoodGaits_Right{i}(2); 
    TimeGoodGaits_Right{i} = [time(StartIndex_Right), time(EndIndex_Right)];
end

% output = struct('Time_Right', TimeGoodGaits_Right, ...
%                 'Time_Left', TimeGoodGaits_Left, ...
%                 'Ind_Right', IndGoodGaits_Right, ...
%                 'Ind_Left', IndGoodGaits_Left);







