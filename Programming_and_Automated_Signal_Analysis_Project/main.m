 clear; close all;

cd('/Users/woutervancaekenberghe/Documents/KU Leuven/Master BReMS/Programming & automated signal analysis/')
tic
%% 1

Filtered_IC = struct;
age = {'young','older'};

% Iterating every subject through the DetectGaitEvents function
% Returning the output of the function in a struct

for i = 1:2
    folder = strcat('DataCodeAssignment2122/osimData_',age{i},'/');
    D = dir(folder);
    for k = 1:length(D)
        if not(D(k).name(1)=='.') % skip hidden folders
            file_loc = strcat(folder,D(k).name,'/RefWalk/GRF/Refwalk_GRF.mot');
            [IndGoodGaits_Left, IndGoodGaits_Right] = DetectGaitEvents(file_loc);
            Filtered_IC.(age{i}).(D(k).name).IndLeft = IndGoodGaits_Left;
            Filtered_IC.(age{i}).(D(k).name).IndRight = IndGoodGaits_Right;
        end
    end
end

%% 2

data = struct;

% You can only choose one side at a time to create the graphs in section 3
% + fill in the side you have chosen

kinetics = {'ground_force_vx', 'ground_force_vy', 'ground_force_vz'};
kinematics = {'hip_flexion_l', 'hip_rotation_l', 'hip_adduction_l', 'knee_angle_l', 'ankle_angle_l'};
ChosenSide = 'Left';
Total = [kinetics kinematics];
                        
% Iterating every goodgait through the GetData function to return the
% interpolated data along with the means and stdev's

for i = 1:2
    folder = strcat('DataCodeAssignment2122/osimData_',age{i},'/');
    D = dir(folder);
    for k = 1:length(D)
        if not(D(k).name(1)=='.') % skip hidden folders
            goodgaits = Filtered_IC.(age{i}).(D(k).name);
            data.(age{i}).(D(k).name) = GetData( ...
                age{i}, D(k).name, goodgaits, ...
                kinetics, kinematics);
        end
    end
end

%% 3

% Iterating all the data through the plotData function to
% return all the necessary plots in a folder. 

for i = 1:2
    folder = strcat('DataCodeAssignment2122/osimData_',age{i},'/');
    D = dir(folder);
    for k = 1:length(D)
        if not(D(k).name(1)=='.')
           plotData(age{i},D(k).name, kinetics, kinematics, data, ChosenSide);
        end
    end
end

%% 4

n_sub = [17,9];
matrix = zeros(101,length(Total));
CompareGroups = struct;

% Calculating mean stdev by adding every subject into an empty matrix and
% dividing by the number of subjects of the group

for i=1:length(age)
    folder = strcat('DataCodeAssignment2122/osimData_',age{i},'/');
    D = dir(folder);
    for k = 1:length(D)
        if not(D(k).name(1)=='.')
           for t = 1:length(Total)
                  matrix(:,t) = matrix(:,t) + (data.(age{i}).(D(k).name).(ChosenSide).stdev{1,t});
                  [Max, IndMax] = max(matrix/n_sub(i));
                  [Min, IndMin] = min(matrix/n_sub(i));
            end
        end
    end
    CompareGroups.(age{i}).MeanStdev = matrix/n_sub(i);
    CompareGroups.(age{i}).IndexMax = IndMax;
    CompareGroups.(age{i}).IndexMin = IndMin;
    CompareGroups.(age{i}).MaxValues = Max;
    CompareGroups.(age{i}).MinValues = Min;
end

% Iterating the previous calculated data that is put in a struct through a modified plotData
% function to return the necessary plots in a folder

for i = 1:2
    plotDatav2(age{i}, kinetics, kinematics,CompareGroups, ChosenSide)
end

% Writing the data to an excel file 

Extrema = {'IndexMax','IndexMin','MaxValues','MinValues'};
variable_names_cellarray = [Total Extrema];
variable_names_table = cell2table(variable_names_cellarray);

cellyoung = struct2cell(CompareGroups.young);
MeanStdev_young = cell2mat(cellyoung);

cellolder = struct2cell(CompareGroups.older);
MeanStdev_older = cell2mat(cellolder);

s1 = 'MeanStdev_Young_Old';
s2 = (ChosenSide);
s3 = '.xls';
filename = append(s1,'_',s2,s3);

writetable(variable_names_table, filename, 'Sheet', 'MeanStdev_young', 'WriteVariableNames', false, 'Range', 'A1')
writematrix(MeanStdev_young(1:101,:), filename,'Sheet','MeanStdev_young','Range','A2');
writematrix(MeanStdev_young(102:end,:).', filename,'Sheet','MeanStdev_young','Range','I2');

writetable(variable_names_table, filename, 'Sheet', 'MeanStdev_older', 'WriteVariableNames', false, 'Range', 'A1')
writematrix(MeanStdev_older(1:101,:), filename,'Sheet','MeanStdev_older','Range','A2');
writematrix(MeanStdev_older(102:end,:).', filename,'Sheet','MeanStdev_older','Range','I2');

%% 5

% Calculate Mean stance duration and std stance duration for all subjects

data_strides = struct; 
data_bar_mean = zeros(1,(n_sub(1)));
data_bar_std = zeros(1,(n_sub(1)));
Data_bar_stance = zeros(70,17);

for w = 1:2
    folder = strcat('DataCodeAssignment2122/osimData_',age{w},'/');
    D = dir(folder);
    count = 0;
    for k = 1:length(D)
        if not(D(k).name(1)=='.')
            count = count + 1;
            data_strides = data.(age{w}).(D(k).name).(ChosenSide).data(:,1);   % put the data in the struct
        
            data.(age{w}).(D(k).name).(ChosenSide).Logical = logical(cell2mat((data_strides.')));        % make a logic to see where there is a GRF_vx (1) and where its value is 0
            data.(age{w}).(D(k).name).(ChosenSide).EndStance = diff(logical(cell2mat((data_strides.'))));        % calculate where we go from 1 to 0, this will give us -1
        
            EndStance = diff(logical(cell2mat((data_strides.'))));
            size_EndStance = size(EndStance);   % find the number of strides
            for z = 1:size_EndStance(1,2)
                StanceDurationStride = find(EndStance(:,z)==-1);    % find the stance duration for every stride (= index of -1)
                Data_bar_stance(:,z) = StanceDurationStride;
                data.(age{w}).(D(k).name).(ChosenSide).Stance(1,z) = StanceDurationStride;
                data.(age{w}).(D(k).name).(ChosenSide).MeanStance = mean(data.(age{w}).(D(k).name).(ChosenSide).Stance,'all'); % calculate mean stance duration
                data.(age{w}).(D(k).name).(ChosenSide).StdStance = std(data.(age{w}).(D(k).name).(ChosenSide).Stance);  % calculate std stance duration
            end
        
            % Compute the data that we will put in the bar graph
         if (D(k).name(1)=='S')
            data_bar_young_mean(1,count) = data.young.(D(k).name).(ChosenSide).MeanStance;  % Mean stance of every subject
            data_bar_young_std(1,count) = data.young.(D(k).name).(ChosenSide).StdStance;    % Std Stance of every subject
         else
            data_bar_older_mean(1,count) = data.older.(D(k).name).(ChosenSide).MeanStance;  % Mean stance of every subject
            data_bar_older_std(1,count) = data.older.(D(k).name).(ChosenSide).StdStance;    % Std Stance of every subject
         end
        end
    end
end

% Make the bar graph

f1 = figure();
hold on;
bar(data_bar_young_mean,'FaceColor',[.5 .5 .5]); 
xlabel('Subjects')
ylabel('Stance Phase Duration [%]')
title('Relative Stance Phase Duration for young subjects')

count = 0;
folder = strcat('DataCodeAssignment2122/osimData_young');
D = dir(folder);
for k = 1:length(D)
        if not(D(k).name(1)=='.')
            count = count + 1;
                C(1,count) = {D(k).name};
        end
end

xticks(1:length(C))
xticklabels(C)

err = data_bar_young_std;
er = errorbar(1:n_sub(1),data_bar_young_mean,err,err);
er.Color = [0 0 0]; 
er.LineStyle = 'none';

count = 0;
folder = strcat('DataCodeAssignment2122/osimData_young');
D = dir(folder);
for k = 1:length(D)
        if not(D(k).name(1)=='.')
            count = count + 1;
plot(count,data.young.(D(k).name).(ChosenSide).Stance, '.');
        end
end
hold off;

str0 = 'figures/';
str1 = 'MeanStanceDuration';
str2 = 'young';
str3 = (ChosenSide);
figuretitel1 = (append(str0,str1,'_',str2,'_',str3));
savefig(f1,figuretitel1,"compact")

close(f1)

f2 = figure();
hold on;
bar(data_bar_older_mean,'FaceColor',[.5 .5 .5]); 
xlabel('Subjects')
ylabel('Stance Phase Duration [%]')
title('Relative Stance Phase Duration for older subjects')

err = data_bar_older_std;
er = errorbar(1:n_sub(2),data_bar_older_mean,err,err);
er.Color = [0 0 0]; 
er.LineStyle = 'none';

count = 0;
folder = strcat('DataCodeAssignment2122/osimData_older');
D = dir(folder);
for k = 1:length(D)
        if not(D(k).name(1)=='.')
            count = count + 1;
                W(1,count) = {D(k).name};
        end
end

xticks(1:length(W))
xticklabels(W)
set(gca,'TickLabelInterpreter','none')

count = 0;
folder = strcat('DataCodeAssignment2122/osimData_older');
D = dir(folder);
for k = 1:length(D)
        if not(D(k).name(1)=='.')
            count = count + 1;
plot(count,data.older.(D(k).name).(ChosenSide).Stance, '.');
        end
end
hold off;

str0 = 'figures/';
str1 = 'MeanStanceDuration';
str2 = 'older';
str3 = (ChosenSide);
figuretitel2 = (append(str0,str1,'_',str2,'_',str3));
savefig(f2,figuretitel2,"compact")

close(f2)

%% Outputs

filename1 = 'FilteredIC.mat';
save(filename1,"Filtered_IC");

filename2 = 'InterpolatedData.mat';
save(filename2,"data");

filename3 = 'CompareGroups.mat';
save(filename3, "CompareGroups");

toc
