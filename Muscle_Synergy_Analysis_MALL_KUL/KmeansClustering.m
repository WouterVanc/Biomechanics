%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Kmeans Clustering of Muscle Synergy Analysis Output %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wouter Van Caekenberghe, November 2022

clearvars;
close all;
clc

%% INPUT VARIABLES

k = 5; % How many clusters do you want? (=max amount of synergies found in previous step)
subnumbs = [203 204 205 207 208 209 210]; % numbers('names') of the subjects you want to cluster
OA_or_CO = 1; % 1 for OsteoArthritis, 0 for control 
IsDataNormalized = 1; % 1 for normalized data, 0 for non normalized data
WhatSide = 1; % 1 for Right side, 0 for Left side

%% INITIALIZED VARIABLES

EMG_clust = [];
EMG_weights = [];
sides = char({'Right', 'Left'});
type = char({'OA', 'CO'});
str = char({'_Normalized' '_NOT_Normalized'});

%% DIRECTORIES

if OA_or_CO==1
    if IsDataNormalized==1
        if WhatSide==1
            path_out = fullfile(cd, ['KmeansClustering_Results_' type(1,:) '_' sides(1,:) '_' str(1,1:11)]);
        elseif WhatSide==0
            path_out = fullfile(cd, ['KmeansClustering_Results_' type(1,:) '_' sides(2,1:4) '_' str(1,1:11)]);
        end   
    elseif IsDataNormalized==0
         if WhatSide==1
            path_out = fullfile(cd, ['KmeansClustering_Results_' type(1,:) '_' sides(1,:) '_' str(2,:)]);
        elseif WhatSide==0
            path_out = fullfile(cd, ['KmeansClustering_Results_' type(1,:) '_' sides(2,1:4) '_' str(2,:)]);
        end 
    end
elseif OA_or_CO==0
    if IsDataNormalized==1
        if WhatSide==1
            path_out = fullfile(cd, ['KmeansClustering_Results_' type(2,:) '_' sides(1,:) '_' str(1,1:11)]);
        elseif WhatSide==0
            path_out = fullfile(cd, ['KmeansClustering_Results_' type(2,:) '_' sides(2,1:4) '_' str(1,1:11)]);
        end   
    elseif IsDataNormalized==0
         if WhatSide==1
            path_out = fullfile(cd, ['KmeansClustering_Results_' type(2,:) '_' sides(1,:) '_' str(2,:)]);
        elseif WhatSide==0
            path_out = fullfile(cd, ['KmeansClustering_Results_' type(2,:) '_' sides(2,1:4) '_' str(2,:)]);
        end 
    end
end
if~exist(path_out, 'dir')
    mkdir(path_out)
end
    
%% Perform Kmeans Clustering

% loading the data
if WhatSide==1
    for i = 1:length(subnumbs)
        data_NMF.(['Subject_' num2str(subnumbs(i))]) = load (fullfile(cd, ['\Output_Matrix_NMF_Right_S' num2str(subnumbs(i)) '\synergy_nmf.mat']));
        data_Synergy.(['Subject_' num2str(subnumbs(i))]) = load (fullfile(cd, ['\Output_Matrix_NMF_Right_S' num2str(subnumbs(i)) '\Synergy_Data_R.mat']));
    end
elseif WhatSide==0
     for i = 1:length(subnumbs)
        data_NMF.(['Subject_' num2str(subnumbs(i))]) = load (fullfile(cd, ['\Output_Matrix_NMF_Left_S' num2str(subnumbs(i)) '\synergy_nmf.mat']));
        data_Synergy.(['Subject_' num2str(subnumbs(i))]) = load (fullfile(cd, ['\Output_Matrix_NMF_Left_S' num2str(subnumbs(i)) '\Synergy_Data_L.mat']));
     end
end

% Creating EMG_cluster: all k synergies per subject in 1 matrix >C
for i = 1:length(subnumbs)
    Matrix_temp = data_NMF.(['Subject_' num2str(subnumbs(i))]).synergy_nmf(k).C;  
    EMG_clust = [EMG_clust; Matrix_temp];
end

% Kmeans clustering
[idx, C, sumd, D] = kmeans(EMG_clust,k, 'Replicates', 100);

% Count the amount of synergies in a cluster
for i = 1:k
    ClusterFreq(i,1) = sum(idx==i);
end
%% Make Bar Plot with weights of synergies 

% Sorting data per muscle 
FML = data_Synergy.(['Subject_' num2str(subnumbs(1))]).MuscleNames; % Final Muscle List for distribution of the rest (= Muscle list first subject)

for i = 1:length(subnumbs)
    Matrix_temp = (data_NMF.(['Subject_' num2str(subnumbs(i))]).synergy_nmf(k).W)';
    SML = data_Synergy.(['Subject_' num2str(subnumbs(i))]).MuscleNames;
    [w,v] = size(Matrix_temp);
    for q = 1:w
        for y = 1:length(FML)
            h = find(FML(y)==SML);
            Matrix_temp(q,y) = Matrix_temp(q,h);
            clear h
        end
    end
    EMG_weights = [EMG_weights; Matrix_temp];
end

% Averaging weights of synergies of the same cluster
Combined_Weights = zeros(k,length(FML));
temp = zeros(1,length(FML));

for i = 1:k
    Combined_Weights_temp = EMG_weights(idx==i,:);
    if size(Combined_Weights_temp,1)>1
        temp = sum(Combined_Weights_temp);
        Combined_Weights(i,:) = temp./ClusterFreq(i);
    else
        Combined_Weights(i,:) = Combined_Weights_temp;
    end
end

% Plot the weights and frequency in bar graphs
cat = categorical(FML);
catorg = reordercats(cat,FML);

if OA_or_CO==1
    group = type(1,:);
elseif OA_or_CO==0
    group = type(2,:);
end

if WhatSide==1
    w=1;
    z=1:5;
elseif WhatSide==0
    w=2;
    z =1:4;
end

t = tiledlayout('flow', 'TileSpacing','loose');
for i = 1:k
    nexttile
    bar(catorg,Combined_Weights(i,1:length(Combined_Weights)), 0.4);
    ax = gca;
    ax.XTickLabelRotation = 90;
    ax.XAxis.FontSize = 7;
    title(['Synergy ' num2str(i)], FontWeight="normal")
    title(t,['Muscle Synergy Analysis ' group ' - ' sides(w,z)] , FontWeight="bold")
    ylabel ('Average Weight', FontSize=8)
    ylim([0 1])
end
f = gcf;
f.WindowState = 'maximized';
saveas(f,fullfile(path_out,'Barplot_KmeansClustering'));
close gcf
hold off

% nexttile('east')
% BarGraph_frequency = bar(1:k,ClusterFreq(1:k), 0.4);
% title('Synergy Frequency')
% xlabel('Synergy', FontSize=11)
% ylabel('Frequency', FontSize=11)

BarGraph_frequency = bar(1:k,ClusterFreq(1:k), 0.4);
title(['Synergy Frequency ' group ' - '  sides(w,z)])
xlabel('Synergy', FontSize=11)
ylabel('Frequency', FontSize=11)
f = gcf;
f.WindowState = 'maximized';
saveas(f,fullfile(path_out,'Barplot_Frequency'));
close gcf

%% Save data

save(fullfile(path_out, 'KmeansClusteringResults_Data'), 'idx', 'C', 'sumd', 'D', 'BarGraph_frequency', "t");

