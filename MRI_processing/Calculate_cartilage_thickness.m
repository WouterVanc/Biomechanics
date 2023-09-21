%% Calculate thickness maps
%Wouter Van Caekenberghe 10/2022
%Ikram Mohout 08/2022

clearvars; close all; clc

%% Define input = path to folder stl files 
tic
path.stl = uigetdir('c:\', 'Select folder that contains STL files');
content = dir(path.stl);
nfiles = size(content,1);

%make structure with al STL files & triangulations 

for file = 3:nfiles  
    if ~content(file).isdir
        
       file_Cartilage  = fullfile(path.stl,content(file).name);
       STLfiles.(['stl_' num2str(file-2)]).TR  = stlread(file_Cartilage);
       STLfiles.(['stl_' num2str(file-2)]).F = faceNormal(STLfiles.(['stl_' num2str(file-2)]).TR);
       STLfiles.(['stl_' num2str(file-2)]).P = incenter(STLfiles.(['stl_' num2str(file-2)]).TR);
    end
end 
%% Plot the loaded stl files 
for i= 1:length(fieldnames(STLfiles))
    %plot STL file 
    figure;
    trisurf(STLfiles.(['stl_' num2str(i)]).TR.ConnectivityList,STLfiles.(['stl_' num2str(i)]).TR.Points(:,1),STLfiles.(['stl_' num2str(i)]).TR.Points(:,2),STLfiles.(['stl_' num2str(i)]).TR.Points(:,3));
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    box off
    grid off
    
    %plot face normal 
    figure;
    quiver3(STLfiles.(['stl_' num2str(i)]).P(:,1),STLfiles.(['stl_' num2str(i)]).P(:,2),STLfiles.(['stl_' num2str(i)]).P(:,3), STLfiles.(['stl_' num2str(i)]).F(:,1),STLfiles.(['stl_' num2str(i)]).F(:,2),STLfiles.(['stl_' num2str(i)]).F(:,3),0.5,'color','r');
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    box off
    grid off
end  
    
%% Thickness calculation algorithm 
%folder to write results to
path_out = fullfile(cd,'STL_raw_Processed'); 

if~exist(path_out)
    mkdir(path_out);
end
%minimal angle between two normals 
MinAngle = 140;

for x= 1:length(fieldnames(STLfiles))
    STLfiles.(['stl_' num2str(x)]).Result.index = [];
    STLfiles.(['stl_' num2str(x)]).Result.angle = [];
    STLfiles.(['stl_' num2str(x)]).Result.distance = [];
for j= 1:length(STLfiles.(['stl_' num2str(x)]).F)
    F1 = STLfiles.(['stl_' num2str(x)]).F(j,:); P1 = STLfiles.(['stl_' num2str(x)]).P(j,:);
    Points.index = [];
    Points.angle = [];
    Points.distance = [];
for i= 1:length(STLfiles.(['stl_' num2str(x)]).F)
    P2 = STLfiles.(['stl_' num2str(x)]).P(i,:); F2 = STLfiles.(['stl_' num2str(x)]).F(i,:);
    dt = dot(F1,F2); CosTheta = dt/(norm(F1)*norm(F2));
    Theta = acosd(CosTheta);
    if MinAngle <= Theta %& pdist2(P1,P2,'euclidean') < MaxDist
        Points.index = [Points.index; i];
        Points.angle = [Points.angle; Theta];
        Points.distance = [Points.distance pdist2(P1,P2,'euclidean')]; 
    end 
end 
if ~isempty(Points.index)
MinDist = min(Points.distance);
IND = find(Points.distance == min(MinDist));
STLfiles.(['stl_' num2str(x)]).Result.index = [STLfiles.(['stl_' num2str(x)]).Result.index; j, IND];
STLfiles.(['stl_' num2str(x)]).Result.angle = [STLfiles.(['stl_' num2str(x)]).Result.angle; Points.angle(IND,1)];
STLfiles.(['stl_' num2str(x)]).Result.distance = [STLfiles.(['stl_' num2str(x)]).Result.distance; MinDist];
end 
end 
%write save raw results 
 Structure = STLfiles.(['stl_' num2str(x)]);
 save(fullfile(path_out,[content((x+2)).name  '.mat']),'Structure'); 
end

str = 'Finished raw thickness' ; 
disp(str)

%%
% Filter out take out bad edges
MaxDist= 7; 

%folder to write results to
path_out_filtered = fullfile(cd,'STL_filtered_Processed'); 

if~exist(path_out_filtered)
    mkdir(path_out_filtered);
end


for x= 1:length(fieldnames(STLfiles))
     load(fullfile(path_out,[content((x+2)).name  '.mat']),'Structure');
     c = 1;
while c < length(Structure.Result.distance)
    if Structure.Result.distance(c,:) > MaxDist
        Structure.Result.distance(c,:) = [];
        Structure.Result.index(c,:) = [];
        Structure.Result.angle(c,:) = [];
    else  
    c = c+1;
    end 
end 
save(fullfile(path_out_filtered,[content((x+2)).name  '.mat']),'Structure');
clear Structure
end

str2 = 'Filtering Parts' ; 
disp(str2)

%% Plots of results
    
for x= 1:length(fieldnames(STLfiles))
    
    load(fullfile(path_out_filtered,[content((x+2)).name  '.mat']),'Structure');
    
    %plot cartilage thickness  
    
    cVect =Structure.Result.distance;
    figure;
%   hh = trisurf(STLfiles.(['stl_' num2str(x)]).TR.ConnectivityList([STLfiles.(['stl_' num2str(x)]).Result.index(:,1)],:),STLfiles.(['stl_' num2str(i)]).TR.Points(:,1),STLfiles.(['stl_' num2str(i)]).TR.Points(:,2),STLfiles.(['stl_' num2str(i)]).TR.Points(:,3));
    hh = trisurf(Structure.TR.ConnectivityList([Structure.Result.index(:,1)],:),Structure.TR.Points(:,1),Structure.TR.Points(:,2),Structure.TR.Points(:,3));

    set(gca,'CLim',[min(cVect), max(cVect)]);
    set(hh,'FaceColor','flat',...
       'FaceVertexCData',cVect,...
       'CDataMapping','scaled');
    colorbar;

 xlabel('x [mm]')
 ylabel('y [mm]')
 zlabel('z [mm]')
 box off
 grid off;
hold on; 
end 

%% Average Thickness and STD per subject
path_out_results = fullfile(cd,'STL_MeanStd'); 

if~exist(path_out_results, 'dir')
    mkdir(path_out_results);
end

for x= 1:length(fieldnames(STLfiles))
    load(fullfile(path_out_filtered,[content((x+2)).name  '.mat']),'Structure');
    MeanDistance = mean(Structure.Result.distance);
    StdDistance = std(Structure.Result.distance);

    save(fullfile(path_out_results,[content((x+2)).name  '.mat']),'MeanDistance', 'StdDistance');
end

%% Computing Duration

t = toc;

t1 = floor(t/86400); % days

t1_temp = rem(t,86400); 

t2 = floor(t1_temp/3600); % hours

t2_temp = mod(t1_temp,3600); 

t3 = floor(t2_temp/60); % minutes

t4 = mod(t2_temp,60); % seconds

disp (['|||||||||| Finished processing in ' num2str(t1) ' Days ' num2str(t2) ' Hours ' num2str(t3) ' Minutes and ' num2str(t4) ' seconds ||||||||||' ]);

