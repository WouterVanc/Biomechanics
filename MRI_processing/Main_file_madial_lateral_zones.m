%%% DEVIDE FEMUR IN MEDIAL & LATERAL ZONES
% Wouter Van Caekenberghe adapted 05/12/2022
clearvars; close all; clc

%Ikram 1/12/2022
% Define input
path.stl = 'D:\WORK_USB\Master Thesis\Masterproef_Data\MRI\Thickness\Control\STL_files\Femur';
MaxDist= 5;

name_cartilage = 's001_Smoothed_Wrapped_Femoral_Cartilage_final.stl';
file_Cartilage = fullfile(path.stl,name_cartilage);
path_out = fullfile(cd, 'Thickness_sides_results');

if ~exist(path_out, 'dir')
    mkdir(path_out)
end

[TR] = stlread(file_Cartilage);
trisurf(TR.ConnectivityList,TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
 xlabel('x [mm]')
 ylabel('y [mm]')
 zlabel('z [mm]')
 box off
 grid off
pause
%select the intercondylar notch
%press enter and then select the point of interest
[X,Y,Z]=ginput(1)
pause
[X2,Y2,Z2]=ginput(1)
v_notch = [X,Y,Z];
v_top = [X2,Y2,Z2];
V = TR.Points;
F = TR.ConnectivityList;

prompt = {'1 is rechts-mediaal, 0 is links-mediaal', 'Subject nummer'};
answer = inputdlg(prompt);
Whatside = str2double(answer(1,1));
Subnr = str2double(answer(2,1));

%% Find medial & lateral side 
% Rechterkant
SideA = V(:,1) <= ((V(:,1)-v_top(1))*(v_notch(2)-v_top(2))-(V(:,2)-v_top(2))*(v_notch(1)-v_top(1))); 
% Linkerkant
SideB = V(:,1) >= ((V(:,1)-v_top(1))*(v_notch(2)-v_top(2))-(V(:,2)-v_top(2))*(v_notch(1)-v_top(1))); 

vertices_sideA = V(SideA,:);
vertices_sideB = V(SideB,:);

Face_sideA = SideA(F); 
if isequal(size(Face_sideA), [3 1]);
    Face_sideA = Face_sideA.';            
end
Face_sideA = any(Face_sideA,2);  

Face_sideB = SideB(F); 
if isequal(size(Face_sideB), [3 1]);
    Face_sideB = Face_sideB.';            
end
Face_sideB = any(Face_sideB,2); 
%% Figures Ikram 
figure
subplot(121);
h = trisurf(TR.ConnectivityList,TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
hold on 
ht1 = trisurf(TR.ConnectivityList(Face_sideA,:),TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
xlabel('x');ylabel('y');zlabel('z')
set(h,'EdgeColor','k')
set(ht1,'EdgeColor','r')
grid off

subplot(122);
h = trisurf(TR.ConnectivityList,TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
hold on   
ht1 = trisurf(TR.ConnectivityList(Face_sideB,:),TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
xlabel('x');ylabel('y');zlabel('z')
set(h,'EdgeColor','k')
set(ht1,'EdgeColor','b')
grid off

%% Keep thickness for only medial lateral zones
load ("D:\WORK_USB\Master Thesis\Masterproef_Data\MRI\Thickness\Control\STL_raw_Processed\s001_Smoothed_Wrapped_Femoral_Cartilage_final.stl.mat")

%keep indexes only for resolts on medial/lateral side 

ThicknessSideA.Result.index = []; ThicknessSideA.Result.angle = []; ThicknessSideA.Result.distance = [];
ThicknessSideB.Result.index = []; ThicknessSideB.Result.angle = []; ThicknessSideB.Result.distance = [];

% Rechterkant
for x=1:length(Face_sideA);
    if Face_sideA(x);  
        ThicknessSideA.Result.index = [ThicknessSideA.Result.index; Structure.Result.index(x,:)];
        ThicknessSideA.Result.angle = [ThicknessSideA.Result.angle; Structure.Result.angle(x,:)];
        ThicknessSideA.Result.distance = [ThicknessSideA.Result.distance; Structure.Result.distance(x,:)];
    end 
end 

% Filter out take out bad edges
 
c = 1;
while c < length(ThicknessSideA.Result.distance)
    if ThicknessSideA.Result.distance(c,:) > MaxDist
        ThicknessSideA.Result.distance(c,:) = [];
        ThicknessSideA.Result.index(c,:) = [];
        ThicknessSideA.Result.angle(c,:) = [];
    else  
    c = c+1;
    end 
end 

% Linkerkant
for x=1:length(Face_sideB);
    if Face_sideB(x);  
        ThicknessSideB.Result.index = [ThicknessSideB.Result.index; Structure.Result.index(x,:)];
        ThicknessSideB.Result.angle = [ThicknessSideB.Result.angle; Structure.Result.angle(x,:)];
        ThicknessSideB.Result.distance = [ThicknessSideB.Result.distance; Structure.Result.distance(x,:)];
    end 
end 

% Filter out take out bad edges
 
c = 1;
while c < length(ThicknessSideB.Result.distance)
    if ThicknessSideB.Result.distance(c,:) > MaxDist
        ThicknessSideB.Result.distance(c,:) = [];
        ThicknessSideB.Result.index(c,:) = [];
        ThicknessSideB.Result.angle(c,:) = [];
    else  
    c = c+1;
    end 
end 

if Whatside == 1
    Medial = ThicknessSideA.Result;
    Lateral = ThicknessSideB.Result;
elseif Whatside == 0
    Medial = ThicknessSideB.Result;
    Lateral = ThicknessSideA.Result;
end

% save(fullfile(path_out,['SideResults_' num2str(Subnr)]),"Medial","Lateral")
%% figures 
cVect =ThicknessSideA.Result.distance;
    figure;
    subplot(121);
    hh = trisurf(Structure.TR.ConnectivityList([ThicknessSideA.Result.index(:,1)],:),Structure.TR.Points(:,1),Structure.TR.Points(:,2),Structure.TR.Points(:,3));

    set(gca,'CLim',[min(cVect), max(cVect)]);
    set(hh,'FaceColor','flat',...
       'FaceVertexCData',cVect,...
       'CDataMapping','scaled');
    colorbar;
    grid off
    
    subplot(122);
cVect2 =ThicknessSideB.Result.distance;
    hh2 = trisurf(Structure.TR.ConnectivityList([ThicknessSideB.Result.index(:,1)],:),Structure.TR.Points(:,1),Structure.TR.Points(:,2),Structure.TR.Points(:,3));

    set(gca,'CLim',[min(cVect2), max(cVect2)]);
    set(hh2,'FaceColor','flat',...
       'FaceVertexCData',cVect2,...
       'CDataMapping','scaled');
    colorbar;
    grid off

% Plot full for comparison

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

cVect = Structure.Result.distance;
    figure;
    hh = trisurf(Structure.TR.ConnectivityList([Structure.Result.index(:,1)],:),Structure.TR.Points(:,1),Structure.TR.Points(:,2),Structure.TR.Points(:,3));
    
    set(gca,'CLim',[min(cVect), max(cVect)]);
    set(hh,'FaceColor','flat',...
       'FaceVertexCData',cVect,...
       'CDataMapping','scaled');
    colormap jet
    
%     colorbar;
    grid off
    axis off

pause

close all;

