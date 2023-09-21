%%% DEVIDE FEMUR IN WEIGHTBEARING ZONES / FIND PEAK AND MEAN THICKNESS
% Wouter Van Caekenberghe adapted 30/12/2022

clearvars; 
close all; 
clc

%Ikram adapted 18/11/2022
% Define input

MaxDist= 7; % Maximum thickness value to remove bad edges
path.stl = 'D:\WORK_USB\Master Thesis\Masterproef_Data\MRI\Thickness\Control\STL_files\Femur';
name_cartilage = 's001_Smoothed_Wrapped_Femoral_Cartilage_final.stl';
file_Cartilage = fullfile(path.stl,name_cartilage);
stlname = ['STL_WeightBearing_fem_s' name_cartilage(2:4)];

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
[X,Y,Z]=ginput(1);
v_notch = [X,Y,Z];
V = TR.Points;
F = TR.ConnectivityList;

prompt = {'1 is rechts-mediaal, 0 is links-mediaal'};
answer = inputdlg(prompt);
Whatside = str2double(answer(1,1));

%% Find most posterior point
%Original code Sam

[value_post,pos_post] = max(V(:,2)); 

distance_yz = sqrt((v_notch(2)-V(pos_post,2))^2+(v_notch(3)-V(pos_post,3))^2);
percentage = 0.60*distance_yz; 

final_point = v_notch(2) + percentage; 

under_y = V(:,2) >= v_notch(2); 
above_y = V(:,2) <=final_point; 
combined = under_y == above_y; 
vertices_under1 = V(combined,:);

average_z = nanmean(vertices_under1(:,3));
dev_z = nanstd(vertices_under1(:,3),0,1);


under_z =  vertices_under1(:,3) <= V(pos_post,3);
vertices_under = vertices_under1(under_z,:);


vertices = ismember(V,vertices_under,'rows');
face_sel = vertices(F);                      
if isequal(size(face_sel), [3 1])
    face_sel = face_sel.';            
end
face_sel = any(face_sel,2);  

vertices2 = ismember(V,vertices_under1,'rows');
face_sel2 = vertices2(F);                      
if isequal(size(face_sel2), [3 1])
    face_sel2 = face_sel2.';            
end
face_sel2 = any(face_sel2,2);  
%% Figures adapted 18/11 Ikram 
figure
% subplot(121);
h = trisurf(TR.ConnectivityList,TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
hold on 

ht1 = trisurf(TR.ConnectivityList(face_sel,:),TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
xlabel('x');ylabel('y');zlabel('z')
set(h,'EdgeColor','k')
set(ht1,'EdgeColor','r')
grid off
axis off

% subplot(122)
% h = trisurf(TR.ConnectivityList,TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
% hold on 
% 
% ht1 = trisurf(TR.ConnectivityList(face_sel2,:),TR.Points(:,1),TR.Points(:,2),TR.Points(:,3));
% xlabel('x');ylabel('y');zlabel('z')
% set(h,'EdgeColor','k')
% set(ht1,'EdgeColor','b')
%% Save thickness and STL, full weight-bearing zone

path_out_results = fullfile(cd,'Thickness_results_WBZ'); 

if~exist(path_out_results, 'dir')
    mkdir(path_out_results);
end

% Find thickness
load (['D:\WORK_USB\Master Thesis\Masterproef_Data\MRI\Thickness\OA\STL_raw_Processed\s' name_cartilage(2:4) '_Smoothed_Wrapped_Femoral_Cartilage_final.stl.mat'])

ThicknessFWZ.Result.index = []; ThicknessFWZ.Result.angle = []; ThicknessFWZ.Result.distance = [];

for x=1:length(face_sel)
    if face_sel(x)  
        ThicknessFWZ.Result.index = [ThicknessFWZ.Result.index; Structure.Result.index(x,:)];
        ThicknessFWZ.Result.angle = [ThicknessFWZ.Result.angle; Structure.Result.angle(x,:)];
        ThicknessFWZ.Result.distance = [ThicknessFWZ.Result.distance; Structure.Result.distance(x,:)];
    end 
end 

% Filter: take out bad edges
c = 1;
while c < length(ThicknessFWZ.Result.distance)
    if ThicknessFWZ.Result.distance(c,:) > MaxDist
        ThicknessFWZ.Result.distance(c,:) = [];
        ThicknessFWZ.Result.index(c,:) = [];
        ThicknessFWZ.Result.angle(c,:) = [];
    else  
    c = c+1;
    end 
end 

% Find mean and peak thickness
MeanThicknessFWZ = mean(ThicknessFWZ.Result.distance);
StdThicknessFWZ = std(ThicknessFWZ.Result.distance);

nvalues = round(length(ThicknessFWZ.Result.distance)*0.1);
topten = sort(ThicknessFWZ.Result.distance, 'descend');
PeakThicknessFWZ = mean(topten(1:nvalues));
StdPeakThicknessFWZ = std(topten(1:nvalues));

% Create new STL
ThicknessFWZ.TR1 = triangulation(Structure.TR.ConnectivityList([ThicknessFWZ.Result.index(:,1)],:),Structure.TR.Points);
figure;

h = trisurf(ThicknessFWZ.TR1.ConnectivityList,ThicknessFWZ.TR1.Points(:,1),ThicknessFWZ.TR1.Points(:,2),ThicknessFWZ.TR1.Points(:,3));


cVect2 =ThicknessFWZ.Result.distance;
    set(gca,'CLim',[min(cVect2), max(cVect2)]);
    set(h,'FaceColor','flat',...
       'FaceVertexCData',cVect2,...
       'CDataMapping','scaled');
    colorbar;
hold on 
 
stlwrite(ThicknessFWZ.TR1,[stlname '.stl'])

pause;

close all;

%% Medial and Lateral compartment 

trisurf(ThicknessFWZ.TR1.ConnectivityList,ThicknessFWZ.TR1.Points(:,1),ThicknessFWZ.TR1.Points(:,2),ThicknessFWZ.TR1.Points(:,3));
 xlabel('x [mm]')
 ylabel('y [mm]')
 zlabel('z [mm]')
 box off
 grid on
pause
%press enter and then select the points of interest
[X,Y,Z]=ginput(1);
pause
[X2,Y2,Z2]=ginput(1);
v_notch = [X,Y,Z];
v_top = [X2,Y2,Z2];
V = ThicknessFWZ.TR1.Points;
F = ThicknessFWZ.TR1.ConnectivityList;

close all;

%% Find medial & lateral side 

% Rechterkant
SideA = V(:,1) <= ((V(:,1)-v_top(1))*(v_notch(2)-v_top(2))-(V(:,2)-v_top(2))*(v_notch(1)-v_top(1)));
% Linkerkant
SideB = V(:,1) >= ((V(:,1)-v_top(1))*(v_notch(2)-v_top(2))-(V(:,2)-v_top(2))*(v_notch(1)-v_top(1))); 

Face_sideA = SideA(F); 
if isequal(size(Face_sideA), [3 1])
    Face_sideA = Face_sideA.';            
end
Face_sideA = any(Face_sideA,2);  

Face_sideB = SideB(F); 
if isequal(size(Face_sideB), [3 1])
    Face_sideB = Face_sideB.';            
end
Face_sideB = any(Face_sideB,2);

%% Figures Ikram 

figure
subplot(121);
h10 = trisurf(ThicknessFWZ.TR1.ConnectivityList,ThicknessFWZ.TR1.Points(:,1),ThicknessFWZ.TR1.Points(:,2),ThicknessFWZ.TR1.Points(:,3));
hold on 
ht1 = trisurf(ThicknessFWZ.TR1.ConnectivityList(Face_sideA,:),ThicknessFWZ.TR1.Points(:,1),ThicknessFWZ.TR1.Points(:,2),ThicknessFWZ.TR1.Points(:,3));
xlabel('x');ylabel('y');zlabel('z')
set(h10,'EdgeColor','k')
set(ht1,'EdgeColor','r')
grid off

subplot(122);
h11 = trisurf(ThicknessFWZ.TR1.ConnectivityList,ThicknessFWZ.TR1.Points(:,1),ThicknessFWZ.TR1.Points(:,2),ThicknessFWZ.TR1.Points(:,3));
hold on 
ht1 = trisurf(ThicknessFWZ.TR1.ConnectivityList(Face_sideB,:),ThicknessFWZ.TR1.Points(:,1),ThicknessFWZ.TR1.Points(:,2),ThicknessFWZ.TR1.Points(:,3));
xlabel('x');ylabel('y');zlabel('z')
set(h11,'EdgeColor','k')
set(ht1,'EdgeColor','b')
grid off

%% Keep thickness for medial & lateral zones

ThicknessSideA.Result.index = ThicknessFWZ.Result.index(Face_sideA,:);
ThicknessSideA.Result.angle = ThicknessFWZ.Result.angle(Face_sideA,:);
ThicknessSideA.Result.distance = ThicknessFWZ.Result.distance(Face_sideA,:);
       
% Filter: take out bad edges
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

ThicknessSideB.Result.index = ThicknessFWZ.Result.index(Face_sideB,:);
ThicknessSideB.Result.angle = ThicknessFWZ.Result.angle(Face_sideB,:);
ThicknessSideB.Result.distance = ThicknessFWZ.Result.distance(Face_sideB,:);

% Filter: take out bad edges
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

% Plot thickness map FWZ to compare

    cVect111 =ThicknessFWZ.Result.distance;
    figure;
    hh = trisurf(ThicknessFWZ.TR1.ConnectivityList,ThicknessFWZ.TR1.Points(:,1),ThicknessFWZ.TR1.Points(:,2),ThicknessFWZ.TR1.Points(:,3));
    set(gca,'CLim',[min(cVect111), max(cVect111)]);
    set(hh,'FaceColor','flat',...
       'FaceVertexCData',cVect111,...
       'CDataMapping','scaled');
    colorbar;

if Whatside == 1
    ThicknessMWZ.Results = ThicknessSideA.Result;
    ThicknessLWZ.Results = ThicknessSideB.Result;
    % Create new STL medial side
    ThicknessMWZ.TR1 = triangulation(ThicknessFWZ.TR1.ConnectivityList(Face_sideA,:),ThicknessFWZ.TR1.Points);
    figure
    h40 = trisurf(ThicknessMWZ.TR1.ConnectivityList,ThicknessMWZ.TR1.Points(:,1),ThicknessMWZ.TR1.Points(:,2),ThicknessMWZ.TR1.Points(:,3));
    cVect2 =ThicknessMWZ.Results.distance;
    set(gca,'CLim',[min(cVect2), max(cVect2)]);
    set(h40,'FaceColor','flat',...
       'FaceVertexCData',cVect2,...
       'CDataMapping','scaled');
    colorbar;
    hold on
    % Create new STL lateral side
    ThicknessLWZ.TR1 = triangulation(ThicknessFWZ.TR1.ConnectivityList(Face_sideB,:),ThicknessFWZ.TR1.Points);
    figure
    h30 = trisurf(ThicknessLWZ.TR1.ConnectivityList,ThicknessLWZ.TR1.Points(:,1),ThicknessLWZ.TR1.Points(:,2),ThicknessLWZ.TR1.Points(:,3));
    cVect23 =ThicknessLWZ.Results.distance;
    set(gca,'CLim',[min(cVect23), max(cVect23)]);
    set(h30,'FaceColor','flat',...
       'FaceVertexCData',cVect23,...
       'CDataMapping','scaled');
    colorbar;
    hold on
elseif Whatside == 0
    ThicknessMWZ.Results = ThicknessSideB.Result;
    ThicknessLWZ.Results = ThicknessSideA.Result;
    % Create new STL medial side
    ThicknessMWZ.TR1 = triangulation(ThicknessFWZ.TR1.ConnectivityList(Face_sideB,:),TR.Points);
    figure
    h220 = trisurf(ThicknessMWZ.TR1.ConnectivityList,ThicknessMWZ.TR1.Points(:,1),ThicknessMWZ.TR1.Points(:,2),ThicknessMWZ.TR1.Points(:,3));
    cVect22 =ThicknessMWZ.Results.distance;
    set(gca,'CLim',[min(cVect22), max(cVect22)]);
    set(h220,'FaceColor','flat',...
       'FaceVertexCData',cVect22,...
       'CDataMapping','scaled');
    colorbar;
    hold on
    % Create new STL lateral side
    ThicknessLWZ.TR1 = triangulation(ThicknessFWZ.TR1.ConnectivityList(Face_sideA,:),TR.Points);
    figure
    h20 = trisurf(ThicknessLWZ.TR1.ConnectivityList,ThicknessLWZ.TR1.Points(:,1),ThicknessLWZ.TR1.Points(:,2),ThicknessLWZ.TR1.Points(:,3));
    cVect2 =ThicknessLWZ.Results.distance;
    set(gca,'CLim',[min(cVect2), max(cVect2)]);
    set(h20,'FaceColor','flat',...
       'FaceVertexCData',cVect2,...
       'CDataMapping','scaled');
    colorbar;
    hold on
end

% Find mean and peak thickness lateral side
MeanThicknessLWZ = mean(ThicknessLWZ.Results.distance);
StdThicknessLWZ = std(ThicknessLWZ.Results.distance);

nvalues = round(length(ThicknessLWZ.Results.distance)*0.1);
topten = sort(ThicknessLWZ.Results.distance, 'descend');
PeakThicknessLWZ = mean(topten(1:nvalues));
StdPeakThicknessLWZ = std(topten(1:nvalues));

% Find mean and peak thickness medial side 
MeanThicknessMWZ = mean(ThicknessMWZ.Results.distance);
StdThicknessMWZ = std(ThicknessMWZ.Results.distance);

nvalues = round(length(ThicknessMWZ.Results.distance)*0.1);
topten = sort(ThicknessMWZ.Results.distance, 'descend');
PeakThicknessMWZ = mean(topten(1:nvalues));
StdPeakThicknessMWZ = std(topten(1:nvalues));

% Save STL medial and lateral side, and variables
% stlwrite(ThicknessLWZ.TR1,[stlname '_lateral.stl'])
% stlwrite(ThicknessMWZ.TR1,[stlname '_medial.stl'])

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Full weightbearing zone thickness = ' num2str(MeanThicknessFWZ)])
disp(['Medial weightbearing zone thickness = ' num2str(MeanThicknessMWZ)])
disp(['Lateral weightbearing zone thickness = ' num2str(MeanThicknessLWZ)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Save(fullfile(path_out_results,[name_cartilage(1:end-10) '.mat']),'MeanThicknessFWZ', 'StdThicknessFWZ', "PeakThicknessFWZ",...
%     "StdPeakThicknessFWZ", 'MeanThicknessMWZ','StdThicknessMWZ', "PeakThicknessMWZ", "StdPeakThicknessMWZ", 'MeanThicknessLWZ',...
%     'StdThicknessLWZ', "PeakThicknessLWZ", "StdPeakThicknessLWZ",...
%     "ThicknessFWZ","ThicknessLWZ", "ThicknessMWZ");

% press enter to close figures
pause

close all;
