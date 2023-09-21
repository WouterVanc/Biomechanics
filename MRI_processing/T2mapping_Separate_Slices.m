clear all
close all
clc
%% Directories & parameters for T2 images
clear all
close all
clc
 
% directory_mask = 'C:\Users\u0139912\Documents\MATLAB\Mapping_Wouter\s001\T1rho_mask'; %for T2 masks from mimics
TEfiles = 'C:\Users\woute\Desktop\Thesis\TeDoen\MRI_203_20072022\ET_files';
% contentM = dir(directory_mask);

path_mask = fullfile(cd,'MATLAB_masksT2_tibia_s203'); 

if~exist(path_mask)
    mkdir(path_mask);
end

% enter slices nr of interest 
SOI=[4 5 6 7 8 9 10 11 12 13 18 19 20 21 22 23 24 25];
%echo times 
TE =  [ 11 22 33 44 55 66 77 88 99 110 121 132];

%number of slices 
nr_slices = 30; %check this number 
nr_TE_total= 12; %don't change
nr_TE = length(TE);
%, should be black and white. ==>  2 == draw ROI, 1 = Masked image, 0 = full image (or anyother value)
masking = 2;

% Read images and mask into 'img' structure 
% All data is read into 1 structure: slices > SLT images + mask
%There should be as many slices as masks in folder, be carefull about name conventions 

%list of the seperate images/slices
files = dir(fullfile(TEfiles)); 
files=files(~ismember({files.name},{'.','..'}));

for i = 1:length(SOI)
    ET_start=((SOI(i)-1)*nr_TE_total)+1;
    ET_stop = ET_start+ nr_TE-1;
    c=0;
    for file = ET_start: ET_stop
        image    = dicomread(fullfile(TEfiles,files(file).name));
        c = c+1;
        img.(['slice_' num2str(SOI(i))]).(['img_' num2str(c)])  = double(image); %(:,:,file)

    end

if masking == 1
BinaryImage = imread(fullfile(directory_mask,contentM((SOI(i)+2)).name)); %+2 because first two rows have no content in directory

BinaryImage = BinaryImage(:,:,1)==254; 
mask = flip(BinaryImage, 2); %dim

img.(['slice_' num2str(SOI(i))]).mask = mask;

elseif masking == 2
    Frame = dicomread(fullfile(TEfiles,files(ET_start).name));%dicomread(TEfiles,"frames",ET_start);
    imshow(Frame, {}, 'Border', 'tight')
    hold on
    set(gcf, 'Position', get(0, 'Screensize'));
    axis tight
    roi = images.roi.AssistedFreehand('Smoothing', 10);
    draw(roi);
    hold off
    mask = createMask(roi);
    img.(['slice_' num2str(SOI(i))]).mask = mask;
    close all
    
    save(fullfile(path_mask,['slice_' num2str(SOI(i)) '.mat']),'mask'); 
    
else
    mask = ones(320);
    img.(['slice_' num2str(SOI(i))]).mask = mask;
end

end
%% Example plot figure 
img_nr = 5; %for which te time you want to plot images
for i = 1:length(SOI);
    figure;
    
    subplot(1,2,1);
    imagesc(img.(['slice_' num2str(SOI(i))]).(['img_' num2str(img_nr)])); %3
    colormap gray;
    axis off ;
    title(['slice ' num2str(SOI(i))]);

    subplot(1,2,2);
    imagesc(img.(['slice_' num2str(SOI(i))]).mask); %double()
    axis off ;
    title(['mask_' num2str(SOI(i))]);
end 

figure
slice = 8; %for which slices you want to see TE times
for i = 1:length(TE);
    
    subplot(3,4,i)
    imagesc(img.(['slice_' num2str(slice)]).(['img_' num2str(i)]));
    colormap gray
    if i == 1
        clim = get(gca,'Clim');
    else
        set(gca,'CLim',clim);
    end
    axis off
    title(['SLT ' num2str(TE(i))]);
end

%% Fitting algorithm
for j = 1:length(SOI)
for i = 1:length(TE)
    Data.(['slice_' num2str(SOI(j))])(:,:,i)  = img.(['slice_' num2str(SOI(j))]).(['img_' num2str(i)]);
end


options2 = optimset('display','off','MaxFunEvals',600,'Maxiter',400,'TolFun',1e-8,'TolX',1e-8,'Algorithm','levenberg-marquardt');
dimensions = size(Data.(['slice_' num2str(SOI(j))]));

[r,c] = find(img.(['slice_' num2str(SOI(j))]).mask==1);

for m = 1:length(r)
        data = (reshape(Data.(['slice_' num2str(SOI(j))])(r(m),c(m),:),dimensions(3),1));
        fun5 = @(x,xdata)[x(1)*exp(-xdata/x(2))];
        x0_5 = [ max(data) 200]; %200 van philips, heeft niet veel invloed 
        lb = [];
        ub =[];
        [x5, resnorm] = lsqcurvefit(fun5,x0_5, unique(TE)', data, lb, ub, options2);
        MapT2.(['slice_' num2str(SOI(j))])(r(m),c(m),1) = x5(2);
        Normresidual.(['slice_' num2str(SOI(j))])(r(m),c(m),1) = norm(((data-fun5(x0_5,double(unique(TE))))))./norm((data));
        
    %display with threshold
    threshold = 200;
    MapT2.(['slice_' num2str(SOI(j))])(MapT2.(['slice_' num2str(SOI(j))]) > threshold) = threshold;
    MapT2.(['slice_' num2str(SOI(j))])(MapT2.(['slice_' num2str(SOI(j))]) <0) = 0;
        
end 
  %% Result figures      
        figure
        subplot(121)
        imagesc(MapT2.(['slice_' num2str(SOI(j))])(:,:,1))
        colorbar
        colormap jet
        set(gca,'CLim',[0 200])
        title('MapT2 map')
        axis off

        subplot(122)
        imagesc(Normresidual.(['slice_' num2str(SOI(j))])(:,:,1))
        colorbar
        colormap jet
        set(gca,'CLim',[0 5])
        title('Normalised error map')
        axis off
        
        display(['ended slice ' num2str(SOI(j))])
end 
%% write results

% [a,b] = fileparts(datapath);
path_out = fullfile(cd,'T2_Processed_tibia_s203'); 

if~exist(path_out)
    mkdir(path_out);
end

for j = 1:length(SOI);
    slicemap = MapT2.(['slice_' num2str(SOI(j))]);
    sliceres = Normresidual.(['slice_' num2str(SOI(j))]);
    save(fullfile(path_out,['MapT2_slice_' num2str(SOI(j)) '.mat']),'slicemap'); 
    save(fullfile(path_out,['Normresidual_slice_' num2str(SOI(j)) '.mat']),'sliceres'); 
end 

%% Filter Slices

% T1rho upper threshold = 130, T2 upper threshold
Upper_Threshold = 100;
Lower_Threshold = 0;

Filtered_MapT2 = struct;
Filtered_NormRes_T2 = struct;

for i = 1:length(SOI)

    Logic_MapT2 = MapT2.(['slice_' num2str(SOI(i))]) > Lower_Threshold & MapT2.(['slice_' num2str(SOI(i))]) < Upper_Threshold;
    % figure; imshow(Logic_MapT1rho)

    Filtered_MapT2.(['slice_' num2str(SOI(i))]) = Logic_MapT2.*MapT2.(['slice_' num2str(SOI(i))]);
    Filtered_NormRes_T2.(['slice_' num2str(SOI(i))]) = Logic_MapT2.*Normresidual.(['slice_' num2str(SOI(i))]);
    % figure; imshow(int8(Filtered_MapT1rho))

end

%% Average and Standard Deviation MapT1rho

cutoff = 14; 

Lat_list = [];
for i = 1:cutoff
    Matrix = Filtered_MapT2.(['slice_' num2str(SOI(i))]);
    Lat_list = [Lat_list ConvertMatrix(Matrix)];
end

Med_list = [];
for i = cutoff+1:length(SOI)
    Matrix = Filtered_MapT2.(['slice_' num2str(SOI(i))]);
    Med_list = [Med_list ConvertMatrix(Matrix)];
end

Mean_subject = mean([Lat_list Med_list]);
Std_subject = std([Lat_list Med_list]);

Mean_subject_lat = mean(Lat_list);
Mean_subject_med = mean(Med_list);

Std_subject_lat = std(Lat_list);
Std_subject_med = std(Med_list);

%% %% Average and Standard Deviation Normresidual

Lat_list_norm = [];
for i = 1:cutoff
    Matrix = Filtered_NormRes_T2.(['slice_' num2str(SOI(i))]);
    Lat_list_norm = [Lat_list_norm ConvertMatrix(Matrix)];
end

Med_list_norm = [];
for i = cutoff+1:length(SOI)
    Matrix = Filtered_NormRes_T2.(['slice_' num2str(SOI(i))]);
    Med_list_norm = [Med_list_norm ConvertMatrix(Matrix)];
end

Mean_subject_norm = mean([Lat_list_norm Med_list_norm]);
Std_subject_norm = std([Lat_list_norm Med_list_norm]);

Mean_subject_lat_norm = mean(Lat_list_norm);
Mean_subject_med_norm = mean(Med_list_norm);

Std_subject_lat_norm = std(Lat_list_norm);
Std_subject_med_norm = std(Med_list_norm);

%% write filtered results 

path_out_f = fullfile(cd,'T2_Filtered_tibia_s203'); 

if~exist(path_out_f)
    mkdir(path_out_f);
end

for j = 1:length(SOI);
    slicemapf = Filtered_MapT2.(['slice_' num2str(SOI(j))]);
    sliceresf = Filtered_NormRes_T2.(['slice_' num2str(SOI(j))]);
    save(fullfile(path_out_f,['Filtered_MapT2_slice_' num2str(SOI(j)) '.mat']),'slicemapf'); 
    save(fullfile(path_out_f,['Filtered_Normresidual_slice_' num2str(SOI(j)) '.mat']),'sliceresf'); 
end

save(fullfile(path_out_f, 'Mean_subject'));
save(fullfile(path_out_f, 'Mean_subject_lat'));
save(fullfile(path_out_f, 'Mean_subject_med'));
save(fullfile(path_out_f, 'Mean_subject_lat_norm'));
save(fullfile(path_out_f, 'Mean_subject_med_norm'));
save(fullfile(path_out_f, 'Mean_subject_norm'));

save(fullfile(path_out_f, 'Std_subject'));
save(fullfile(path_out_f, 'Std_subject_lat'));
save(fullfile(path_out_f, 'Std_subject_med'));
save(fullfile(path_out_f, 'Std_subject_lat_norm'));
save(fullfile(path_out_f, 'Std_subject_med_norm'));
save(fullfile(path_out_f, 'Std_subject_norm'));

