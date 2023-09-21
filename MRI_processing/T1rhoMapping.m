clear all
close all
clc
%% Directories & parameters 
datapath = 'E:\Biomedical Research in Movement Sciences\Master Thesis\Masterproef_Data\MRI\Mapping\OA\MRI_S201\DICOMT1rho';
% directory_mask = 'C:\Users\woute\Desktop\s001\s001\T1rho_mask_femur_cart\Femur';
SLfiles = "E:\Biomedical Research in Movement Sciences\Master Thesis\Masterproef_Data\MRI\Mapping\OA\MRI_S201\SL_file\SL_10_s201";
% contentM = dir(directory_mask);
content = dir(datapath);
nfiles = size(content,1);

path_mask = fullfile(cd,'MATLAB_masksT1rho'); 

if~exist(path_mask)
    mkdir(path_mask);
end

% enter slices nr of interest 
% SOI = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
SOI=[10];

%Spin lock times
SLT =  [ 0 10 20 30 40 50 60];

%, should be black and white. ==>  2 == draw ROI, 1 = Masked image, 0 = full image (or anyother value)
masking = 0;

%% Read images and mask into 'img' structure 
% All data is read into 1 structure: slices > SLT images + mask
%There should be as many slices as masks in folder, be carefull about name conventions 

for i = 1:length(SOI)
    c=0;
for file = 3:nfiles 
    if ~content(file).isdir
        
        file_in  = fullfile(datapath,content(file).name);
        c = c+1;
        temp    = dicomread(file_in);
        info = dicominfo(file_in);
        img.(['slice_' num2str(SOI(i))]).(['img_' num2str(c)])  = double(temp(:,:,SOI(i)));
        clear temp
    end
end
if masking == 1
BinaryImage = imread(fullfile(directory_mask,contentM((SOI(i)+2)).name)); %+2 because first two rows have no content in directory

BinaryImage = BinaryImage(:,:,1)==254; 
mask = flip(BinaryImage, 2); %dim

img.(['slice_' num2str(SOI(i))]).mask = mask;

elseif masking == 2
    Frame = dicomread(SLfiles,"frames",SOI(i));
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

% %% Example plot figure 
% img_nr = 10; %for which SLT time you want to plot images
% for i = 1:length(SOI);
%     figure;
%     
%     subplot(1,2,1);
%     imagesc(img.(['slice_' num2str(SOI(i))]).(['img_' num2str(img_nr)])); %3
%     colormap gray;
%     axis off ;
%     title(['slice ' num2str(SOI(i))]);
% 
%     subplot(1,2,2);
%     imagesc(img.(['slice_' num2str(SOI(i))]).mask); %double()
%     axis off ;
%     title(['mask_' num2str(SOI(i))]);
% end 
% 
% figure
% slice = 10; %for which slices you want to see SLT times
% for i = 1:length(SLT);
%     
%     subplot(3,3,i)
%     imagesc(img.(['slice_' num2str(slice)]).(['img_' num2str(i)]));
%     colormap gray
%     if i == 1
%         clim = get(gca,'Clim');
%     else
%         set(gca,'CLim',clim);
%     end
%     axis off
%     title(['SLT ' num2str(SLT(i))]);
% end

%% Fitting algorithm
for j = 1:length(SOI)
for i = 1:length(SLT)
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
        [x5, resnorm] = lsqcurvefit(fun5,x0_5, unique(SLT)', data, lb, ub, options2);
        MapT1rho.(['slice_' num2str(SOI(j))])(r(m),c(m),1) = x5(2);
        Normresidual.(['slice_' num2str(SOI(j))])(r(m),c(m),1) = norm(((data-fun5(x0_5,double(unique(SLT))))))./norm((data));
        
    %display with threshold
    threshold = 200;
    MapT1rho.(['slice_' num2str(SOI(j))])(MapT1rho.(['slice_' num2str(SOI(j))]) > threshold) = threshold;
    MapT1rho.(['slice_' num2str(SOI(j))])(MapT1rho.(['slice_' num2str(SOI(j))]) <0) = 0;
        
end 
  %% Result figures      
        figure
        subplot(121)
        imagesc(MapT1rho.(['slice_' num2str(SOI(j))])(:,:,1))
        colorbar
        colormap jet
        set(gca,'CLim',[0 200])
        title('T1rho map')
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
path_out = fullfile(cd,'Processed'); 

if~exist(path_out)
    mkdir(path_out);
end

for j = 1:length(SOI);
    slicemap = MapT1rho.(['slice_' num2str(SOI(j))]);
    sliceres = Normresidual.(['slice_' num2str(SOI(j))]);
    save(fullfile(path_out,['MapT1rho_slice_' num2str(SOI(j)) '.mat']),'slicemap'); 
    save(fullfile(path_out,['Normresidual_slice_' num2str(SOI(j)) '.mat']),'sliceres'); 
end 

%% Filter Slices

% T1rho upper threshold = 130, T2 upper threshold
Upper_Threshold = 130;
Lower_Threshold = 0;

Filtered_MapT1rho = struct;
Filtered_NormRes_T1rho = struct;

for i = 1:length(SOI)

    Logic_MapT1rho = MapT1rho.(['slice_' num2str(SOI(i))]) > Lower_Threshold & MapT1rho.(['slice_' num2str(SOI(i))]) < Upper_Threshold;
    % figure; imshow(Logic_MapT1rho)

    Filtered_MapT1rho.(['slice_' num2str(SOI(i))]) = Logic_MapT1rho.*MapT1rho.(['slice_' num2str(SOI(i))]);
    Filtered_NormRes_T1rho.(['slice_' num2str(SOI(i))]) = Logic_MapT1rho.*Normresidual.(['slice_' num2str(SOI(i))]);
    % figure; imshow(int8(Filtered_MapT1rho))

end

%% Average and Standard Deviation MapT1rho

cutoff = 9; % last lateral slice

Lat_list = [];
for i = 1:cutoff
    Matrix = Filtered_MapT1rho.(['slice_' num2str(SOI(i))]);
    Lat_list = [Lat_list ConvertMatrix(Matrix)];
end

Med_list = [];
for i = cutoff+1:length(SOI)
    Matrix = Filtered_MapT1rho.(['slice_' num2str(SOI(i))]);
    Med_list = [Med_list ConvertMatrix(Matrix)];
end

Mean_subject = mean([Lat_list Med_list]);
Std_subject = std([Lat_list Med_list]);

Mean_subject_lat = mean(Lat_list);
Mean_subject_med = mean(Med_list);

Std_subject_lat = std(Lat_list);
Std_subject_med = std(Med_list);

%% %% Average and Standard Deviation Normresidual

cutoff = 9; % last lateral slice

Lat_list_norm = [];
for i = 1:cutoff
    Matrix = Filtered_NormRes_T1rho.(['slice_' num2str(SOI(i))]);
    Lat_list_norm = [Lat_list_norm ConvertMatrix(Matrix)];
end

Med_list_norm = [];
for i = cutoff+1:length(SOI)
    Matrix = Filtered_NormRes_T1rho.(['slice_' num2str(SOI(i))]);
    Med_list_norm = [Med_list_norm ConvertMatrix(Matrix)];
end

Mean_subject_norm = mean([Lat_list_norm Med_list_norm]);
Std_subject_norm = std([Lat_list_norm Med_list_norm]);

Mean_subject_lat_norm = mean(Lat_list_norm);
Mean_subject_med_norm = mean(Med_list_norm);

Std_subject_lat_norm = std(Lat_list_norm);
Std_subject_med_norm = std(Med_list_norm);

%% write filtered results 

path_out_f = fullfile(cd,'Filtered'); 

if~exist(path_out_f)
    mkdir(path_out_f);
end

for j = 1:length(SOI);
    slicemapf = Filtered_MapT1rho.(['slice_' num2str(SOI(j))]);
    sliceresf = Filtered_NormRes_T1rho.(['slice_' num2str(SOI(j))]);
    save(fullfile(path_out_f,['Filtered_MapT1rho_slice_' num2str(SOI(j)) '.mat']),'slicemapf'); 
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
