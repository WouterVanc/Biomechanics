clearvars 
close all
clc
%% Directories & parameters 
% subject = 210;
% 
% % enter slices nr of interest 
% SOI = [8,9,10,11,12,13];

% side = 1; %1=medial 0=lateral
% if side==1
%     lol1 = 'Medial';
%     lol2 = 'medial';
%     lol3 = 'MED';
% elseif side==0
%     lol1 = 'Lateral';
%     lol2 = 'lateral';
%     lol3 = 'LAT';
% end

load slices.mat
for w = 1:length(wouter)
    clearvars -except wouter w
    subject = wouter(w,1);
    SOI = wouter(w,2):wouter(w,3);
datapath = (['E:\Biomedical Research in Movement Sciences\Master Thesis\Masterproef_Data\MRI\Mapping\OA\MRI_S' num2str(subject) '\DICOMT1rho']);
directory_mask = (['E:\Biomedical Research in Movement Sciences\Master Thesis\Masterproef_Data\MRI\Mapping\FULL_NON_Weightbearing_zones\Masks_NONFWZ\T1\s' num2str(subject)]);
contentM = dir(directory_mask);
content = dir(datapath);

path_mask = fullfile(cd,['MATLAB_masksT1rho_FNWZ__s' num2str(subject)]); 

if~exist(path_mask)
    mkdir(path_mask);
end

% Spin lock times
if w==3
    SLT =  [ 0 10 20 30 40 50];
else
    SLT =  [ 0 10 20 30 40 50 60];
end

%, should be black and white. ==>  1 = Masked image, 0 = full image (or anyother value)
masking = 1;

folders = dir(fullfile(datapath)); 
folders=folders(~ismember({folders.name},{'.','..'}));
nfiles = size(folders,1);

%% Read images and mask into 'img' structure 
% All data is read into 1 structure: slices > SLT images + mask
%There should be as many slices as masks in folder, be carefull about name conventions 

for i = 1:length(SOI)
    c=0;
for file = 1:nfiles 
        files = dir(fullfile(datapath,folders(file).name)); 
        files=files(~ismember({files.name},{'.','..'}));
    
        file_in  = fullfile(files(SOI(i)).folder,files(SOI(i)).name);
        c = c+1;
        temp    = dicomread(file_in);
        info = dicominfo(file_in);
        img.(['slice_' num2str(SOI(i))]).(['img_' num2str(c)])  = double(temp); %(:,:,SOI(i))
        clear temp
end
if masking == 1
BinaryImage = imread(fullfile(directory_mask,contentM((SOI(i)+2)).name)); %+2 because first two rows have no content in directory

BinaryImage = BinaryImage(:,:,1)==254; 
mask = flip(BinaryImage, 2); %dim

img.(['slice_' num2str(SOI(i))]).(['Mask_' contentM((SOI(i)+2)).name(2:(end-4))]) = mask;

save(fullfile(path_mask,['slice_' num2str(SOI(i)) '.mat']),'mask'); 
else
    mask = ones(320);
    img.(['slice_' num2str(SOI(i))]).(['Mask_' contentM((SOI(i)+2)).name(2:(end-4))]) = mask;
end
end 


%% Example plot figure 
img_nr = 2; %for which SLT time you want to plot images
for i = 1:length(SOI);
    figure;
    
    subplot(1,2,1);
    imagesc(img.(['slice_' num2str(SOI(i))]).(['img_' num2str(img_nr)])); %3
    colormap gray;
    axis off ;
    title(['slice ' num2str(SOI(i))]);

    subplot(1,2,2);
    imagesc(img.(['slice_' num2str(SOI(i))]).(['Mask_' contentM((SOI(i)+2)).name(2:(end-4))])); %double()
    axis off ;
    title(['mask ' (contentM((SOI(i)+2)).name(2:(end-4)))]);
end 

% figure
% slice = 5; %for which slices you want to see SLT times
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

[r,c] = find(img.(['slice_' num2str(SOI(j))]).(['Mask_' contentM((SOI(j)+2)).name(2:(end-4))])==1);
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
path_out = fullfile(cd,['Processed_FNWZ__s' num2str(subject)]); 

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

listRT=[];
for i = 1:length(SOI)
    Matrix = Filtered_MapT1rho.(['slice_' num2str(SOI(i))]);
    listRT = [listRT ConvertMatrix(Matrix)];
end

MeanRT = mean(listRT);
stdRT = std(listRT);

%% write filtered results 

path_out_f = fullfile(cd,['Filtered_FNWZ__s' num2str(subject)]); 

if~exist(path_out_f)
    mkdir(path_out_f);
end

for j = 1:length(SOI)
    slicemapf = Filtered_MapT1rho.(['slice_' num2str(SOI(j))]);
    sliceresf = Filtered_NormRes_T1rho.(['slice_' num2str(SOI(j))]);
    save(fullfile(path_out_f,['Filtered_MapT1rho_slice_' num2str(SOI(j)) '.mat']),'slicemapf'); 
    save(fullfile(path_out_f,['Filtered_Normresidual_slice_' num2str(SOI(j)) '.mat']),'sliceresf'); 
end

save(fullfile(cd,['T1RHO_FNWZ_S' num2str(subject)]), 'MeanRT','stdRT')

disp('%%%%%%%%%%%%%%%%%%%%%%%')
disp(num2str(subject))
disp(num2str(MeanRT))
disp(num2str(stdRT))
disp('%%%%%%%%%%%%%%%%%%%%%%%')
close all;
end

 




