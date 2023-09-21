clearvars 
close all
clc
%% Directories & parameters 
% subject = 206;
% side = 1; %1=medial 0=lateral
% % enter slices nr of interest 
% SOI = 5:13;
% 
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
    clearvars -except w wouter
    load slices.mat
    subject = wouter(w,1);
    SOI = wouter(w,2):wouter(w,3);

TEfiles = (['E:\Biomedical Research in Movement Sciences\Master Thesis\Masterproef_Data\MRI\Mapping\Control\MRI_S00' num2str(subject) '\ET_files_s00' num2str(subject) '\ET_files_s00' num2str(subject)]);
directory_mask = (['E:\Biomedical Research in Movement Sciences\Master Thesis\Masterproef_Data\MRI\Mapping\FULL_NON_Weightbearing_zones\Masks_NONFWZ\T2\S00' num2str(subject)]);
contentM = dir(directory_mask);


path_mask = fullfile(cd,['MATLAB_masksT2_FNWBZ_s00' num2str(subject)]); 

if~exist(path_mask)
    mkdir(path_mask);
end

if w==5
    SOI = [5:7, 9, 12:20];
end


%echo times 
TE =  [ 11 22 33 44 55 66 77 88 99 110 121 132];
% TE =  [ 11 22 33 44 55 66 77 88 99];

%number of slices 
nr_slices = 26;
nr_TE_total= 12; %don't change
nr_TE = length(TE);
% nr_TE = 12;
%, should be black and white. ==>  2 == draw ROI, 1 = Masked image, 0 = full image (or anyother value)
masking = 1;

%% Read images and mask into 'img' structure 
% All data is read into 1 structure: slices > SLT images + mask
%There should be as many slices as masks in folder, be carefull about name conventions 

image    = dicomread(TEfiles);
info = dicominfo(TEfiles);

for i = 1:length(SOI)
    ET_start=((SOI(i)-1)*nr_TE_total)+1;
    ET_stop = ET_start+ nr_TE-1;
    c=0;
    for file = ET_start: ET_stop
        c = c+1;
        img.(['slice_' num2str(SOI(i))]).(['img_' num2str(c)])  = double(image(:,:,file));

    end

if masking == 1
BinaryImage = imread(fullfile(directory_mask,contentM((SOI(i)+2)).name)); %+2 because first two rows have no content in directory

BinaryImage = BinaryImage(:,:,1)==254; 
mask = flip(BinaryImage, 2); %dim

img.(['slice_' num2str(SOI(i))]).(['Mask_' contentM((SOI(i)+2)).name(2:(end-4))]) = mask;

save(fullfile(path_mask,['slice_' num2str(SOI(i)) '.mat']),'mask');

elseif masking == 2
    Frame = dicomread(TEfiles,"frames",ET_start);
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
img_nr = 2; %for which te time you want to plot images
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
    title(['mask_' num2str(SOI(i))]);
end 

figure
slice = SOI(4); %for which slices you want to see TE times
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

[r,c] = find(img.(['slice_' num2str(SOI(j))]).(['Mask_' contentM((SOI(j)+2)).name(2:(end-4))])==1);

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
path_out = fullfile(cd,['Processed_FNWBZ_s00' num2str(subject)]); 

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

listRT=[];
for i = 1:length(SOI)
    Matrix = Filtered_MapT2.(['slice_' num2str(SOI(i))]);
    listRT = [listRT ConvertMatrix(Matrix)];
end

MeanRT = mean(listRT);
stdRT = std(listRT);

%% write filtered results 

path_out_f = fullfile(cd,['Filtered_FNWBZ_s00' num2str(subject)]); 

if~exist(path_out_f)
    mkdir(path_out_f);
end

for j = 1:length(SOI)
    slicemapf = Filtered_MapT2.(['slice_' num2str(SOI(j))]);
    sliceresf = Filtered_NormRes_T2.(['slice_' num2str(SOI(j))]);
    save(fullfile(path_out_f,['Filtered_MapT2_slice_' num2str(SOI(j)) '.mat']),'slicemapf'); 
    save(fullfile(path_out_f,['Filtered_Normresidual_slice_' num2str(SOI(j)) '.mat']),'sliceresf'); 
end

save(fullfile(cd,['T2_FNWBZ_S00' num2str(subject)]), 'MeanRT','stdRT')

disp('%%%%%%%%%%%%%%%%%%%%%%%')
disp(num2str(subject))
disp(num2str(MeanRT))
disp(num2str(stdRT))
disp('%%%%%%%%%%%%%%%%%%%%%%%')

close all;

end


