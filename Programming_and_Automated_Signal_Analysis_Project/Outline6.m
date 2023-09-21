%% 6

% Load the data

load("DataCodeAssignment2122/rawGRF/GRF_filtered.mat");
load("DataCodeAssignment2122/rawGRF/GRF_raw.mat");

% Analyzing the frequency spectrum using the function 'make_freq_spect' from Maarten Afschrift

[data_ft,frequency_out,magnitude_out] = make_freq_spect(F1,FrameRate,1);
xlim([-1 20]);
[data_ft,frequency_out,magnitude_out] = make_freq_spect(F1_filt,FrameRate,1);
xlim([-1 20]);

% Measuring the variability 

V_GRFraw = var(F1);
V_GRFfilt = var(F1_filt);
V_diff = V_GRFraw - V_GRFfilt;

% Filtering the raw data with the right cut off score and filter 

fc = 10; 
fs = FrameRate; 
n = 4; 
Wn = fc/(fs/2); 
Y_raw = F1;

[b,a] = butter(n,Wn,'low');

Y_filt = filtfilt(b,a,Y_raw);

figure(1)
plot(Y_filt)
figure(2)
make_freq_spect(Y_filt,FrameRate,1)
xlim([-1 20]);

% Visualizing the influence of different cut off scores on the variability 

n = 4;
Y_raw = F1;
VariabilityDependance = zeros(200,3);

for i= 1:200
    
    [b,a] = butter(n,((0.1*i)/(FrameRate/2)),'low');

    Y_filt = filtfilt(b,a,Y_raw);
    
    VariabilityDependance(i, 1:3) = var(Y_filt);

end

fig = figure();
subplot(1,3,1)
plot(linspace(1,20,200), VariabilityDependance(:,1))
title('GRFx')
subplot(1,3,2)
plot(linspace(1,20,200), VariabilityDependance(:,2))
title('GRFy')
subplot(1,3,3)
plot(linspace(1,20,200), VariabilityDependance(:,3))
title('GRFz')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Variability');
xlabel(han,'Cut-off frequency');

str0 = 'figures/';
str1 = 'Variability_Cut-off_frequency';
figuretitel = (append(str0,str1));
savefig(fig,figuretitel,"compact")

close(fig)


