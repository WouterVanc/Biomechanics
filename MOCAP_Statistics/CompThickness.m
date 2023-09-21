% Compare Thickness OA - CO

clear all;
close all;
clc

set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultAxesFontName', 'Arial')

load ThicknessResultsFinal.mat

thicknesspart = fieldnames(ThicknessResults);
type = fieldnames(ThicknessResults.FWZ);
group = fieldnames(ThicknessResults.FWZ.normal);

% Table info
rownames = thicknesspart;
columnames = {'Part','Mean Thickness OA', 'Std OA', 'Mean Thickness CO', 'Std CO', 'Pval',...
    'Peak Thickness OA', 'Std Peak OA', 'Peak Thickness CO', 'Std Peak CO', 'Pval_peak'};

% Results directories
path_out = fullfile(cd,'Results_Thickness_Comp');
path_out_p1 = fullfile(path_out, 'Figures');

if ~exist("path_out", 'dir')
    mkdir(path_out)
end

if ~exist("path_out_p1", 'dir')
    mkdir(path_out_p1)
end

hh = [];
h2 = [];
hhp = [];
h2p = [];
for q = 1:length(thicknesspart)
        % OA
        mean_oa = mean(ThicknessResults.(thicknesspart{q}).normal.OA);
        std_oa = std(ThicknessResults.(thicknesspart{q}).normal.OA);
        peak_mean_oa = mean(ThicknessResults.(thicknesspart{q}).peak.OA);
        peak_std_oa = std(ThicknessResults.(thicknesspart{q}).peak.OA);
        % CO
        mean_co = mean(ThicknessResults.(thicknesspart{q}).normal.CO);
        std_co = std(ThicknessResults.(thicknesspart{q}).normal.CO);
        peak_mean_co = mean(ThicknessResults.(thicknesspart{q}).peak.CO);
        peak_std_co = std(ThicknessResults.(thicknesspart{q}).peak.CO);
        
        % Normal
        [p,h] = ranksum(ThicknessResults.(thicknesspart{q}).normal.OA, ThicknessResults.(thicknesspart{q}).normal.CO);

        % Peak
        [p2,h9] = ranksum(ThicknessResults.(thicknesspart{q}).peak.OA, ThicknessResults.(thicknesspart{q}).peak.CO);
        
        % Normal
        matrix(q,1) = mean_oa;
        matrix(q,2) = std_oa;
        matrix(q,3) = mean_co;
        matrix(q,4) = std_co;
        matrix(q,5) = p;
        
        % Peak
        matrix(q,6) = peak_mean_oa;
        matrix(q,7) = peak_std_oa;
        matrix(q,8) = peak_mean_co;
        matrix(q,9) = peak_std_co;
        matrix(q,10) = p2;

        % Combine data
        y_temp = ThicknessResults.(thicknesspart{q}).normal.CO;
        y_temp(6:10,1) = NaN;

        y_temp_peak = ThicknessResults.(thicknesspart{q}).peak.CO;
        y_temp_peak(6:10,1) = NaN;
        
        % Normal
        hh = [hh; ThicknessResults.(thicknesspart{q}).normal.OA];
        h2 = [h2; y_temp];

        % Peak
        hhp = [hhp; ThicknessResults.(thicknesspart{q}).peak.OA];
        h2p = [h2p; y_temp_peak];
end
h3 = [hh;h2];
h3p = [hhp; h2p];

group = zeros(100,1);
group(1:50,1) = 2;
group(51:100,1) = 1;

part(1:10,1) = 1; % fullfemur
part(11:20,1) = 2; % FWBZ
part(21:30,1) = 3; % FNWBZ
part(31:40,1) = 4; % LatWBZ
part(41:50,1) = 5; % LatWBZ
part(51:60,1) = 1; % fullfemur
part(61:70,1) = 2; % FWBZ
part(71:80,1) = 3; % FNWBZ
part(81:90,1) = 4; % LatWBZ
part(91:100,1) = 5; % LatWBZ

value = h3;
value2 = h3p;

% Mean
Finalmatrix = [group part value];
colnames = {'group', 'part', 'value'};
table = array2table(Finalmatrix, 'VariableNames',colnames);

% Boxplot
boxchart(table.part,table.value,'GroupByColor',table.group)
title('Mean Thickness')
ylabel('Thickness [mm]')
xticks([1 2 3 4 5])
xticklabels({'FULLFEM','FWBZ', 'FNWBZ', 'LATWBZ', 'MEDWBZ'})
legend('Control', 'OA', Location='southeast')

plotname = (['MeanThicknessBoxplots' '.png']);
saveas(gcf, fullfile(path_out_p1,plotname))

close gcf

% Peak
Finalmatrix = [group part value2];
colnames = {'group', 'part', 'value'};
table = array2table(Finalmatrix, 'VariableNames',colnames);

% Boxplot
boxchart(table.part,table.value,'GroupByColor',table.group)
title('Peak Thickness')
ylabel('Thickness [mm]')
xticks([1 2 3 4 5])
xticklabels({'FULLFEM','FWBZ', 'FNWBZ', 'LATWBZ', 'MEDWBZ'})
legend('Control', 'OA', Location='southeast')

plotname = (['PeakThicknessBoxplots' '.png']);
saveas(gcf, fullfile(path_out_p1,plotname))

close gcf

% Combine data
matrix_total = [rownames num2cell(matrix)];

% Convert to table 
sTable_pval = cell2table(matrix_total, 'VariableNames',columnames);

% Write to Excel file
filename = 'CompThickness.xlsx';

writetable(sTable_pval,fullfile(path_out, filename),"Sheet",'OA_vs_Control')
    


