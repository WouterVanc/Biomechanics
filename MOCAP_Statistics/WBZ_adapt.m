% Weightbearing vs non-Weightbearing % 

clear all 
close all
clc

set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultAxesFontName', 'Arial')

load ThicknessResultsFinal.mat
load RTresults.mat

% Fieldnames
thicknesspart = {'FWZ', 'FNWZ'};
relaxpart = {'FWBZ', 'FNWBZ'};

% Results directories
path_out = fullfile(cd,'Results_WBZ_adapt');
path_out_p1 = fullfile(path_out, 'Figures');

if ~exist("path_out", 'dir')
    mkdir(path_out)
end

if ~exist("path_out_p1", 'dir')
    mkdir(path_out_p1)
end

% THICKNESS

% Average thickness
[p,h] = signrank(ThicknessResults.FWZ.normal.OA,ThicknessResults.FNWZ.normal.OA);
[p2,h2] = signrank(ThicknessResults.FWZ.normal.CO,ThicknessResults.FNWZ.normal.CO);

% peak thickness
[p3,h3] = signrank(ThicknessResults.FWZ.peak.OA,ThicknessResults.FNWZ.peak.OA);
[p4,h4] = signrank(ThicknessResults.FWZ.peak.CO,ThicknessResults.FNWZ.peak.CO);

% RELAXATION TIMES

% T1
[p5,h5] = signrank(RTresults.FWBZ.T1.OA,RTresults.FNWBZ.T1.OA);
[p6,h6] = signrank(RTresults.FWBZ.T1.CO,RTresults.FNWBZ.T1.CO);

% T2
[p7,h7] = signrank(RTresults.FWBZ.T2.OA,RTresults.FNWBZ.T2.OA);
[p8,h8] = signrank(RTresults.FWBZ.T2.CO,RTresults.FNWBZ.T2.CO);

% Combine data
pval_matrix_oa = [p;p3;p5;p7];
pval_matrix_co = [p2;p4;p6;p8];

%%%%%%%%% THICKNESS %%%%%%%%%
%%%%%%% OA %%%%%%%
% Averages / std
% FWZ
OA_mean_thick_avg_FWZ = mean(ThicknessResults.FWZ.normal.OA);
OA_mean_thick_peak_FWZ = mean(ThicknessResults.FWZ.peak.OA);

OA_std_thick_avg_FWZ = std(ThicknessResults.FWZ.normal.OA);
OA_std_thick_peak_FWZ = std(ThicknessResults.FWZ.peak.OA);

%FNWZ
OA_mean_thick_avg_FNWZ = mean(ThicknessResults.FNWZ.normal.OA);
OA_mean_thick_peak_FNWZ = mean(ThicknessResults.FNWZ.peak.OA);

OA_std_thick_avg_FNWZ = std(ThicknessResults.FNWZ.normal.OA);
OA_std_thick_peak_FNWZ = std(ThicknessResults.FNWZ.peak.OA);

%%%%%%% CO %%%%%%%
% Averages / std
% FWZ
CO_mean_thick_avg_FWZ = mean(ThicknessResults.FWZ.normal.CO);
CO_mean_thick_peak_FWZ = mean(ThicknessResults.FWZ.peak.CO);

CO_std_thick_avg_FWZ = std(ThicknessResults.FWZ.normal.CO);
CO_std_thick_peak_FWZ = std(ThicknessResults.FWZ.peak.CO);

%FNWZ
CO_mean_thick_avg_FNWZ = mean(ThicknessResults.FNWZ.normal.CO);
CO_mean_thick_peak_FNWZ = mean(ThicknessResults.FNWZ.peak.CO);

CO_std_thick_avg_FNWZ = std(ThicknessResults.FNWZ.normal.CO);
CO_std_thick_peak_FNWZ = std(ThicknessResults.FNWZ.peak.CO);

%%%%%%%%% T1rho relax times %%%%%%%%%
%%%% OA %%%%
OA_mean_t1rho_FWZ = mean(RTresults.FWBZ.T1.OA);
OA_mean_t1rho_FNWZ = mean(RTresults.FNWBZ.T1.OA);

OA_std_t1rho_FWZ = std(RTresults.FWBZ.T1.OA);
OA_std_t1rho_FNWZ = std(RTresults.FNWBZ.T1.OA);

%%%% CO %%%%
CO_mean_t1rho_FWZ = mean(RTresults.FWBZ.T1.CO);
CO_mean_t1rho_FNWZ = mean(RTresults.FNWBZ.T1.CO);

CO_std_t1rho_FWZ = std(RTresults.FWBZ.T1.CO);
CO_std_t1rho_FNWZ = std(RTresults.FNWBZ.T1.CO);

%%%%%%%%% T2 relax times %%%%%%%%%
%%%% OA %%%%
OA_mean_T2_FWZ = mean(RTresults.FWBZ.T2.OA);
OA_mean_T2_FNWZ = mean(RTresults.FNWBZ.T2.OA);

OA_std_T2_FWZ = std(RTresults.FWBZ.T2.OA);
OA_std_T2_FNWZ = std(RTresults.FNWBZ.T2.OA);

%%%% CO %%%%
CO_mean_T2_FWZ = mean(RTresults.FWBZ.T2.CO);
CO_mean_T2_FNWZ = mean(RTresults.FNWBZ.T2.CO);

CO_std_T2_FWZ = std(RTresults.FWBZ.T2.CO);
CO_std_T2_FNWZ = std(RTresults.FNWBZ.T2.CO);

% OA: combine mean and std
WBZ_mean_matrix_oa = [OA_mean_thick_avg_FWZ;OA_mean_thick_peak_FWZ;OA_mean_t1rho_FWZ;OA_mean_T2_FWZ];
WBZ_std_matrix_oa = [OA_std_thick_avg_FWZ;OA_std_thick_peak_FWZ;OA_std_t1rho_FWZ;OA_std_T2_FWZ];
NWBZ_mean_matrix_oa = [OA_mean_thick_avg_FNWZ;OA_mean_thick_peak_FNWZ;OA_mean_t1rho_FNWZ;OA_mean_T2_FNWZ];
NWBZ_std_matrix_oa = [OA_std_thick_avg_FNWZ;OA_std_thick_peak_FNWZ;OA_std_t1rho_FNWZ;OA_std_T2_FNWZ];

rownames = {'Mean Thickness', 'Peak Thickness', 'T1rho relaxation time', 'T2 relaxation time'}';
colnames = {'Variable', 'WBZ mean', 'WBZ std' 'NWBZ mean', 'NWBZ std', 'Pval'};

matrix_oa = [WBZ_mean_matrix_oa WBZ_std_matrix_oa NWBZ_mean_matrix_oa NWBZ_std_matrix_oa pval_matrix_oa];

% Combine data
matrix_total = [rownames num2cell(matrix_oa)];

matrix_total_ns = matrix_total;

for col = 6
    for row = 1:4
        if ~iscellstr(matrix_total_ns(row,col)) && cell2mat(matrix_total_ns(row,col)) > 0.05
            matrix_total_ns(row,col) = {'n.s'};
        end
    end
end

% Convert to table 
sTable_pval = cell2table(matrix_total, 'VariableNames',colnames);
sTable_pval_ns = cell2table(matrix_total_ns, 'VariableNames',colnames);

% Write to Excel file
filename = 'WBZ_adap.xlsx';

writetable(sTable_pval,fullfile(path_out, filename),"Sheet",'OA')
writetable(sTable_pval_ns,fullfile(path_out, filename),"Sheet",'OA_ns')

% CO: combine mean and std
WBZ_mean_matrix_CO = [CO_mean_thick_avg_FWZ;CO_mean_thick_peak_FWZ;CO_mean_t1rho_FWZ;CO_mean_T2_FWZ];
WBZ_std_matrix_CO = [CO_std_thick_avg_FWZ;CO_std_thick_peak_FWZ;CO_std_t1rho_FWZ;CO_std_T2_FWZ];
NWBZ_mean_matrix_CO = [CO_mean_thick_avg_FNWZ;CO_mean_thick_peak_FNWZ;CO_mean_t1rho_FNWZ;CO_mean_T2_FNWZ];
NWBZ_std_matrix_CO = [CO_std_thick_avg_FNWZ;CO_std_thick_peak_FNWZ;CO_std_t1rho_FNWZ;CO_std_T2_FNWZ];

rownames = {'Mean Thickness', 'Peak Thickness', 'T1rho relaxation time', 'T2 relaxation time'}';
colnames = {'Variable', 'WBZ mean', 'WBZ std' 'NWBZ mean', 'NWBZ std', 'Pval'};

matrix_CO = [WBZ_mean_matrix_CO WBZ_std_matrix_CO NWBZ_mean_matrix_CO NWBZ_std_matrix_CO pval_matrix_co];

% Combine data
matrix_total = [rownames num2cell(matrix_CO)];

matrix_total_ns = matrix_total;

for col = 6
    for row = 1:4
        if ~iscellstr(matrix_total_ns(row,col)) && cell2mat(matrix_total_ns(row,col)) > 0.05
            matrix_total_ns(row,col) = {'n.s'};
        end
    end
end

% Convert to table 
sTable_pval = cell2table(matrix_total, 'VariableNames',colnames);
sTable_pval_ns = cell2table(matrix_total_ns, 'VariableNames',colnames);

% Write to Excel file
filename = 'WBZ_adap.xlsx';

writetable(sTable_pval,fullfile(path_out, filename),"Sheet",'CO')
writetable(sTable_pval_ns,fullfile(path_out, filename),"Sheet",'CO_ns')

%% OA Figure

%Mean thickness
y(:,1)=ThicknessResults.FWZ.normal.OA;
y(:,2)=ThicknessResults.FNWZ.normal.OA;

f1 = figure;
    boxchart(y)
    xticklabels({'WBZ', 'NWBZ'})
    ylabel('Mean Thickness [mm]', Interpreter="none")
    title('OA')
    dim = [.73 .01 .3 .3];
    if p<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p,4)) symbol];
    annotation('Textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['MeanThicknessOA' '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

%Peak thickness
y(:,1)=ThicknessResults.FWZ.peak.OA;
y(:,2)=ThicknessResults.FNWZ.peak.OA;

f1 = figure;
    boxchart(y)
    xticklabels({'WBZ', 'NWBZ'})
    ylabel('Peak Thickness [mm]', Interpreter="none")
    title('OA')
    dim = [.73 .01 .3 .3];
    if p3<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p3,4)) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['PeakThicknessOA' '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

% T1 relax time
y(:,1)=RTresults.FWBZ.T1.OA;
y(:,2)=RTresults.FNWBZ.T1.OA;

f1 = figure;
    boxchart(y)
    xticklabels({'WBZ', 'NWBZ'})
    ylabel('T1\rho relaxation time [ms]')
    title('OA')
    dim = [.73 .01 .3 .3];
    if p5<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p5,4)) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['T1_relax_OA' '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

% T2 relax time
y(:,1)=RTresults.FWBZ.T2.OA;
y(:,2)=RTresults.FNWBZ.T2.OA;

f1 = figure;
    boxchart(y)
    xticklabels({'WBZ', 'NWBZ'})
    ylabel('T2 relaxation time [ms]')
    title('OA')
    dim = [.73 .01 .3 .3];
    if p7<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p7,4)) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['T2_relax_OA' '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

%% CO Figure

clear y

%Mean thickness
y(:,1)=ThicknessResults.FWZ.normal.CO;
y(:,2)=ThicknessResults.FNWZ.normal.CO;

f1 = figure;
    boxchart(y)
    xticklabels({'WBZ', 'NWBZ'})
    ylabel('Mean Thickness [mm]', Interpreter="none")
    title('Control')
    dim = [.73 .01 .3 .3];
    if p2<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p2,4)) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['MeanThicknessCO' '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

%Peak thickness
y(:,1)=ThicknessResults.FWZ.peak.CO;
y(:,2)=ThicknessResults.FNWZ.peak.CO;

f1 = figure;
    boxchart(y)
    xticklabels({'WBZ', 'NWBZ'})
    ylabel('Peak Thickness [mm]', Interpreter="none")
    title('Control')
    dim = [.73 .01 .3 .3];
    if p4<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p4,4)) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['PeakThicknessCO' '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

% T1 relax time
y(:,1)=RTresults.FWBZ.T1.CO;
y(:,2)=RTresults.FNWBZ.T1.CO;

f1 = figure;
    boxchart(y)
    xticklabels({'WBZ', 'NWBZ'})
    ylabel('T1\rho relaxation time [ms]')
    title('Control')
    dim = [.73 .01 .3 .3];
    if p6<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p6,4)) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['T1_relax_CO' '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

% T2 relax time
y(:,1)=RTresults.FWBZ.T2.CO;
y(:,2)=RTresults.FNWBZ.T2.CO;

f1 = figure;
    boxchart(y)
    xticklabels({'WBZ', 'NWBZ'})
    ylabel('T2 relaxation time [ms]')
    title('Control')
    dim = [.73 .01 .3 .3];
    if p8<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p8,4)) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['T2_relax_CO' '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))


