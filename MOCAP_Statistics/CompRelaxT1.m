% Compare relax times T1

clear all;
close all;
clc

set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultAxesFontName', 'Arial')

load RTresults.mat

side = fieldnames(RTresults);

% Results directories
path_out = fullfile(cd,'Results_relax_T1_Comp');
path_out_p1 = fullfile(path_out, 'Figures');

if ~exist("path_out", 'dir')
    mkdir(path_out)
end

if ~exist("path_out_p1", 'dir')
    mkdir(path_out_p1)
end

% Table info
rownames = side;
columnames = {'Part', 'OA mean','OA std', 'Control mean', 'Control std', 'Pval'};

c = -1;
c2 = 0;
hh = [];
h2 = [];
for q = 1:length(side) 
    % Peak 1
    mean_oa = mean(RTresults.(side{q}).T1.OA);
    std_oa = std(RTresults.(side{q}).T1.OA);
    mean_co = mean(RTresults.(side{q}).T1.CO);
    std_co = std(RTresults.(side{q}).T1.CO);

    pval_matrix_p1(q,1) = mean_oa;
    pval_matrix_p1(q,2) = std_oa;
    pval_matrix_p1(q,3) = mean_co;
    pval_matrix_p1(q,4) = std_co;

    [p,h] = ranksum(RTresults.(side{q}).T1.OA,RTresults.(side{q}).T1.CO);

    pval_matrix_p1(q,5) = p;
    
    % Iterators 
    c = c+2;
    c2 = c2+2;
    
    % Combine data
    y(:,c) = RTresults.(side{q}).T1.OA;
    y_temp = RTresults.(side{q}).T1.CO;
    y_temp(6:10,1) = NaN;
    y(:,c2) = y_temp;   

    hh = [hh; RTresults.(side{q}).T1.OA];
    h2 = [h2; y_temp];
end
h3 = [hh;h2];

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

Finalmatrix = [group part value];
colnames = {'group', 'part', 'value'};
table = array2table(Finalmatrix, 'VariableNames',colnames);


% Boxplot
boxchart(table.part,table.value,'GroupByColor',table.group)
ylabel('T1\rho relaxation time [ms]')
xticks([1 2 3 4 5])
xticklabels({'FULLFEM','FWBZ', 'FNWBZ', 'LATWBZ', 'MEDWBZ'})
legend('Control', 'OA', Location='southeast')

plotname = (['Relaxation times T1' '.png']);
saveas(gcf, fullfile(path_out_p1,plotname))

close gcf

% Combine data
pvalmatrix_temp = pval_matrix_p1;
pvalmatrix_total = [rownames num2cell(pvalmatrix_temp)];

pvalmatrix_total_ns = pvalmatrix_total;

for col = 6
    for row = 1:5
        if ~iscellstr(pvalmatrix_total_ns(row,col)) && cell2mat(pvalmatrix_total_ns(row,col)) > 0.05
            pvalmatrix_total_ns(row,col) = {'n.s'};
        end
    end
end

% Convert to table 
sTable_pval = cell2table(pvalmatrix_total, 'VariableNames',columnames);
sTable_pval_ns = cell2table(pvalmatrix_total_ns, 'VariableNames',columnames);

% Write to Excel file
filename = 'CompRelaxT1.xlsx';

writetable(sTable_pval,fullfile(path_out, filename),"Sheet",'OA_vs_Control')
writetable(sTable_pval_ns,fullfile(path_out, filename),"Sheet",'OA_vs_Control_ns')