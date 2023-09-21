% Statistics thesis
clear all; 
close all
clc

set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultAxesFontName', 'Arial')

% load in data
load Final_loadvars.mat
load RTresults.mat

side = fieldnames(RTresults);
loadvars_perpeak = fieldnames(loading_results.oa.p1);
loadvars_imp = fieldnames(loading_results.oa.imp);
loadvars_stance = fieldnames(loading_results.oa.stance);

% Results directories
path_out = fullfile(cd,'Results_load_relaxT2');
path_out_p1 = fullfile(path_out, 'p1');
path_out_p2 = fullfile(path_out, 'p2');
path_out_imp = fullfile(path_out, 'imp');
path_out_stance = fullfile(path_out, 'stance');

if ~exist("path_out", 'dir')
    mkdir(path_out)
end

if ~exist("path_out_p1", 'dir')
    mkdir(path_out_p1)
end

if ~exist("path_out_p2", 'dir')
    mkdir(path_out_p2)
end

if ~exist("path_out_imp", 'dir')
    mkdir(path_out_imp)
end

if ~exist("path_out_stance", 'dir')
    mkdir(path_out_stance)
end

% Table info
rownames = [loadvars_perpeak; loadvars_perpeak; loadvars_imp; loadvars_stance];
columnames = {'Loadingvars', 'Total', 'Medial', 'Lateral'};

% Loading variables x Thickness(mean, peak) (lat/med) (1 group = OA)

for y = [2,4,5]
    for q = 1:length(loadvars_perpeak) % loading variables per peak
        if y==2
            c=1;
        elseif y==4
            c=3;
        elseif y==5
            c=2;
        end

        if q<=12
            suffix = '[N]';
        else 
            suffix = '[MPa]';
        end 

        if y==2
            prefix = 'Total';
        elseif y==4
            prefix = 'Lateral';
        elseif y==5
            prefix = 'Medial';
        end

        RTresults.(side{y}).T2.OA
        
        %%%%% OA %%%%%
        % peak 1
        [r,pval] = corr(loading_results.oa.p1.(loadvars_perpeak{q}), RTresults.(side{y}).T2.OA,'Type','Spearman', 'tail','both');
        coeffmatrix_oa_p1(q,c) = r;
        pvalmatrix_oa_p1(q,c) = pval;
        % peak 2
        [r2,pval2] = corr(loading_results.oa.p2.(loadvars_perpeak{q}), RTresults.(side{y}).T2.OA,'Type','Spearman', 'tail','both');
        coeffmatrix_oa_p2(q,c) = r2;
        pvalmatrix_oa_p2(q,c) = pval2;

        %%%%% CO %%%%%
        [r3,pval3] = corr(loading_results.co.p1.(loadvars_perpeak{q}), RTresults.(side{y}).T2.CO,'Type','Spearman', 'tail','both');
        coeffmatrix_p1(q,c) = r3;
        pvalmatrix_p1(q,c) = pval3;
        % peak 2
        [r4,pval4] = corr(loading_results.co.p2.(loadvars_perpeak{q}), RTresults.(side{y}).T2.CO,'Type','Spearman', 'tail','both');
        coeffmatrix_p2(q,c) = r4;
        pvalmatrix_p2(q,c) = pval4;
        
        if q>12
            factor = 1000000;
        else
            factor = 1;
        end

        if contains(loadvars_perpeak{q}, 'med') && contains(side{y}, 'L')
            continue
        end

        if contains(loadvars_perpeak{q}, 'lat') && contains(side{y}, 'M')
            continue
        end
        % Scatter plot with fitted line
        %%% peak 1 %%%
        % OA
        y_mean = mean(loading_results.oa.p1.(loadvars_perpeak{q})/factor);
        x_mean = mean(RTresults.(side{y}).T2.OA);
        y_std = std(loading_results.oa.p1.(loadvars_perpeak{q})/factor);
        x_std = std(RTresults.(side{y}).T2.OA);
        slope = r*y_std/x_std;
        yfit = y_mean + slope*(RTresults.(side{y}).T2.OA - x_mean);
        
        % CO
        y_mean2 = mean(loading_results.co.p1.(loadvars_perpeak{q})/factor);
        x_mean2 = mean(RTresults.(side{y}).T2.CO);
        y_std2 = std(loading_results.co.p1.(loadvars_perpeak{q})/factor);
        x_std2 = std(RTresults.(side{y}).T2.CO);
        slope2 = r3*y_std2/x_std2;
        yfit2 = y_mean2 + slope2*(RTresults.(side{y}).T2.CO - x_mean2);
        
        f1 = figure;
        hold on;
        s = scatter(RTresults.(side{y}).T2.OA, loading_results.oa.p1.(loadvars_perpeak{q})/factor,'rsq','filled');
        p = plot(RTresults.(side{y}).T2.OA, yfit, 'color', 'r', 'LineWidth', 1, 'DisplayName', 'OA');
        s2 = scatter(RTresults.(side{y}).T2.CO, loading_results.co.p1.(loadvars_perpeak{q})/factor,'bsq','filled');
        p2 = plot(RTresults.(side{y}).T2.CO, yfit2, 'color', 'b', 'LineWidth', 1, 'DisplayName','CO');
        
        title('Peak 1')
        xlabel([prefix ' T2\rho relaxation time [ms]'])
        ylabel([(loadvars_perpeak{q}) ' ' suffix], 'Interpreter','none')
        dim = [0.53 .09 .10 .10];
        dim2 = [0.53 0.16 .10 .10];
        str = ['p = ' num2str(pval) ' rho = ' num2str(r)];
        str2 = ['p = ' num2str(pval3) ' rho = ' num2str(r3)];
        a = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor', '#f9f9f9', FaceAlpha=1);
        a2 = annotation('textbox',dim2,'String',str2,'FitBoxToText','on', 'BackgroundColor', '#f9f9f9', FaceAlpha=1);
        a.Color = 'red';
        a2.Color = 'blue';
        hold off;
        
        plotname = (['P1 ' (loadvars_perpeak{q}) ' - ' (side{y}) 'relaxation time' '.png']);
        saveas(gcf, fullfile(path_out_p1,plotname))

        %%% peak 2 %%%
        % OA
        y_mean = mean(loading_results.oa.p2.(loadvars_perpeak{q})/factor);
        x_mean = mean(RTresults.(side{y}).T2.OA);
        y_std = std(loading_results.oa.p2.(loadvars_perpeak{q})/factor);
        x_std = std(RTresults.(side{y}).T2.OA);
        slope = r2*y_std/x_std;
        yfit = y_mean + slope*(RTresults.(side{y}).T2.OA - x_mean);
        
        % CO
        y_mean2 = mean(loading_results.co.p2.(loadvars_perpeak{q})/factor);
        x_mean2 = mean(RTresults.(side{y}).T2.CO);
        y_std2 = std(loading_results.co.p2.(loadvars_perpeak{q})/factor);
        x_std2 = std(RTresults.(side{y}).T2.CO);
        slope2 = r4*y_std2/x_std2;
        yfit2 = y_mean2 + slope2*(RTresults.(side{y}).T2.CO - x_mean2);
        
        f2 = figure;
        hold on;
        s = scatter(RTresults.(side{y}).T2.OA, loading_results.oa.p2.(loadvars_perpeak{q})/factor,'rsq','filled');
        p = plot(RTresults.(side{y}).T2.OA, yfit, 'color', 'r', 'LineWidth', 1, 'DisplayName', 'OA');
        s2 = scatter(RTresults.(side{y}).T2.CO, loading_results.co.p2.(loadvars_perpeak{q})/factor,'bsq','filled');
        p2 = plot(RTresults.(side{y}).T2.CO, yfit2, 'color', 'b', 'LineWidth', 1, 'DisplayName','CO');
        
        title('Peak 2')
        xlabel([prefix ' T2\rho relaxation time [ms]'])
        ylabel([(loadvars_perpeak{q}) ' ' suffix], 'Interpreter','none')
        dim = [0.53 .09 .10 .10];
        dim2 = [0.53 0.16 .10 .10];
        str = ['p = ' num2str(pval2) ' rho = ' num2str(r2)];
        str2 = ['p = ' num2str(pval4) ' rho = ' num2str(r4)];
        a = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor', '#f9f9f9', FaceAlpha=1);
        a2 = annotation('textbox',dim2,'String',str2,'FitBoxToText','on', 'BackgroundColor', '#f9f9f9', FaceAlpha=1);
        a.Color = 'red';
        a2.Color = 'blue';
        hold off;
        
        plotname = (['P2 ' (loadvars_perpeak{q}) ' - ' (side{y}) 'relaxation time' '.png']);
        saveas(gcf, fullfile(path_out_p2,plotname))
        
        close all
    end
end

for y = [2,4,5]
    for q = 1:length(loadvars_imp) % loading variables impulse
        if y==2
            c=1;
        elseif y==4
            c=3;
        elseif y==5
            c=2;
        end

        if y==2
            prefix = 'Total';
        elseif y==4
            prefix = 'Lateral';
        elseif y==5
            prefix = 'Medial';
        end
        %%%%% OA %%%%%
        % Impulse 
        [r5,pval5] = corr(loading_results.oa.imp.(loadvars_imp{q}), RTresults.(side{y}).T2.OA,'Type','Spearman', 'tail','both');
        coeffmatrix_oa_imp(q,c) = r5;
        pvalmatrix_oa_imp(q,c) = pval5;

        %%%%% CO %%%%%
        [r6,pval6] = corr(loading_results.co.imp.(loadvars_imp{q}), RTresults.(side{y}).T2.CO,'Type','Spearman', 'tail','both');
        coeffmatrix_imp(q,c) = r6;
        pvalmatrix_imp(q,c) = pval6;
        
        if contains(loadvars_imp{q}, 'med') && contains(side{y}, 'L')
            continue
        end

        if contains(loadvars_imp{q}, 'lat') && contains(side{y}, 'M')
            continue
        end

        % Scatter plot with fitted line
        % OA
        y_mean = mean(loading_results.oa.imp.(loadvars_imp{q}));
        x_mean = mean(RTresults.(side{y}).T2.OA);
        y_std = std(loading_results.oa.imp.(loadvars_imp{q}));
        x_std = std(RTresults.(side{y}).T2.OA);
        slope = r5*y_std/x_std;
        yfit = y_mean + slope*(RTresults.(side{y}).T2.OA - x_mean);
        
        % CO
        y_mean2 = mean(loading_results.co.imp.(loadvars_imp{q}));
        x_mean2 = mean(RTresults.(side{y}).T2.CO);
        y_std2 = std(loading_results.co.imp.(loadvars_imp{q}));
        x_std2 = std(RTresults.(side{y}).T2.CO);
        slope2 = r6*y_std2/x_std2;
        yfit2 = y_mean2 + slope2*(RTresults.(side{y}).T2.CO - x_mean2);
        
        f1 = figure;
        hold on;
        s = scatter(RTresults.(side{y}).T2.OA, loading_results.oa.imp.(loadvars_imp{q}),'rsq','filled');
        p = plot(RTresults.(side{y}).T2.OA, yfit, 'color', 'r', 'LineWidth', 1, 'DisplayName', 'OA');
        s2 = scatter(RTresults.(side{y}).T2.CO, loading_results.co.imp.(loadvars_imp{q}),'bsq','filled');
        p2 = plot(RTresults.(side{y}).T2.CO, yfit2, 'color', 'b', 'LineWidth', 1, 'DisplayName','CO');
        
        title('Stance Phase')
        xlabel([prefix ' T2\rho relaxation time [ms]'])
        ylabel(['Impulse ' (loadvars_imp{q}) ' [N*s]'], 'Interpreter','none')
        dim = [0.53 .80 .10 .10];
        dim2 = [0.53 0.72 .10 .10];
        str = ['p = ' num2str(pval5) ' rho = ' num2str(r5)];
        str2 = ['p = ' num2str(pval6) ' rho = ' num2str(r6)];
        a = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor', '#f9f9f9', FaceAlpha=1);
        a2 = annotation('textbox',dim2,'String',str2,'FitBoxToText','on', 'BackgroundColor', '#f9f9f9', FaceAlpha=1);
        a.Color = 'red';
        a2.Color = 'blue';
        hold off;
        
        plotname = (['Impulse ' (loadvars_imp{q}) ' - ' (side{y}) 'relaxation time' '.png']);
        saveas(gcf, fullfile(path_out_imp,plotname))
        close all
    end
end


for y = [2,4,5]   
    for q = 1:length(loadvars_stance) % loading variables stance
        if y==2
            c=1;
        elseif y==4
            c=3;
        elseif y==5
            c=2;
        end

        if y==2
            prefix = 'Total';
        elseif y==4
            prefix = 'Lateral';
        elseif y==5
            prefix = 'Medial';
        end

        if contains(loadvars_stance{q}, 'med') && contains(side{y}, 'L')
            continue
        end

        if contains(loadvars_stance{q}, 'lat') && contains(side{y}, 'M')
            continue
        end
        %%%% OA %%%%%
        [r7,pval7] = corr(loading_results.oa.stance.(loadvars_stance{q}), RTresults.(side{y}).T2.OA,'Type','Spearman', 'tail','both');
        coeffmatrix_oa_stance(q,c) = r7;
        pvalmatrix_oa_stance(q,c) = pval7;

        %%%% CO %%%%%            
        [r8,pval8] = corr(loading_results.co.stance.(loadvars_stance{q}), RTresults.(side{y}).T2.CO,'Type','Spearman', 'tail','both');
        coeffmatrix_stance(q,c) = r8;
        pvalmatrix_stance(q,c) = pval8;

        % Scatter plot with fitted line
        % OA
        y_mean = mean(loading_results.oa.stance.(loadvars_stance{q})/1000000);
        x_mean = mean(RTresults.(side{y}).T2.OA);
        y_std = std(loading_results.oa.stance.(loadvars_stance{q})/1000000);
        x_std = std(RTresults.(side{y}).T2.OA);
        slope = r7*y_std/x_std;
        yfit = y_mean + slope*(RTresults.(side{y}).T2.OA - x_mean);
        
        % CO
        y_mean2 = mean(loading_results.co.stance.(loadvars_stance{q})/1000000);
        x_mean2 = mean(RTresults.(side{y}).T2.CO);
        y_std2 = std(loading_results.co.stance.(loadvars_stance{q})/1000000);
        x_std2 = std(RTresults.(side{y}).T2.CO);
        slope2 = r8*y_std2/x_std2;
        yfit2 = y_mean2 + slope2*(RTresults.(side{y}).T2.CO - x_mean2);
        
        f1 = figure;
        hold on;
        s = scatter(RTresults.(side{y}).T2.OA, loading_results.oa.stance.(loadvars_stance{q})/1000000,'rsq','filled');
        p = plot(RTresults.(side{y}).T2.OA, yfit, 'color', 'r', 'LineWidth', 1, 'DisplayName', 'OA');
        s2 = scatter(RTresults.(side{y}).T2.CO, loading_results.co.stance.(loadvars_stance{q})/1000000,'bsq','filled');
        p2 = plot(RTresults.(side{y}).T2.CO, yfit2, 'color', 'b', 'LineWidth', 1, 'DisplayName','CO');
        
        title('Stance Phase')
        xlabel([prefix ' T2\rho relaxation time [ms]'])
        ylabel(['stance ' (loadvars_stance{q}) ' [N*s]'], 'Interpreter','none')
        dim = [0.53 .09 .10 .10];
        dim2 = [0.53 0.16 .10 .10];
        str = ['p = ' num2str(pval7) ' rho = ' num2str(r7)];
        str2 = ['p = ' num2str(pval8) ' rho = ' num2str(r8)];
        a = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor', '#f9f9f9', FaceAlpha=1);
        a2 = annotation('textbox',dim2,'String',str2,'FitBoxToText','on', 'BackgroundColor', '#f9f9f9', FaceAlpha=1);
        a.Color = 'red';
        a2.Color = 'blue';
        hold off;
       
        plotname = (['stance ' (loadvars_stance{q}) ' - ' (side{y}) 'relaxation time' '.png']);
        saveas(gcf, fullfile(path_out_stance,plotname))
        close all
    end
end


%%%% OA %%%%%
% Combine data
pvalmatrix_temp_oa = [pvalmatrix_oa_p1; pvalmatrix_oa_p2; pvalmatrix_oa_imp; pvalmatrix_oa_stance];
pvalmatrix_total_oa = [rownames num2cell(pvalmatrix_temp_oa)];

coeffmatrix_temp_oa = [coeffmatrix_oa_p1; coeffmatrix_oa_p2; coeffmatrix_oa_imp; coeffmatrix_oa_stance];
coeffmatrix_total_oa = [rownames num2cell(coeffmatrix_temp_oa)];

% Convert to table 
sTable_pval_oa = cell2table(pvalmatrix_total_oa, 'VariableNames',columnames);
sTable_coeff_oa = cell2table(coeffmatrix_total_oa, 'VariableNames',columnames);

%%%% CO %%%%%
% Combine data
pvalmatrix_temp = [pvalmatrix_p1; pvalmatrix_p2; pvalmatrix_imp; pvalmatrix_stance];
pvalmatrix_total = [rownames num2cell(pvalmatrix_temp)];

coeffmatrix_temp = [coeffmatrix_p1; coeffmatrix_p2; coeffmatrix_imp; coeffmatrix_stance];
coeffmatrix_total = [rownames num2cell(coeffmatrix_temp)];

% Filter out med-lat, lat-med correlations
for rows = 1:length(rownames)
    for columns = 2:length(columnames)
        if contains(rownames{rows}, 'med') && contains(columnames{columns}, 'Lateral')
            pvalmatrix_total(rows,columns) = {'n.a.'};
            pvalmatrix_total_oa(rows,columns) = {'n.a.'};
            coeffmatrix_total(rows,columns) = {'n.a.'};
            coeffmatrix_total_oa(rows,columns) = {'n.a.'};
        elseif contains(rownames{rows}, 'lat') && contains(columnames{columns}, 'Medial')
            pvalmatrix_total(rows,columns) = {'n.a.'};
            pvalmatrix_total_oa(rows,columns) = {'n.a.'};
            coeffmatrix_total(rows,columns) = {'n.a.'};
            coeffmatrix_total_oa(rows,columns) = {'n.a.'};
        end
    end
end


pvalmatrix_total_ns = pvalmatrix_total;
pvalmatrix_total_oa_ns = pvalmatrix_total_oa;

for col = 2:4
    for row = 1:length(pvalmatrix_total_ns)
        if ~iscellstr(pvalmatrix_total_ns(row,col)) && cell2mat(pvalmatrix_total_ns(row,col)) > 0.05
            pvalmatrix_total_ns(row,col) = {'n.s'};
        end
    end
end

for col = 2:4
    for row = 1:length(pvalmatrix_total_oa_ns)
        if ~iscellstr(pvalmatrix_total_oa_ns(row,col)) && cell2mat(pvalmatrix_total_oa_ns(row,col)) > 0.05
            pvalmatrix_total_oa_ns(row,col) = {'n.s'};
        end
    end
end

% Convert to table 
sTable_pval = cell2table(pvalmatrix_total, 'VariableNames',columnames);
sTable_pval_ns = cell2table(pvalmatrix_total_ns, 'VariableNames',columnames);
sTable_coeff = cell2table(coeffmatrix_total, 'VariableNames',columnames);
sTable_pval_oa = cell2table(pvalmatrix_total_oa, 'VariableNames',columnames);
sTable_pval_oa_ns = cell2table(pvalmatrix_total_oa_ns, 'VariableNames',columnames);
sTable_coeff_oa = cell2table(coeffmatrix_total_oa, 'VariableNames',columnames);
                
filename = 'Loading_x_relaxT2.xlsx';

writetable(sTable_pval,fullfile(path_out, filename),"Sheet",'pval_controls')
writetable(sTable_pval_ns,fullfile(path_out, filename),"Sheet",'pval_controls_ns')
writetable(sTable_coeff,fullfile(path_out,filename),"Sheet",'coeff_controls')
writetable(sTable_pval_oa,fullfile(path_out,filename),"Sheet",'pval_OA')
writetable(sTable_pval_oa_ns,fullfile(path_out,filename),"Sheet",'pval_OA_ns')
writetable(sTable_coeff_oa,fullfile(path_out,filename),"Sheet",'coeff_OA')
  

          

            









